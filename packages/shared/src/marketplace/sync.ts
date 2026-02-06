/**
 * Marketplace Metadata Sync
 *
 * 管理应用市场元数据的同步
 * - 定期同步（4小时间隔）
 * - 只同步元数据，不下载文件
 * - 支持后台更新
 */

import { existsSync, readFileSync, writeFileSync, mkdirSync } from 'fs';
import { join, dirname } from 'path';
import { getMarketplaceCacheDir, saveMarketplaceCache } from './storage.ts';
import { listSkills, listApps, listCategories } from './api.ts';
import { debug } from '../utils/debug.ts';
import type { MarketplaceSkill, MarketplaceApp } from './types.ts';

// ============================================================
// 常量
// ============================================================

const SYNC_INTERVAL_MS = 4 * 60 * 60 * 1000; // 4 小时
const LAST_SYNC_FILE = 'last-sync.json';

// ============================================================
// 类型定义
// ============================================================

interface LastSyncData {
  lastSync: string;
  nextSync: string;
}

// ============================================================
// 同步时间追踪
// ============================================================

/**
 * 获取上次同步时间文件路径
 */
function getLastSyncFilePath(): string {
  return join(getMarketplaceCacheDir(), LAST_SYNC_FILE);
}

/**
 * 读取上次同步时间
 */
function readLastSyncTime(): Date | null {
  const filePath = getLastSyncFilePath();
  if (!existsSync(filePath)) {
    return null;
  }

  try {
    const content = readFileSync(filePath, 'utf-8');
    const data = JSON.parse(content) as LastSyncData;
    return new Date(data.lastSync);
  } catch (error) {
    debug('[readLastSyncTime] 读取同步时间失败:', error);
    return null;
  }
}

/**
 * 写入同步时间
 */
function writeLastSyncTime(): void {
  const cacheDir = getMarketplaceCacheDir();
  if (!existsSync(cacheDir)) {
    mkdirSync(cacheDir, { recursive: true });
  }

  const now = new Date();
  const nextSync = new Date(now.getTime() + SYNC_INTERVAL_MS);

  const data: LastSyncData = {
    lastSync: now.toISOString(),
    nextSync: nextSync.toISOString(),
  };

  const filePath = getLastSyncFilePath();
  writeFileSync(filePath, JSON.stringify(data, null, 2));
}

/**
 * 检查是否需要同步
 */
function shouldSync(): boolean {
  const lastSync = readLastSyncTime();
  if (!lastSync) {
    return true; // 从未同步过
  }

  const timeSinceLastSync = Date.now() - lastSync.getTime();
  return timeSinceLastSync >= SYNC_INTERVAL_MS;
}

// ============================================================
// 同步操作
// ============================================================

/**
 * 同步应用市场元数据
 *
 * 从云端获取技能、应用和分类信息并保存到本地缓存
 * 如果距离上次同步不足 4 小时，则跳过
 */
export async function syncMarketplaceMetadata(): Promise<void> {
  // 检查是否需要同步
  if (!shouldSync()) {
    debug('[syncMarketplaceMetadata] 距离上次同步不足 4 小时，跳过');
    return;
  }

  debug('[syncMarketplaceMetadata] 开始同步应用市场元数据...');

  try {
    // 分页获取所有数据（API 限制每页最多 200 条）
    const PAGE_SIZE = 200;

    // 获取技能列表（分页）
    const allSkills: MarketplaceSkill[] = [];
    let skillPage = 1;
    let hasMoreSkills = true;

    while (hasMoreSkills) {
      const skillsResponse = await listSkills({ page: skillPage, size: PAGE_SIZE });
      allSkills.push(...skillsResponse.items);

      // 计算总页数
      const totalPages = Math.ceil(skillsResponse.total / PAGE_SIZE);
      hasMoreSkills = skillPage < totalPages;
      skillPage++;

      debug(`[syncMarketplaceMetadata] 已获取技能: ${allSkills.length}/${skillsResponse.total}`);
    }

    // 获取应用列表（分页）
    const allApps: MarketplaceApp[] = [];
    let appPage = 1;
    let hasMoreApps = true;

    while (hasMoreApps) {
      const appsResponse = await listApps({ page: appPage, size: PAGE_SIZE });
      allApps.push(...appsResponse.items);

      // 计算总页数
      const totalPages = Math.ceil(appsResponse.total / PAGE_SIZE);
      hasMoreApps = appPage < totalPages;
      appPage++;

      debug(`[syncMarketplaceMetadata] 已获取应用: ${allApps.length}/${appsResponse.total}`);
    }

    // 获取分类列表（通常不多，一次获取）
    const categoriesResponse = await listCategories();

    // 验证数据有效性
    if (!allSkills || !allApps || !categoriesResponse.items) {
      throw new Error('同步响应数据格式无效');
    }

    // 只有成功获取数据才保存到缓存
    saveMarketplaceCache({
      skills: allSkills,
      apps: allApps,
      categories: categoriesResponse.items,
    });

    // 更新同步时间（只在成功时更新）
    writeLastSyncTime();

    debug(
      `[syncMarketplaceMetadata] 同步完成: ${allSkills.length} 个技能, ` +
      `${allApps.length} 个应用, ${categoriesResponse.items.length} 个分类`
    );
  } catch (error) {
    debug('[syncMarketplaceMetadata] 同步失败:', error);
    // 不更新同步时间，允许下次启动重试
    throw error;
  }
}

/**
 * 强制同步（忽略时间间隔）
 */
export async function forceSyncMarketplaceMetadata(): Promise<void> {
  debug('[forceSyncMarketplaceMetadata] 强制同步应用市场元数据...');

  try {
    // 分页获取所有数据（API 限制每页最多 200 条）
    const PAGE_SIZE = 200;

    // 获取技能列表（分页）
    const allSkills: MarketplaceSkill[] = [];
    let skillPage = 1;
    let hasMoreSkills = true;

    while (hasMoreSkills) {
      const skillsResponse = await listSkills({ page: skillPage, size: PAGE_SIZE });
      allSkills.push(...skillsResponse.items);

      // 计算总页数
      const totalPages = Math.ceil(skillsResponse.total / PAGE_SIZE);
      hasMoreSkills = skillPage < totalPages;
      skillPage++;
    }

    // 获取应用列表（分页）
    const allApps: MarketplaceApp[] = [];
    let appPage = 1;
    let hasMoreApps = true;

    while (hasMoreApps) {
      const appsResponse = await listApps({ page: appPage, size: PAGE_SIZE });
      allApps.push(...appsResponse.items);

      // 计算总页数
      const totalPages = Math.ceil(appsResponse.total / PAGE_SIZE);
      hasMoreApps = appPage < totalPages;
      appPage++;
    }

    // 获取分类列表
    const categoriesResponse = await listCategories();

    // 保存到缓存
    saveMarketplaceCache({
      skills: allSkills,
      apps: allApps,
      categories: categoriesResponse.items,
    });

    // 更新同步时间
    writeLastSyncTime();

    debug(
      `[forceSyncMarketplaceMetadata] 同步完成: ${allSkills.length} 个技能, ` +
      `${allApps.length} 个应用, ${categoriesResponse.items.length} 个分类`
    );
  } catch (error) {
    debug('[forceSyncMarketplaceMetadata] 同步失败:', error);
    throw error;
  }
}

/**
 * 获取下次同步时间
 */
export function getNextSyncTime(): Date | null {
  const filePath = getLastSyncFilePath();
  if (!existsSync(filePath)) {
    return null;
  }

  try {
    const content = readFileSync(filePath, 'utf-8');
    const data = JSON.parse(content) as LastSyncData;
    return new Date(data.nextSync);
  } catch (error) {
    debug('[getNextSyncTime] 读取下次同步时间失败:', error);
    return null;
  }
}

/**
 * 获取上次同步时间
 */
export function getLastSyncTime(): Date | null {
  return readLastSyncTime();
}
