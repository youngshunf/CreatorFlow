/**
 * Global Plugins System
 *
 * 管理全局插件（存储在 ~/.sprouty-ai/plugins/）
 * 全局插件在首次启动时从内置插件复制，并自动加载到所有工作区
 *
 * 插件结构：{name}/.claude-plugin/plugin.json + commands/ + skills/ 等
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  readdirSync,
  copyFileSync,
} from 'fs';
import { join, dirname } from 'path';
import { homedir } from 'os';
import { fileURLToPath } from 'url';
import { debug } from '../utils/debug.ts';

// ============================================================
// 路径工具
// ============================================================

/**
 * 获取全局插件目录路径
 */
export function getGlobalPluginsDir(): string {
  return join(homedir(), '.sprouty-ai', 'plugins');
}

/**
 * 获取全局插件元数据文件路径
 */
function getGlobalPluginsMetaPath(): string {
  return join(getGlobalPluginsDir(), '.meta.json');
}

/**
 * 获取内置插件源目录路径
 * 支持开发环境和打包环境
 */
function getBundledPluginsDir(): string {
  const possiblePaths: string[] = [];

  // 1. ESM 上下文（开发环境 - 源代码目录）
  try {
    const currentFile = fileURLToPath(import.meta.url);
    possiblePaths.push(join(dirname(currentFile), '..', 'apps', 'bundled-plugins'));
  } catch {
    // ESM 不可用
  }

  // 2. CJS 上下文（开发环境 - 源代码目录）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, '..', 'apps', 'bundled-plugins'));
  }

  // 3. dist/assets/bundled-plugins（Electron 开发模式）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'assets', 'bundled-plugins'));
    possiblePaths.push(join(__dirname, '..', 'assets', 'bundled-plugins'));
  }

  // 4. 打包的 Electron 应用路径（resources/bundled-plugins）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'resources', 'bundled-plugins'));
    possiblePaths.push(join(__dirname, '..', 'resources', 'bundled-plugins'));
  }

  // 5. process.resourcesPath（Electron 特定）
  if (typeof process !== 'undefined' && (process as any).resourcesPath) {
    possiblePaths.push(join((process as any).resourcesPath, 'bundled-plugins'));
  }

  // 返回第一个存在的路径
  for (const p of possiblePaths) {
    if (existsSync(p)) {
      debug(`[getBundledPluginsDir] 找到内置插件目录: ${p}`);
      return p;
    }
  }

  // 如果所有路径都不存在，抛出错误
  const fallbackPath = possiblePaths[0];
  if (!fallbackPath) {
    throw new Error('无法确定内置插件目录路径，请检查打包配置');
  }

  debug(`[getBundledPluginsDir] 警告：使用降级路径 ${fallbackPath}，但该路径不存在`);
  return fallbackPath;
}

// ============================================================
// 元数据类型
// ============================================================

interface GlobalPluginMeta {
  version: string;
  source: 'bundled' | 'marketplace';
  installedAt: string;
}

interface GlobalPluginsMetadata {
  version: string;
  lastSync: string;
  plugins: Record<string, GlobalPluginMeta>;
}

// ============================================================
// 元数据操作
// ============================================================

/**
 * 读取全局插件元数据
 */
function readGlobalPluginsMeta(): GlobalPluginsMetadata | null {
  const metaPath = getGlobalPluginsMetaPath();
  if (!existsSync(metaPath)) {
    return null;
  }

  try {
    const content = readFileSync(metaPath, 'utf-8');
    return JSON.parse(content) as GlobalPluginsMetadata;
  } catch (error) {
    debug('[readGlobalPluginsMeta] 读取元数据失败:', error);
    return null;
  }
}

/**
 * 写入全局插件元数据
 */
function writeGlobalPluginsMeta(meta: GlobalPluginsMetadata): void {
  const metaPath = getGlobalPluginsMetaPath();
  const dir = dirname(metaPath);

  if (!existsSync(dir)) {
    mkdirSync(dir, { recursive: true });
  }

  writeFileSync(metaPath, JSON.stringify(meta, null, 2));
}

// ============================================================
// 文件复制工具
// ============================================================

/**
 * 递归复制目录
 */
function copyDirectoryRecursive(source: string, target: string): void {
  if (!existsSync(target)) {
    mkdirSync(target, { recursive: true });
  }

  const entries = readdirSync(source, { withFileTypes: true });

  for (const entry of entries) {
    const sourcePath = join(source, entry.name);
    const targetPath = join(target, entry.name);

    if (entry.isDirectory()) {
      copyDirectoryRecursive(sourcePath, targetPath);
    } else {
      copyFileSync(sourcePath, targetPath);
    }
  }
}

// ============================================================
// 全局插件初始化
// ============================================================

/**
 * 初始化全局插件系统
 * 首次启动时将内置插件复制到全局插件目录 ~/.sprouty-ai/plugins/
 */
export async function initializeGlobalPlugins(): Promise<void> {
  const globalPluginsDir = getGlobalPluginsDir();
  const metaFile = getGlobalPluginsMetaPath();

  // 检查是否已初始化
  if (existsSync(metaFile)) {
    debug('[initializeGlobalPlugins] 全局插件已初始化，跳过');
    return;
  }

  debug('[initializeGlobalPlugins] 开始初始化全局插件...');

  // 获取内置插件目录
  const bundledPluginsDir = getBundledPluginsDir();
  if (!bundledPluginsDir || !existsSync(bundledPluginsDir)) {
    debug('[initializeGlobalPlugins] 未找到内置插件目录，跳过初始化');
    return;
  }

  // 创建全局插件目录
  if (!existsSync(globalPluginsDir)) {
    mkdirSync(globalPluginsDir, { recursive: true });
  }

  // 初始化元数据
  const meta: GlobalPluginsMetadata = {
    version: '1.0.0',
    lastSync: new Date().toISOString(),
    plugins: {},
  };

  // 复制所有内置插件
  try {
    const entries = readdirSync(bundledPluginsDir, { withFileTypes: true });
    let copiedCount = 0;

    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const sourceDir = join(bundledPluginsDir, entry.name);
      const targetDir = join(globalPluginsDir, entry.name);

      // 检查是否有 .claude-plugin/plugin.json
      const pluginFile = join(sourceDir, '.claude-plugin', 'plugin.json');
      if (!existsSync(pluginFile)) {
        debug(`[initializeGlobalPlugins] 跳过 ${entry.name}：缺少 .claude-plugin/plugin.json`);
        continue;
      }

      // 如果目标目录已存在，跳过（避免覆盖用户修改）
      if (existsSync(targetDir)) {
        debug(`[initializeGlobalPlugins] 跳过 ${entry.name}：目录已存在`);
        continue;
      }

      // 读取插件版本
      let pluginVersion = '1.0.0';
      try {
        const pluginJson = JSON.parse(readFileSync(pluginFile, 'utf-8'));
        pluginVersion = pluginJson.version || '1.0.0';
      } catch {
        // 使用默认版本
      }

      // 复制插件目录
      try {
        copyDirectoryRecursive(sourceDir, targetDir);

        // 记录到元数据
        meta.plugins[entry.name] = {
          version: pluginVersion,
          source: 'bundled',
          installedAt: new Date().toISOString(),
        };

        copiedCount++;
        debug(`[initializeGlobalPlugins] 已复制插件: ${entry.name}`);
      } catch (error) {
        debug(`[initializeGlobalPlugins] 复制插件 ${entry.name} 失败:`, error);
      }
    }

    // 保存元数据
    writeGlobalPluginsMeta(meta);
    debug(`[initializeGlobalPlugins] 初始化完成，共复制 ${copiedCount} 个插件`);
  } catch (error) {
    debug('[initializeGlobalPlugins] 初始化失败:', error);
    throw error;
  }
}

// ============================================================
// 工具函数
// ============================================================

/**
 * 检查全局插件是否已初始化
 */
export function isGlobalPluginsInitialized(): boolean {
  return existsSync(getGlobalPluginsMetaPath());
}

/**
 * 获取全局插件元数据
 */
export function getGlobalPluginsMeta(): GlobalPluginsMetadata | null {
  return readGlobalPluginsMeta();
}
