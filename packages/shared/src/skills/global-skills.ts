/**
 * Global Skills System
 *
 * 管理全局技能（存储在 ~/.creator-flow/skills/）
 * 全局技能在首次启动时从内置技能复制，并自动加载到所有工作区
 *
 * 兼容 Claude 技能系统：使用相同的目录结构 ~/.creator-flow/skills/
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  readdirSync,
  copyFileSync,
  statSync,
} from 'fs';
import { join, dirname } from 'path';
import { homedir } from 'os';
import { fileURLToPath } from 'url';
import type { LoadedSkill } from './types.ts';
import { debug } from '../utils/debug.ts';

// ============================================================
// 路径工具
// ============================================================

/**
 * 获取全局技能目录路径
 * 使用 ~/.creator-flow/skills/ 以兼容 Claude 技能系统
 */
export function getGlobalSkillsDir(): string {
  return join(homedir(), '.creator-flow', 'skills');
}

/**
 * 获取全局技能元数据文件路径
 */
function getGlobalSkillsMetaPath(): string {
  return join(getGlobalSkillsDir(), '.meta.json');
}

/**
 * 获取内置技能源目录路径
 * 支持开发环境和打包环境
 */
function getBundledSkillsDir(): string {
  const possiblePaths: string[] = [];

  // 1. ESM 上下文（开发环境 - 源代码目录）
  try {
    const currentFile = fileURLToPath(import.meta.url);
    possiblePaths.push(join(dirname(currentFile), '..', 'apps', 'bundled-skills'));
  } catch {
    // ESM 不可用
  }

  // 2. CJS 上下文（开发环境 - 源代码目录）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, '..', 'apps', 'bundled-skills'));
  }

  // 3. dist/assets/bundled-skills（Electron 开发模式）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'assets', 'bundled-skills'));
    possiblePaths.push(join(__dirname, '..', 'assets', 'bundled-skills'));
  }

  // 4. 打包的 Electron 应用路径（resources/bundled-skills）
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'resources', 'bundled-skills'));
    possiblePaths.push(join(__dirname, '..', 'resources', 'bundled-skills'));
  }

  // 5. process.resourcesPath（Electron 特定）
  if (typeof process !== 'undefined' && (process as any).resourcesPath) {
    possiblePaths.push(join((process as any).resourcesPath, 'bundled-skills'));
  }

  // 返回第一个存在的路径
  for (const p of possiblePaths) {
    if (existsSync(p)) {
      debug(`[getBundledSkillsDir] 找到内置技能目录: ${p}`);
      return p;
    }
  }

  // 如果所有路径都不存在，抛出错误
  const fallbackPath = possiblePaths[0];
  if (!fallbackPath) {
    throw new Error('无法确定内置技能目录路径，请检查打包配置');
  }

  debug(`[getBundledSkillsDir] 警告：使用降级路径 ${fallbackPath}，但该路径不存在`);
  return fallbackPath;
}

// ============================================================
// 元数据类型
// ============================================================

interface GlobalSkillMeta {
  version: string;
  source: 'bundled' | 'marketplace';
  installedAt: string;
}

interface GlobalSkillsMetadata {
  version: string;
  lastSync: string;
  skills: Record<string, GlobalSkillMeta>;
}

// ============================================================
// 元数据操作
// ============================================================

/**
 * 读取全局技能元数据
 */
function readGlobalSkillsMeta(): GlobalSkillsMetadata | null {
  const metaPath = getGlobalSkillsMetaPath();
  if (!existsSync(metaPath)) {
    return null;
  }

  try {
    const content = readFileSync(metaPath, 'utf-8');
    return JSON.parse(content) as GlobalSkillsMetadata;
  } catch (error) {
    debug('[readGlobalSkillsMeta] 读取元数据失败:', error);
    return null;
  }
}

/**
 * 写入全局技能元数据
 */
function writeGlobalSkillsMeta(meta: GlobalSkillsMetadata): void {
  const metaPath = getGlobalSkillsMetaPath();
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
// 全局技能初始化
// ============================================================

/**
 * 初始化全局技能系统
 * 首次启动时将内置技能复制到全局技能目录 ~/.creator-flow/skills/
 */
export async function initializeGlobalSkills(): Promise<void> {
  const globalSkillsDir = getGlobalSkillsDir();
  const metaFile = getGlobalSkillsMetaPath();

  // 检查是否已初始化
  if (existsSync(metaFile)) {
    debug('[initializeGlobalSkills] 全局技能已初始化，跳过');
    return;
  }

  debug('[initializeGlobalSkills] 开始初始化全局技能...');

  // 获取内置技能目录
  const bundledSkillsDir = getBundledSkillsDir();
  if (!bundledSkillsDir || !existsSync(bundledSkillsDir)) {
    debug('[initializeGlobalSkills] 未找到内置技能目录，跳过初始化');
    return;
  }

  // 创建全局技能目录
  if (!existsSync(globalSkillsDir)) {
    mkdirSync(globalSkillsDir, { recursive: true });
  }

  // 初始化元数据
  const meta: GlobalSkillsMetadata = {
    version: '1.0.0',
    lastSync: new Date().toISOString(),
    skills: {},
  };

  // 复制所有内置技能
  try {
    const entries = readdirSync(bundledSkillsDir, { withFileTypes: true });
    let copiedCount = 0;

    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const sourceDir = join(bundledSkillsDir, entry.name);
      const targetDir = join(globalSkillsDir, entry.name);

      // 检查是否有 SKILL.md
      const skillFile = join(sourceDir, 'SKILL.md');
      if (!existsSync(skillFile)) {
        debug(`[initializeGlobalSkills] 跳过 ${entry.name}：缺少 SKILL.md`);
        continue;
      }

      // 如果目标目录已存在（可能是用户的 Claude 技能），跳过
      if (existsSync(targetDir)) {
        debug(`[initializeGlobalSkills] 跳过 ${entry.name}：目录已存在（可能是 Claude 技能）`);
        continue;
      }

      // 复制技能目录
      try {
        copyDirectoryRecursive(sourceDir, targetDir);

        // 记录到元数据
        meta.skills[entry.name] = {
          version: '1.0.0',
          source: 'bundled',
          installedAt: new Date().toISOString(),
        };

        copiedCount++;
        debug(`[initializeGlobalSkills] 已复制技能: ${entry.name}`);
      } catch (error) {
        debug(`[initializeGlobalSkills] 复制技能 ${entry.name} 失败:`, error);
      }
    }

    // 保存元数据
    writeGlobalSkillsMeta(meta);
    debug(`[initializeGlobalSkills] 初始化完成，共复制 ${copiedCount} 个技能`);
  } catch (error) {
    debug('[initializeGlobalSkills] 初始化失败:', error);
    throw error;
  }
}

// ============================================================
// 加载全局技能
// ============================================================

/**
 * 从目录加载单个技能
 */
function loadSkillFromDir(skillDir: string, slug: string): LoadedSkill | null {
  const skillFile = join(skillDir, 'SKILL.md');

  // 检查 SKILL.md 是否存在
  if (!existsSync(skillFile)) {
    return null;
  }

  try {
    const content = readFileSync(skillFile, 'utf-8');

    // 解析 frontmatter
    const frontmatterMatch = content.match(/^---\n([\s\S]*?)\n---\n([\s\S]*)$/);
    if (!frontmatterMatch) {
      return null;
    }

    const frontmatter = frontmatterMatch[1];
    const body = frontmatterMatch[2];

    // 提取必需字段
    const nameMatch = frontmatter.match(/name:\s*(.+)/);
    const descMatch = frontmatter.match(/description:\s*(.+)/);

    if (!nameMatch || !descMatch) {
      return null;
    }

    // 提取可选字段
    const iconMatch = frontmatter.match(/icon:\s*(.+)/);
    const globsMatch = frontmatter.match(/globs:\s*\[([\s\S]*?)\]/);
    const alwaysAllowMatch = frontmatter.match(/alwaysAllow:\s*\[([\s\S]*?)\]/);

    const metadata: any = {
      name: nameMatch[1].trim(),
      description: descMatch[1].trim(),
    };

    if (iconMatch) {
      metadata.icon = iconMatch[1].trim();
    }

    if (globsMatch) {
      metadata.globs = globsMatch[1]
        .split(',')
        .map(s => s.trim().replace(/['"]/g, ''))
        .filter(s => s);
    }

    if (alwaysAllowMatch) {
      metadata.alwaysAllow = alwaysAllowMatch[1]
        .split(',')
        .map(s => s.trim().replace(/['"]/g, ''))
        .filter(s => s);
    }

    // 查找图标文件
    let iconPath: string | undefined;
    const iconExtensions = ['.svg', '.png', '.jpg', '.jpeg'];
    for (const ext of iconExtensions) {
      const iconFile = join(skillDir, `icon${ext}`);
      if (existsSync(iconFile)) {
        iconPath = iconFile;
        break;
      }
    }

    return {
      slug,
      metadata,
      content: body,
      iconPath,
      path: skillDir,
    };
  } catch (error) {
    debug(`[loadSkillFromDir] 加载技能 ${slug} 失败:`, error);
    return null;
  }
}

/**
 * 加载所有全局技能
 */
export function loadGlobalSkills(): LoadedSkill[] {
  const globalSkillsDir = getGlobalSkillsDir();

  if (!existsSync(globalSkillsDir)) {
    return [];
  }

  const skills: LoadedSkill[] = [];

  try {
    const entries = readdirSync(globalSkillsDir, { withFileTypes: true });

    for (const entry of entries) {
      // 跳过隐藏文件和非目录
      if (!entry.isDirectory() || entry.name.startsWith('.')) {
        continue;
      }

      const skillDir = join(globalSkillsDir, entry.name);
      const skill = loadSkillFromDir(skillDir, entry.name);

      if (skill) {
        skills.push(skill);
      }
    }
  } catch (error) {
    debug('[loadGlobalSkills] 加载全局技能失败:', error);
  }

  return skills;
}

// ============================================================
// 工具函数
// ============================================================

/**
 * 检查全局技能是否已初始化
 */
export function isGlobalSkillsInitialized(): boolean {
  return existsSync(getGlobalSkillsMetaPath());
}

/**
 * 获取全局技能元数据
 */
export function getGlobalSkillsMeta(): GlobalSkillsMetadata | null {
  return readGlobalSkillsMeta();
}
