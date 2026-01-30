/**
 * Bundled Applications
 *
 * Built-in applications that ship with CreatorFlow.
 * These are registered at startup and provide default functionality.
 * 
 * Applications are loaded from bundled-apps/ directory (manifest.json files).
 * Each app can have an AGENTS.md file that provides context for the AI agent.
 */

import { existsSync, mkdirSync, readdirSync, copyFileSync, readFileSync, writeFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';
import type { AppManifest } from './types.ts';
import { registerBundledApp, getAppPath, ensureAppsDir } from './storage.ts';
import { debug } from '../utils/debug.ts';

// ============================================================
// Bundled Apps Directory Discovery
// ============================================================

/**
 * Get the path to bundled-apps directory.
 * Works in both development and packaged (Electron) contexts.
 */
function getBundledAppsSourceDir(): string {
  const possiblePaths: string[] = [];
  
  // 1. Try ESM context (development - source directory)
  try {
    const currentFile = fileURLToPath(import.meta.url);
    possiblePaths.push(join(dirname(currentFile), 'bundled-apps'));
  } catch {
    // ESM not available
  }
  
  // 2. Try CJS context (development - source directory)
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'bundled-apps'));
  }
  
  // 3. Try dist/assets/bundled-apps (Electron dev mode after copy-assets)
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'assets', 'bundled-apps'));
    possiblePaths.push(join(__dirname, '..', 'assets', 'bundled-apps'));
  }
  
  // 4. Try packaged Electron app paths (resources/bundled-apps)
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'resources', 'bundled-apps'));
    possiblePaths.push(join(__dirname, '..', 'resources', 'bundled-apps'));
  }
  
  // 5. Try process.resourcesPath (Electron specific)
  if (typeof process !== 'undefined' && (process as any).resourcesPath) {
    possiblePaths.push(join((process as any).resourcesPath, 'bundled-apps'));
  }
  
  // Return the first path that exists
  for (const p of possiblePaths) {
    if (existsSync(p)) {
      debug(`[getBundledAppsSourceDir] Found bundled-apps at: ${p}`);
      return p;
    }
  }
  
  // Fallback to first possible path
  return possiblePaths[0] || '';
}

/**
 * Load an app manifest from a directory.
 * Returns null if manifest doesn't exist or is invalid.
 */
function loadAppManifestFromDir(appDir: string): AppManifest | null {
  const manifestPath = join(appDir, 'manifest.json');
  
  if (!existsSync(manifestPath)) {
    return null;
  }
  
  try {
    const content = readFileSync(manifestPath, 'utf-8');
    const manifest = JSON.parse(content) as AppManifest;
    
    // Validate required fields
    if (!manifest.id || !manifest.name || !manifest.version) {
      debug(`[loadAppManifestFromDir] Invalid manifest at ${manifestPath}: missing required fields`);
      return null;
    }
    
    // Set default type
    if (!manifest.type) {
      manifest.type = manifest.id.startsWith('plugin.') ? 'plugin' : 'app';
    }
    
    return manifest;
  } catch (error) {
    debug(`[loadAppManifestFromDir] Error loading manifest from ${manifestPath}:`, error);
    return null;
  }
}

/**
 * Load all bundled apps from the bundled-apps directory.
 * Returns an array of app manifests with their source paths.
 */
function loadBundledAppsFromDir(): Array<{ manifest: AppManifest; sourcePath: string }> {
  const bundledAppsDir = getBundledAppsSourceDir();
  const apps: Array<{ manifest: AppManifest; sourcePath: string }> = [];
  
  if (!bundledAppsDir || !existsSync(bundledAppsDir)) {
    debug(`[loadBundledAppsFromDir] Bundled apps directory not found: ${bundledAppsDir}`);
    return apps;
  }
  
  try {
    const entries = readdirSync(bundledAppsDir, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;
      
      const appDir = join(bundledAppsDir, entry.name);
      const manifest = loadAppManifestFromDir(appDir);
      
      if (manifest) {
        apps.push({ manifest, sourcePath: appDir });
        debug(`[loadBundledAppsFromDir] Loaded app: ${manifest.id} from ${appDir}`);
      }
    }
  } catch (error) {
    debug(`[loadBundledAppsFromDir] Error reading bundled apps directory:`, error);
  }
  
  return apps;
}

/**
 * Get the path to an app's AGENTS.md file in the source directory.
 * Returns null if the file doesn't exist.
 */
export function getAppAgentsFilePath(appSourcePath: string): string | null {
  const agentsPath = join(appSourcePath, 'AGENTS.md');
  return existsSync(agentsPath) ? agentsPath : null;
}

// ============================================================
// Fallback Hard-coded App Manifests (for backward compatibility)
// ============================================================

/**
 * General purpose app for workspaces without a specific application.
 * Used for backward compatibility and custom workspaces.
 * @deprecated Use manifest.json in bundled-apps/app-general/ instead
 */
export const GENERAL_APP: AppManifest = {
  id: 'app.general',
  name: '通用工作区',
  version: '1.0.0',
  type: 'app',
  author: {
    name: 'CreatorFlow Team',
    email: 'team@creatorflow.com',
  },
  description: '通用 AI 助手工作区，适合各种任务场景',
  license: 'MIT',

  pricing: {
    type: 'free',
  },

  capabilities: {
    skills: [],
  },

  workspace: {
    defaultDirectory: '~/CreatorFlow',
    directoryStructure: {},
    defaultSettings: {
      autoSave: true,
    },
  },

  permissions: ['filesystem', 'network'],

  compatibility: {
    platform: '>=1.0.0',
    os: ['macos', 'windows', 'linux'],
  },

  marketplace: {
    category: 'general',
    tags: ['通用', 'AI助手'],
  },
};

/**
 * Content creation app for social media creators.
 * Supports Xiaohongshu, Douyin, WeChat Official Accounts, etc.
 * @deprecated Use manifest.json in bundled-apps/app-creator-media/ instead
 */
export const CREATOR_MEDIA_APP: AppManifest = {
  id: 'app.creator-media',
  name: '自媒体创作',
  version: '1.0.0',
  type: 'app',
  author: {
    name: 'CreatorFlow Team',
    email: 'team@creatorflow.com',
  },
  description: '完整的自媒体创作工作流，支持小红书、抖音、公众号等平台',
  license: 'MIT',
  homepage: 'https://creatorflow.com/apps/creator-media',

  pricing: {
    type: 'free',
  },

  capabilities: {
    skills: [
      'material-organize',
      'script-create',
      'topic-research',
      'content-review',
      'platform-adaptation',
      'viral-content-analysis',
      'topic-to-publish-flow',
      'data-review',
    ],
  },

  workspace: {
    defaultDirectory: '~/CreatorFlow/Media',
    directoryStructure: {
      '素材库': {
        '图片': {},
        '视频': {},
        '文本': {},
        '音频': {},
      },
      '脚本库': {
        '草稿': {},
        '待审核': {},
        '已发布': {},
      },
      '发布记录': {
        '小红书': {},
        '抖音': {},
        '公众号': {},
      },
      '数据分析': {},
    },
    defaultSettings: {
      autoSave: true,
      autoBackup: true,
      backupInterval: 3600,
      language: 'zh-CN',
    },
    presetData: {
      labels: [
        { name: '待创作', color: '#3b82f6' },
        { name: '创作中', color: '#f59e0b' },
        { name: '待审核', color: '#8b5cf6' },
        { name: '已发布', color: '#10b981' },
      ],
      statuses: [
        { name: '草稿', type: 'draft' },
        { name: '待审核', type: 'review' },
        { name: '已发布', type: 'published' },
      ],
    },
  },

  permissions: ['filesystem', 'network', 'clipboard'],

  compatibility: {
    platform: '>=1.0.0',
    os: ['macos', 'windows', 'linux'],
  },

  marketplace: {
    category: 'content-creation',
    tags: ['自媒体', '创作', 'AI', '小红书', '抖音', '公众号', '内容创作'],
    changelog: `## v1.0.0
- 初始版本发布
- 支持小红书、抖音、公众号、B站等主流平台
- 完整的创作工作流支持
- 内置 8 个核心技能：素材整理、脚本创作、选题研究、内容审核、平台适配、爆款分析、选题到发布流程、数据复盘`,
  },

  plugins: {
    compatible: [
      'plugin.data-analytics',
      'plugin.seo-optimizer',
    ],
    recommended: [
      'plugin.data-analytics',
    ],
  },
};

/** Fallback app manifests (used when bundled-apps directory is not available) */
const FALLBACK_APPS: AppManifest[] = [GENERAL_APP, CREATOR_MEDIA_APP];

// ============================================================
// Skill Copying Utilities
// ============================================================

/**
 * Get the path to bundled-skills directory
 * Works in both development and packaged (Electron) contexts
 */
function getBundledSkillsDir(): string {
  const possiblePaths: string[] = [];
  
  // 1. Try ESM context (development)
  try {
    const currentFile = fileURLToPath(import.meta.url);
    possiblePaths.push(join(dirname(currentFile), 'bundled-skills'));
  } catch {
    // ESM not available
  }
  
  // 2. Try CJS context (development)
  if (typeof __dirname !== 'undefined') {
    possiblePaths.push(join(__dirname, 'bundled-skills'));
  }
  
  // 3. Try packaged Electron app paths (resources/bundled-skills)
  if (typeof __dirname !== 'undefined') {
    // In packaged app, __dirname might be in dist/, look for resources/bundled-skills
    possiblePaths.push(join(__dirname, 'resources', 'bundled-skills'));
    possiblePaths.push(join(__dirname, '..', 'resources', 'bundled-skills'));
  }
  
  // 4. Try process.resourcesPath (Electron specific)
  if (typeof process !== 'undefined' && (process as any).resourcesPath) {
    possiblePaths.push(join((process as any).resourcesPath, 'bundled-skills'));
  }
  
  // Return the first path that exists
  for (const p of possiblePaths) {
    if (existsSync(p)) {
      return p;
    }
  }
  
  // Fallback to first possible path (will fail gracefully in copyBundledSkillsToApp)
  return possiblePaths[0] || '';
}

/**
 * Copy a directory recursively
 */
function copyDirectoryRecursive(source: string, target: string): void {
  if (!existsSync(source)) return;
  
  mkdirSync(target, { recursive: true });
  
  const entries = readdirSync(source, { withFileTypes: true });
  for (const entry of entries) {
    const sourcePath = join(source, entry.name);
    const targetPath = join(target, entry.name);
    
    if (entry.isFile()) {
      copyFileSync(sourcePath, targetPath);
    } else if (entry.isDirectory()) {
      copyDirectoryRecursive(sourcePath, targetPath);
    }
  }
}

/**
 * Copy bundled skills to app directory
 */
function copyBundledSkillsToApp(appId: string, skillSlugs: string[]): void {
  if (skillSlugs.length === 0) return;
  
  const bundledSkillsDir = getBundledSkillsDir();
  const appPath = getAppPath(appId, true);
  const appSkillsDir = join(appPath, 'skills');
  
  // Ensure skills directory exists
  mkdirSync(appSkillsDir, { recursive: true });
  
  for (const skillSlug of skillSlugs) {
    const sourceSkillPath = join(bundledSkillsDir, skillSlug);
    const targetSkillPath = join(appSkillsDir, skillSlug);
    
    if (existsSync(sourceSkillPath)) {
      copyDirectoryRecursive(sourceSkillPath, targetSkillPath);
    }
  }
}

// ============================================================
// Registration
// ============================================================

/** Cache for loaded bundled apps (with source paths) */
let loadedBundledApps: Array<{ manifest: AppManifest; sourcePath: string }> | null = null;

/**
 * Register all bundled apps.
 * Called at application startup.
 * Loads apps from bundled-apps directory (config files), falls back to hard-coded manifests.
 * Also copies bundled skills and AGENTS.md to app directories.
 */
export function registerBundledApps(): void {
  ensureAppsDir();
  
  // Try to load apps from bundled-apps directory first
  loadedBundledApps = loadBundledAppsFromDir();
  
  if (loadedBundledApps.length > 0) {
    // Register apps loaded from config files
    for (const { manifest, sourcePath } of loadedBundledApps) {
      registerBundledApp(manifest);
      copyBundledSkillsToApp(manifest.id, manifest.capabilities?.skills || []);
      
      // Copy AGENTS.md to app directory if it exists
      copyAppAgentsFile(manifest.id, sourcePath);
    }
    debug(`[registerBundledApps] Registered ${loadedBundledApps.length} apps from config files`);
  } else {
    // Fallback to hard-coded manifests
    debug(`[registerBundledApps] Using fallback hard-coded manifests`);
    
    for (const manifest of FALLBACK_APPS) {
      registerBundledApp(manifest);
      copyBundledSkillsToApp(manifest.id, manifest.capabilities?.skills || []);
    }
  }
}

/**
 * Copy AGENTS.md from app source directory to installed app directory.
 */
function copyAppAgentsFile(appId: string, sourcePath: string): void {
  const agentsSourcePath = join(sourcePath, 'AGENTS.md');
  if (!existsSync(agentsSourcePath)) {
    return;
  }
  
  const appPath = getAppPath(appId, true);
  const agentsTargetPath = join(appPath, 'AGENTS.md');
  
  try {
    copyFileSync(agentsSourcePath, agentsTargetPath);
    debug(`[copyAppAgentsFile] Copied AGENTS.md for ${appId}`);
  } catch (error) {
    debug(`[copyAppAgentsFile] Error copying AGENTS.md for ${appId}:`, error);
  }
}

/**
 * Get all bundled app manifests.
 * Returns apps loaded from config files, or fallback hard-coded manifests.
 */
export function getBundledAppManifests(): AppManifest[] {
  if (loadedBundledApps && loadedBundledApps.length > 0) {
    return loadedBundledApps.map(({ manifest }) => manifest);
  }
  return FALLBACK_APPS;
}

/**
 * Get the source path for a bundled app (where AGENTS.md and other resources are located).
 * Returns null if the app is not found or was loaded from fallback.
 */
export function getBundledAppSourcePath(appId: string): string | null {
  if (!loadedBundledApps) {
    return null;
  }
  
  const app = loadedBundledApps.find(({ manifest }) => manifest.id === appId);
  return app?.sourcePath || null;
}
