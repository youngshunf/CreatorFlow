/**
 * Bundled Applications
 *
 * Built-in applications that ship with CreatorFlow.
 * These are registered at startup and provide default functionality.
 */

import { existsSync, mkdirSync, readdirSync, copyFileSync, readFileSync, writeFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';
import type { AppManifest } from './types.ts';
import { registerBundledApp, getAppPath, ensureAppsDir } from './storage.ts';

// ============================================================
// General App (Default)
// ============================================================

/**
 * General purpose app for workspaces without a specific application.
 * Used for backward compatibility and custom workspaces.
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

// ============================================================
// Creator Media App (自媒体创作)
// ============================================================

/**
 * Content creation app for social media creators.
 * Supports Xiaohongshu, Douyin, WeChat Official Accounts, etc.
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

// ============================================================
// Skill Copying Utilities
// ============================================================

/**
 * Get the path to bundled-skills directory
 * Works in both ESM and CJS contexts
 */
function getBundledSkillsDir(): string {
  // Try to get current file directory
  try {
    // ESM context
    const currentFile = fileURLToPath(import.meta.url);
    return join(dirname(currentFile), 'bundled-skills');
  } catch {
    // Fallback for CJS or other contexts
    return join(__dirname, 'bundled-skills');
  }
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

/**
 * Register all bundled apps.
 * Called at application startup.
 * Also copies bundled skills to app directories.
 */
export function registerBundledApps(): void {
  ensureAppsDir();
  
  // Register General App
  registerBundledApp(GENERAL_APP);
  copyBundledSkillsToApp(GENERAL_APP.id, GENERAL_APP.capabilities?.skills || []);
  
  // Register Creator Media App
  registerBundledApp(CREATOR_MEDIA_APP);
  copyBundledSkillsToApp(CREATOR_MEDIA_APP.id, CREATOR_MEDIA_APP.capabilities?.skills || []);
}

/**
 * Get all bundled app manifests
 */
export function getBundledAppManifests(): AppManifest[] {
  return [GENERAL_APP, CREATOR_MEDIA_APP];
}
