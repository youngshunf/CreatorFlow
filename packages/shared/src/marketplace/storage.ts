/**
 * Marketplace Local Storage
 *
 * Manages local cache and version tracking for marketplace items.
 * Handles:
 * - Marketplace metadata cache (~/.creator-flow/marketplace/cache/)
 * - Skill meta files ({workspace}/skills/{skill_id}/.skill-meta.json)
 * - Installed package tracking
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  readdirSync,
  statSync,
} from 'fs';
import { join, dirname } from 'path';
import { homedir } from 'os';
import type {
  SkillMeta,
  MarketplaceCacheData,
  MarketplaceSkill,
  MarketplaceApp,
  MarketplaceCategory,
  InstalledSkillInfo,
} from './types.ts';

// ============================================================
// Path Utilities
// ============================================================

/**
 * Get the root config directory
 */
export function getConfigRoot(): string {
  return join(homedir(), '.creator-flow');
}

/**
 * Get the marketplace directory
 */
export function getMarketplaceDir(): string {
  return join(getConfigRoot(), 'marketplace');
}

/**
 * Get the marketplace cache directory
 */
export function getMarketplaceCacheDir(): string {
  return join(getMarketplaceDir(), 'cache');
}

/**
 * Get the marketplace skills cache directory
 */
export function getMarketplaceSkillsDir(): string {
  return join(getMarketplaceDir(), 'skills');
}

/**
 * Get the marketplace apps cache directory
 */
export function getMarketplaceAppsDir(): string {
  return join(getMarketplaceDir(), 'apps');
}

/**
 * Ensure directory exists
 */
function ensureDir(dir: string): void {
  if (!existsSync(dir)) {
    mkdirSync(dir, { recursive: true });
  }
}

/**
 * Expand ~ to home directory in path
 */
function expandPath(path: string): string {
  if (path.startsWith('~')) {
    return path.replace(/^~/, homedir());
  }
  return path;
}

// ============================================================
// Marketplace Cache Operations
// ============================================================

const CACHE_FILE = 'meta.json';
const CACHE_TTL_MS = 4 * 60 * 60 * 1000; // 4 hours

/**
 * Get cached marketplace data
 */
export function getMarketplaceCache(): MarketplaceCacheData | null {
  const cacheFile = join(getMarketplaceCacheDir(), CACHE_FILE);
  if (!existsSync(cacheFile)) {
    return null;
  }

  try {
    const content = readFileSync(cacheFile, 'utf-8');
    const data = JSON.parse(content) as MarketplaceCacheData;

    // Check if cache is expired
    const cacheTime = new Date(data.cachedAt).getTime();
    if (Date.now() - cacheTime > CACHE_TTL_MS) {
      return null; // Cache expired
    }

    return data;
  } catch {
    return null;
  }
}

/**
 * Save marketplace data to cache
 */
export function saveMarketplaceCache(data: {
  skills: MarketplaceSkill[];
  apps: MarketplaceApp[];
  categories: MarketplaceCategory[];
}): void {
  const cacheDir = getMarketplaceCacheDir();
  ensureDir(cacheDir);

  const cacheData: MarketplaceCacheData = {
    cachedAt: new Date().toISOString(),
    skills: data.skills,
    apps: data.apps,
    categories: data.categories,
  };

  const cacheFile = join(cacheDir, CACHE_FILE);
  writeFileSync(cacheFile, JSON.stringify(cacheData, null, 2));
}

/**
 * Clear marketplace cache
 */
export function clearMarketplaceCache(): void {
  const cacheFile = join(getMarketplaceCacheDir(), CACHE_FILE);
  if (existsSync(cacheFile)) {
    const { rmSync } = require('fs');
    rmSync(cacheFile);
  }
}

// ============================================================
// Skill Meta Operations
// ============================================================

const SKILL_META_FILE = '.skill-meta.json';

/**
 * Get skill meta file path
 */
export function getSkillMetaPath(skillDir: string): string {
  return join(skillDir, SKILL_META_FILE);
}

/**
 * Read skill meta from a workspace skill directory
 */
export function readSkillMeta(skillDir: string): SkillMeta | null {
  const metaPath = getSkillMetaPath(skillDir);
  if (!existsSync(metaPath)) {
    return null;
  }

  try {
    const content = readFileSync(metaPath, 'utf-8');
    return JSON.parse(content) as SkillMeta;
  } catch {
    return null;
  }
}

/**
 * Write skill meta to a workspace skill directory
 */
export function writeSkillMeta(skillDir: string, meta: SkillMeta): void {
  ensureDir(skillDir);
  const metaPath = getSkillMetaPath(skillDir);
  writeFileSync(metaPath, JSON.stringify(meta, null, 2));
}

/**
 * Create initial skill meta for a newly installed skill
 */
export function createSkillMeta(
  sourceId: string,
  version: string
): SkillMeta {
  return {
    sourceId,
    baseVersion: version,
    localVersion: version,
    isModified: false,
    lastSynced: new Date().toISOString(),
    isPrivate: false,
  };
}

/**
 * Mark skill as modified locally
 */
export function markSkillModified(skillDir: string): void {
  const meta = readSkillMeta(skillDir);
  if (meta) {
    // Increment local version
    const parts = meta.localVersion.split('-local.');
    const baseVersion = parts[0];
    const localNum = parts.length > 1 ? parseInt(parts[1], 10) + 1 : 1;
    
    meta.localVersion = `${baseVersion}-local.${localNum}`;
    meta.isModified = true;
    writeSkillMeta(skillDir, meta);
  }
}

// ============================================================
// Installed Skills Tracking
// ============================================================

/**
 * Get all installed skills with their meta info from a workspace
 * @param workspaceRoot - Absolute path to workspace root
 */
export function getInstalledSkills(
  workspaceRoot: string
): InstalledSkillInfo[] {
  const expandedRoot = expandPath(workspaceRoot);
  const skillsDir = join(expandedRoot, '.creator-flow', 'skills');
  if (!existsSync(skillsDir)) {
    return [];
  }

  const result: InstalledSkillInfo[] = [];

  try {
    const entries = readdirSync(skillsDir, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const skillDir = join(skillsDir, entry.name);
      const skillFile = join(skillDir, 'SKILL.md');

      // Skip if no SKILL.md
      if (!existsSync(skillFile)) continue;

      // Read meta if exists
      const meta = readSkillMeta(skillDir);
      
      // Parse SKILL.md for name and description
      let name = entry.name;
      let description = '';
      try {
        const content = readFileSync(skillFile, 'utf-8');
        const frontmatterMatch = content.match(/^---\n([\s\S]*?)\n---/);
        if (frontmatterMatch) {
          const frontmatter = frontmatterMatch[1];
          const nameMatch = frontmatter.match(/name:\s*(.+)/);
          const descMatch = frontmatter.match(/description:\s*(.+)/);
          if (nameMatch) name = nameMatch[1].trim();
          if (descMatch) description = descMatch[1].trim();
        }
      } catch {
        // Ignore parsing errors
      }

      // Get version: prefer meta, fallback to config.yaml
      let version = meta?.localVersion;
      if (!version) {
        const configPath = join(skillDir, 'config.yaml');
        if (existsSync(configPath)) {
          try {
            const configContent = readFileSync(configPath, 'utf-8');
            const versionMatch = configContent.match(/^version:\s*(.+)$/m);
            if (versionMatch) version = versionMatch[1].trim();
          } catch {
            // Ignore parsing errors
          }
        }
      }

      result.push({
        skillId: meta?.sourceId || entry.name,
        name,
        description,
        version: version || 'unknown',
        isModified: meta?.isModified || false,
        hasUpdate: false, // Will be updated by sync
        path: skillDir,
      });
    }
  } catch {
    // Ignore errors
  }

  return result;
}

/**
 * Check if a skill is installed in a workspace
 */
export function isSkillInstalled(
  workspaceRoot: string,
  skillId: string
): boolean {
  const skillDir = join(workspaceRoot, '.creator-flow', 'skills', skillId);
  const skillFile = join(skillDir, 'SKILL.md');
  return existsSync(skillFile);
}

/**
 * Get installed skill version
 */
export function getInstalledSkillVersion(
  workspaceRoot: string,
  skillId: string
): string | null {
  const skillDir = join(workspaceRoot, '.creator-flow', 'skills', skillId);
  const meta = readSkillMeta(skillDir);
  return meta?.baseVersion || null;
}

// ============================================================
// Downloaded Package Cache
// ============================================================

/**
 * Get path to cached skill package
 */
export function getCachedSkillPackagePath(
  skillId: string,
  version: string
): string {
  return join(getMarketplaceSkillsDir(), skillId, version);
}

/**
 * Check if skill package is cached
 */
export function isSkillPackageCached(
  skillId: string,
  version: string
): boolean {
  const packageDir = getCachedSkillPackagePath(skillId, version);
  return existsSync(packageDir) && existsSync(join(packageDir, 'SKILL.md'));
}

/**
 * Get path to cached app package
 */
export function getCachedAppPackagePath(
  appId: string,
  version: string
): string {
  return join(getMarketplaceAppsDir(), appId, version);
}

/**
 * Check if app package is cached
 */
export function isAppPackageCached(appId: string, version: string): boolean {
  const packageDir = getCachedAppPackagePath(appId, version);
  return existsSync(packageDir) && existsSync(join(packageDir, 'manifest.json'));
}
