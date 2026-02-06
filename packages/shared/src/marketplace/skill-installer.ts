/**
 * Skill Installer
 *
 * Handles downloading and installing skills from the marketplace.
 * Supports:
 * - Download skill packages from cloud
 * - Extract and validate packages
 * - Install to workspace
 * - Version tracking
 */

import {
  existsSync,
  mkdirSync,
  writeFileSync,
  readFileSync,
  readdirSync,
  copyFileSync,
  rmSync,
  createWriteStream,
} from 'fs';
import { join, basename, dirname } from 'path';
import { pipeline } from 'stream/promises';
import { createHash } from 'crypto';
import { Extract } from 'unzipper';
import type {
  InstallProgress,
  InstallProgressCallback,
  InstallResult,
  DownloadResponse,
} from './types.ts';
import { getSkillDownload } from './api.ts';
import {
  getCachedSkillPackagePath,
  isSkillPackageCached,
  createSkillMeta,
  writeSkillMeta,
  getMarketplaceSkillsDir,
} from './storage.ts';

// ============================================================
// File Search Utilities
// ============================================================

/**
 * 递归查找文件
 * @param dir - 搜索目录
 * @param filename - 目标文件名
 * @param maxDepth - 最大搜索深度
 * @returns 文件路径或 null
 */
function findFileRecursive(dir: string, filename: string, maxDepth: number): string | null {
  if (maxDepth <= 0) return null;

  try {
    const entries = readdirSync(dir, { withFileTypes: true });

    for (const entry of entries) {
      const fullPath = join(dir, entry.name);

      if (entry.isFile() && entry.name === filename) {
        return fullPath;
      }

      if (entry.isDirectory()) {
        const found = findFileRecursive(fullPath, filename, maxDepth - 1);
        if (found) return found;
      }
    }
  } catch (error) {
    // 忽略读取错误
  }

  return null;
}

// ============================================================
// Download Utilities
// ============================================================

/**
 * Download a file from URL to local path
 */
async function downloadFile(
  url: string,
  destPath: string,
  onProgress?: (percent: number) => void
): Promise<void> {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Download failed: ${response.status} ${response.statusText}`);
  }

  const contentLength = response.headers.get('content-length');
  const totalSize = contentLength ? parseInt(contentLength, 10) : 0;

  const fileStream = createWriteStream(destPath);
  const reader = response.body?.getReader();
  
  if (!reader) {
    throw new Error('Unable to read response body');
  }

  let downloadedSize = 0;

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    fileStream.write(value);
    downloadedSize += value.length;

    if (totalSize > 0 && onProgress) {
      const percent = Math.round((downloadedSize / totalSize) * 100);
      onProgress(percent);
    }
  }

  fileStream.end();
  
  // Wait for file to be fully written
  await new Promise<void>((resolve, reject) => {
    fileStream.on('finish', resolve);
    fileStream.on('error', reject);
  });
}

/**
 * Calculate SHA256 hash of a file
 */
function calculateFileHash(filePath: string): string {
  const content = readFileSync(filePath);
  return createHash('sha256').update(content).digest('hex');
}

/**
 * Extract a zip file to a directory
 */
async function extractZip(zipPath: string, destDir: string): Promise<void> {
  // Ensure destination exists
  if (!existsSync(destDir)) {
    mkdirSync(destDir, { recursive: true });
  }

  // Extract using unzipper
  const extract = Extract({ path: destDir });
  const readStream = require('fs').createReadStream(zipPath);
  
  await pipeline(readStream, extract);
}

/**
 * Copy directory recursively (overwrite existing files)
 */
function copyDir(src: string, dest: string): void {
  if (!existsSync(dest)) {
    mkdirSync(dest, { recursive: true });
  }

  const entries = readdirSync(src, { withFileTypes: true });
  for (const entry of entries) {
    const srcPath = join(src, entry.name);
    const destPath = join(dest, entry.name);

    if (entry.isDirectory()) {
      copyDir(srcPath, destPath);
    } else {
      copyFileSync(srcPath, destPath);
    }
  }
}

/**
 * Merge directory - copy new files and overwrite existing, but keep files not in source
 */
function mergeDir(src: string, dest: string): void {
  if (!existsSync(dest)) {
    mkdirSync(dest, { recursive: true });
  }

  const entries = readdirSync(src, { withFileTypes: true });
  for (const entry of entries) {
    const srcPath = join(src, entry.name);
    const destPath = join(dest, entry.name);

    if (entry.isDirectory()) {
      mergeDir(srcPath, destPath);
    } else {
      // Overwrite existing file or create new
      copyFileSync(srcPath, destPath);
    }
  }
}

// ============================================================
// Skill Installation
// ============================================================

/**
 * Download and cache a skill package
 */
export async function downloadSkillPackage(
  skillId: string,
  version: string = 'latest',
  onProgress?: InstallProgressCallback
): Promise<{ packageDir: string; actualVersion: string }> {
  onProgress?.({
    stage: 'downloading',
    percent: 0,
    message: `正在获取技能信息: ${skillId}`,
    skillId,
  });

  // Get download URL from API
  const downloadInfo = await getSkillDownload(skillId, version);
  
  if (!downloadInfo.download_url) {
    throw new Error('无法获取下载链接');
  }

  // Determine actual version from URL or use 'latest'
  const actualVersion = version === 'latest' 
    ? extractVersionFromUrl(downloadInfo.download_url) || 'latest'
    : version;

  // Check if already cached
  if (isSkillPackageCached(skillId, actualVersion)) {
    onProgress?.({
      stage: 'complete',
      percent: 100,
      message: '已使用缓存版本',
      skillId,
    });
    return {
      packageDir: getCachedSkillPackagePath(skillId, actualVersion),
      actualVersion,
    };
  }

  // Create cache directory - ensure all parent directories exist
  const cacheDir = getMarketplaceSkillsDir();
  if (!existsSync(cacheDir)) {
    mkdirSync(cacheDir, { recursive: true });
  }
  const skillCacheDir = join(cacheDir, skillId);
  const versionDir = join(skillCacheDir, actualVersion);
  const zipPath = join(skillCacheDir, `${actualVersion}.zip`);

  if (!existsSync(skillCacheDir)) {
    mkdirSync(skillCacheDir, { recursive: true });
  }

  // Download the package
  onProgress?.({
    stage: 'downloading',
    percent: 10,
    message: `正在下载: ${skillId}@${actualVersion}`,
    skillId,
  });

  await downloadFile(downloadInfo.download_url, zipPath, (percent) => {
    onProgress?.({
      stage: 'downloading',
      percent: 10 + Math.round(percent * 0.6), // 10-70%
      message: `正在下载: ${percent}%`,
      skillId,
    });
  });

  // Verify hash if provided
  if (downloadInfo.file_hash) {
    const actualHash = calculateFileHash(zipPath);
    if (actualHash !== downloadInfo.file_hash) {
      rmSync(zipPath);
      throw new Error('文件校验失败，下载可能已损坏');
    }
  }

  // Extract package
  onProgress?.({
    stage: 'extracting',
    percent: 75,
    message: '正在解压...',
    skillId,
  });

  // Ensure version directory exists before extraction
  if (!existsSync(versionDir)) {
    mkdirSync(versionDir, { recursive: true });
  }

  await extractZip(zipPath, versionDir);

  // Clean up zip file
  rmSync(zipPath);

  // Validate extracted content with recursive search
  const skillFile = join(versionDir, 'SKILL.md');
  if (!existsSync(skillFile)) {
    // Recursively search for SKILL.md (max 2 levels deep)
    const found = findFileRecursive(versionDir, 'SKILL.md', 2);

    if (found) {
      // Move contents from nested directory to root
      const parentDir = dirname(found);
      if (parentDir !== versionDir) {
        const entries = readdirSync(parentDir);
        for (const entry of entries) {
          const src = join(parentDir, entry);
          const dest = join(versionDir, entry);
          require('fs').renameSync(src, dest);
        }
        // Remove empty parent directory
        rmSync(parentDir, { recursive: true, force: true });
      }
    }
  }

  // Final validation
  if (!existsSync(join(versionDir, 'SKILL.md'))) {
    rmSync(versionDir, { recursive: true });
    throw new Error('无效的技能包：缺少 SKILL.md 文件');
  }

  onProgress?.({
    stage: 'complete',
    percent: 100,
    message: '下载完成',
    skillId,
  });

  return { packageDir: versionDir, actualVersion };
}

/**
 * Install a skill to a workspace
 */
export async function installSkill(
  workspaceRoot: string,
  skillId: string,
  version: string = 'latest',
  onProgress?: InstallProgressCallback
): Promise<InstallResult> {
  try {
    // Download package
    const { packageDir, actualVersion } = await downloadSkillPackage(
      skillId,
      version,
      onProgress
    );

    // Install to workspace
    onProgress?.({
      stage: 'installing',
      percent: 85,
      message: '正在安装到工作区...',
      skillId,
    });

    const skillsDir = join(workspaceRoot, '.creator-flow', 'skills');
    const targetDir = join(skillsDir, skillId);

    // Ensure skills directory exists
    if (!existsSync(skillsDir)) {
      mkdirSync(skillsDir, { recursive: true });
    }

    // Merge with existing skill directory (overwrite same-name files, keep others)
    // This preserves user customizations that don't conflict with the package
    mergeDir(packageDir, targetDir);

    // Write skill meta
    const meta = createSkillMeta(skillId, actualVersion);
    writeSkillMeta(targetDir, meta);

    onProgress?.({
      stage: 'complete',
      percent: 100,
      message: '安装完成',
      skillId,
    });

    return {
      success: true,
      skillId,
      version: actualVersion,
      path: targetDir,
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : '未知错误';
    
    onProgress?.({
      stage: 'error',
      percent: 0,
      message: errorMessage,
      skillId,
      error: errorMessage,
    });

    return {
      success: false,
      skillId,
      version,
      error: errorMessage,
    };
  }
}

/**
 * Install multiple skills (for app dependencies)
 */
export async function installSkills(
  workspaceRoot: string,
  skills: Array<{ skillId: string; version?: string }>,
  onProgress?: InstallProgressCallback
): Promise<InstallResult[]> {
  const results: InstallResult[] = [];
  const total = skills.length;

  for (let i = 0; i < skills.length; i++) {
    const skill = skills[i];
    if (!skill) continue;
    
    const { skillId, version } = skill;
    
    const result = await installSkill(
      workspaceRoot,
      skillId,
      version || 'latest',
      (progress) => {
        // Adjust progress for batch install
        const batchPercent = ((i + progress.percent / 100) / total) * 100;
        onProgress?.({
          ...progress,
          percent: Math.round(batchPercent),
          message: `[${i + 1}/${total}] ${progress.message}`,
        });
      }
    );

    results.push(result);
  }

  return results;
}

/**
 * Update an installed skill to a new version
 */
export async function updateSkill(
  workspaceRoot: string,
  skillId: string,
  targetVersion: string = 'latest',
  onProgress?: InstallProgressCallback
): Promise<InstallResult> {
  // Simply reinstall - the installer handles existing skill replacement
  return installSkill(workspaceRoot, skillId, targetVersion, onProgress);
}

// ============================================================
// Utility Functions
// ============================================================

/**
 * Extract version from download URL
 * e.g., "https://xxx/skills/material-organize/1.0.0.zip" -> "1.0.0"
 */
function extractVersionFromUrl(url: string): string | null {
  const match = url.match(/\/(\d+\.\d+\.\d+)\.zip/);
  return match ? match[1] : null;
}
