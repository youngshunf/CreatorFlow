/**
 * App Installer
 *
 * Handles downloading and installing apps from the marketplace.
 * Supports:
 * - Download app packages from cloud
 * - Extract and validate packages
 * - Parse skill dependencies
 * - Download and install required skills
 * - Initialize workspace with app
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
import { join } from 'path';
import { pipeline } from 'stream/promises';
import { createHash } from 'crypto';
import { Extract } from 'unzipper';
import type {
  InstallProgress,
  AppInstallProgress,
  InstallProgressCallback,
  AppInstallProgressCallback,
  AppDownloadResponse,
  SkillDependencyDownload,
} from './types.ts';
import { getAppDownload } from './api.ts';
import {
  getCachedAppPackagePath,
  isAppPackageCached,
  getMarketplaceAppsDir,
} from './storage.ts';
import { installSkill } from './skill-installer.ts';
import { getGlobalSkillsDir } from '../skills/global-skills.ts';

// ============================================================
// Types
// ============================================================

export interface AppInstallResult {
  success: boolean;
  appId: string;
  version: string;
  appPath?: string;
  skillResults?: Array<{
    skillId: string;
    success: boolean;
    error?: string;
  }>;
  error?: string;
  /** Backup path if existing app was backed up */
  backupPath?: string;
}

/**
 * Installed app information
 */
export interface InstalledAppInfo {
  appId: string;
  name: string;
  version: string;
  description?: string;
  installedAt?: string;
  skillDependencies?: string[];
}

/**
 * App install options
 */
export interface AppInstallOptions {
  /** Force install even if app already exists (will backup first) */
  force?: boolean;
  /** Merge install - keep existing files and only add/update from package */
  merge?: boolean;
  /** Skip backup when force installing */
  skipBackup?: boolean;
  /** Progress callback */
  onProgress?: AppInstallProgressCallback;
}

/**
 * App uninstall options
 */
export interface AppUninstallOptions {
  /** Also remove skills that were installed with the app */
  removeSkills?: boolean;
  /** Create backup before uninstalling */
  backup?: boolean;
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
 * Copy directory recursively
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
      copyFileSync(srcPath, destPath);
    }
  }
}

// ============================================================
// App Status & Management
// ============================================================

/**
 * Check if an app is already installed in the workspace
 * @returns InstalledAppInfo if app exists, null otherwise
 */
export function checkInstalledApp(workspaceRoot: string): InstalledAppInfo | null {
  const creatorFlowDir = join(workspaceRoot, '.creator-flow');
  const manifestPath = join(creatorFlowDir, 'app-manifest.json');

  if (!existsSync(manifestPath)) {
    return null;
  }

  try {
    const manifest = JSON.parse(readFileSync(manifestPath, 'utf-8'));
    return {
      appId: manifest.id || 'unknown',
      name: manifest.name || manifest.id || 'Unknown App',
      version: manifest.version || 'unknown',
      description: manifest.description,
      installedAt: manifest.installed_at,
      skillDependencies: manifest.skill_dependencies,
    };
  } catch {
    return null;
  }
}

/**
 * Backup app data before reinstalling
 * @returns Backup directory path
 */
export function backupAppData(workspaceRoot: string): string | null {
  const creatorFlowDir = join(workspaceRoot, '.creator-flow');
  
  if (!existsSync(creatorFlowDir)) {
    return null;
  }

  // Create backup directory with timestamp
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  const backupDir = join(workspaceRoot, '.creator-flow-backup', timestamp);
  
  if (!existsSync(join(workspaceRoot, '.creator-flow-backup'))) {
    mkdirSync(join(workspaceRoot, '.creator-flow-backup'), { recursive: true });
  }

  // Copy directories that might have user customizations
  const dirsToBackup = ['prompts', 'guides', 'sources', 'skills'];
  let hasBackup = false;

  for (const dir of dirsToBackup) {
    const sourceDir = join(creatorFlowDir, dir);
    if (existsSync(sourceDir)) {
      const targetDir = join(backupDir, dir);
      copyDir(sourceDir, targetDir);
      hasBackup = true;
    }
  }

  // Also backup manifest
  const manifestPath = join(creatorFlowDir, 'app-manifest.json');
  if (existsSync(manifestPath)) {
    if (!existsSync(backupDir)) {
      mkdirSync(backupDir, { recursive: true });
    }
    copyFileSync(manifestPath, join(backupDir, 'app-manifest.json'));
    hasBackup = true;
  }

  return hasBackup ? backupDir : null;
}

/**
 * Uninstall an app from workspace
 */
export function uninstallApp(
  workspaceRoot: string,
  options: AppUninstallOptions = {}
): { success: boolean; backupPath?: string; error?: string } {
  const creatorFlowDir = join(workspaceRoot, '.creator-flow');

  // Check if app exists
  const installedApp = checkInstalledApp(workspaceRoot);
  if (!installedApp) {
    return { success: false, error: '没有找到已安装的应用' };
  }

  try {
    let backupPath: string | undefined;

    // Backup if requested
    if (options.backup) {
      backupPath = backupAppData(workspaceRoot) || undefined;
    }

    // Remove app-specific directories
    const dirsToRemove = ['prompts', 'guides', 'sources'];
    for (const dir of dirsToRemove) {
      const dirPath = join(creatorFlowDir, dir);
      if (existsSync(dirPath)) {
        rmSync(dirPath, { recursive: true });
      }
    }

    // Remove manifest
    const manifestPath = join(creatorFlowDir, 'app-manifest.json');
    if (existsSync(manifestPath)) {
      rmSync(manifestPath);
    }

    // Optionally remove skills
    if (options.removeSkills && installedApp.skillDependencies) {
      const skillsDir = join(creatorFlowDir, 'skills');
      for (const skillId of installedApp.skillDependencies) {
        const skillDir = join(skillsDir, skillId);
        if (existsSync(skillDir)) {
          rmSync(skillDir, { recursive: true });
        }
      }
    }

    return { success: true, backupPath };
  } catch (error) {
    return {
      success: false,
      error: error instanceof Error ? error.message : '卸载失败',
    };
  }
}

/**
 * Restore app data from backup
 */
export function restoreAppData(workspaceRoot: string, backupPath: string): boolean {
  const creatorFlowDir = join(workspaceRoot, '.creator-flow');

  if (!existsSync(backupPath)) {
    return false;
  }

  try {
    // Restore directories
    const dirsToRestore = ['prompts', 'guides', 'sources', 'skills'];
    for (const dir of dirsToRestore) {
      const sourceDir = join(backupPath, dir);
      if (existsSync(sourceDir)) {
        const targetDir = join(creatorFlowDir, dir);
        mergeDir(sourceDir, targetDir);
      }
    }

    // Restore manifest
    const manifestPath = join(backupPath, 'app-manifest.json');
    if (existsSync(manifestPath)) {
      copyFileSync(manifestPath, join(creatorFlowDir, 'app-manifest.json'));
    }

    return true;
  } catch {
    return false;
  }
}

// ============================================================
// App Installation
// ============================================================

/**
 * Download and cache an app package
 */
export async function downloadAppPackage(
  appId: string,
  version: string = 'latest',
  onProgress?: AppInstallProgressCallback
): Promise<{
  packageDir: string;
  actualVersion: string;
  skillDependencies: SkillDependencyDownload[] | null;
}> {
  onProgress?.({
    stage: 'downloading',
    percent: 0,
    message: `正在获取应用信息: ${appId}`,
    appId,
  });

  // Get download URL from API
  const downloadInfo = await getAppDownload(appId, version);

  if (!downloadInfo.download_url) {
    throw new Error('无法获取下载链接');
  }

  // Determine actual version from URL or use 'latest'
  const actualVersion =
    version === 'latest'
      ? extractVersionFromUrl(downloadInfo.download_url) || 'latest'
      : version;

  // Check if already cached
  if (isAppPackageCached(appId, actualVersion)) {
    onProgress?.({
      stage: 'downloading',
      percent: 30,
      message: '已使用缓存版本',
      appId,
    });
    return {
      packageDir: getCachedAppPackagePath(appId, actualVersion),
      actualVersion,
      skillDependencies: downloadInfo.skill_dependencies,
    };
  }

  // Create cache directory
  const cacheDir = getMarketplaceAppsDir();
  if (!existsSync(cacheDir)) {
    mkdirSync(cacheDir, { recursive: true });
  }
  const appCacheDir = join(cacheDir, appId);
  const versionDir = join(appCacheDir, actualVersion);
  const zipPath = join(appCacheDir, `${actualVersion}.zip`);

  if (!existsSync(appCacheDir)) {
    mkdirSync(appCacheDir, { recursive: true });
  }

  // Download the package
  onProgress?.({
    stage: 'downloading',
    percent: 5,
    message: `正在下载应用: ${appId}@${actualVersion}`,
    appId,
  });

  await downloadFile(downloadInfo.download_url, zipPath, (percent) => {
    onProgress?.({
      stage: 'downloading',
      percent: 5 + Math.round(percent * 0.2), // 5-25%
      message: `正在下载应用: ${percent}%`,
      appId,
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
    percent: 28,
    message: '正在解压应用...',
    appId,
  });

  if (!existsSync(versionDir)) {
    mkdirSync(versionDir, { recursive: true });
  }

  await extractZip(zipPath, versionDir);

  // Clean up zip file
  rmSync(zipPath);

  // Validate extracted content - check for manifest.json
  const manifestFile = join(versionDir, 'manifest.json');
  if (!existsSync(manifestFile)) {
    // Check if content is in a subdirectory
    const entries = readdirSync(versionDir, { withFileTypes: true });
    const subDir = entries.find((e) => e.isDirectory());
    if (subDir) {
      const subDirPath = join(versionDir, subDir.name);
      if (existsSync(join(subDirPath, 'manifest.json'))) {
        // Move contents up one level
        const subEntries = readdirSync(subDirPath);
        for (const entry of subEntries) {
          const src = join(subDirPath, entry);
          const dest = join(versionDir, entry);
          require('fs').renameSync(src, dest);
        }
        rmSync(subDirPath, { recursive: true });
      }
    }
  }

  // Final validation
  if (!existsSync(join(versionDir, 'manifest.json'))) {
    rmSync(versionDir, { recursive: true });
    throw new Error('无效的应用包：缺少 manifest.json 文件');
  }

  onProgress?.({
    stage: 'downloading',
    percent: 30,
    message: '应用下载完成',
    appId,
  });

  return {
    packageDir: versionDir,
    actualVersion,
    skillDependencies: downloadInfo.skill_dependencies,
  };
}

/**
 * Install skill dependencies for an app
 */
async function installSkillDependencies(
  workspaceRoot: string,
  dependencies: SkillDependencyDownload[],
  onProgress?: AppInstallProgressCallback
): Promise<Array<{ skillId: string; success: boolean; error?: string; skipped?: boolean }>> {
  const results: Array<{ skillId: string; success: boolean; error?: string; skipped?: boolean }> = [];
  const total = dependencies.length;

  // 获取全局技能目录，用于检查是否已存在相同技能
  const globalSkillsDir = getGlobalSkillsDir();

  onProgress?.({
    stage: 'installing-skills',
    percent: 35,
    message: `正在安装 ${total} 个技能依赖...`,
    totalSkills: total,
    installedSkills: 0,
  });

  for (let i = 0; i < dependencies.length; i++) {
    const dep = dependencies[i];
    if (!dep) continue;

    // 检查全局技能目录中是否已存在相同ID的技能
    const globalSkillPath = join(globalSkillsDir, dep.id);
    const globalSkillFile = join(globalSkillPath, 'SKILL.md');
    if (existsSync(globalSkillFile)) {
      // 全局技能已存在，跳过安装
      onProgress?.({
        stage: 'installing-skills',
        percent: 35 + Math.round(((i + 1) / total) * 50),
        message: `技能 ${dep.id} 已存在于全局技能中，跳过安装`,
        currentSkill: dep.id,
        totalSkills: total,
        installedSkills: i + 1,
      });

      results.push({
        skillId: dep.id,
        success: true,
        skipped: true,
      });
      continue;
    }

    onProgress?.({
      stage: 'installing-skills',
      percent: 35 + Math.round((i / total) * 50), // 35-85%
      message: `正在安装技能 ${i + 1}/${total}`,
      currentSkill: dep.id,
      totalSkills: total,
      installedSkills: i,
    });

    try {
      // Install skill with progress callback
      await installSkill(
        workspaceRoot,
        dep.id,
        dep.version,
        (skillProgress) => {
          onProgress?.({
            stage: 'installing-skills',
            percent: 35 + Math.round((i / total) * 50) + Math.round((skillProgress.percent / 100) * (50 / total)),
            message: `正在安装技能 ${i + 1}/${total}: ${skillProgress.message}`,
            currentSkill: dep.id,
            totalSkills: total,
            installedSkills: i,
            skillProgress: skillProgress.percent,
          });
        }
      );

      results.push({
        skillId: dep.id,
        success: true,
      });
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : '未知错误';
      results.push({
        skillId: dep.id,
        success: false,
        error: errorMessage,
      });
    }
  }

  // Check for failures and report
  const failedSkills = results.filter((r) => !r.success);
  if (failedSkills.length > 0) {
    onProgress?.({
      stage: 'installing-skills',
      percent: 85,
      message: `警告：${failedSkills.length} 个技能安装失败，应用可能无法正常工作`,
      totalSkills: total,
      installedSkills: results.filter(r => r.success).length,
    });
  } else {
    onProgress?.({
      stage: 'installing-skills',
      percent: 85,
      message: `所有技能安装完成 (${total}/${total})`,
      totalSkills: total,
      installedSkills: total,
    });
  }

  return results;
}

/**
 * Install an app to a workspace
 * This downloads the app, its skill dependencies, and copies everything to the workspace
 * @param workspaceRoot - Workspace root directory
 * @param appId - App ID to install
 * @param version - Version to install (default: 'latest')
 * @param options - Install options (force, skipBackup, onProgress)
 */
export async function installApp(
  workspaceRoot: string,
  appId: string,
  version: string = 'latest',
  options: AppInstallOptions | AppInstallProgressCallback = {}
): Promise<AppInstallResult> {
  // Support legacy signature: installApp(workspaceRoot, appId, version, onProgress)
  const opts: AppInstallOptions =
    typeof options === 'function' ? { onProgress: options } : options;
  const { force = false, merge = false, skipBackup = false, onProgress } = opts;

  try {
    // Check for existing installation
    const existingApp = checkInstalledApp(workspaceRoot);
    if (existingApp) {
      // merge mode allows installation without force flag
      if (!force && !merge) {
        return {
          success: false,
          appId,
          version,
          error: `工作区已安装应用 "${existingApp.name}" (${existingApp.version})。如需覆盖安装，请使用 force 选项。`,
        };
      }

      onProgress?.({
        stage: 'installing-app',
        percent: 2,
        message: merge ? '检测到已有应用，正在准备合并安装...' : '检测到已有应用，正在准备覆盖安装...',
        appId,
      });
    }

    // Backup existing data if force installing (not merge) and not skipping backup
    let backupPath: string | undefined;
    if (existingApp && force && !merge && !skipBackup) {
      onProgress?.({
        stage: 'installing-app',
        percent: 5,
        message: '正在备份现有数据...',
        appId,
      });
      backupPath = backupAppData(workspaceRoot) || undefined;
    }

    // Download app package
    const { packageDir, actualVersion, skillDependencies } = await downloadAppPackage(
      appId,
      version,
      onProgress
    );

    // Install skill dependencies if any
    let skillResults: Array<{ skillId: string; success: boolean; error?: string }> = [];
    if (skillDependencies && skillDependencies.length > 0) {
      skillResults = await installSkillDependencies(
        workspaceRoot,
        skillDependencies,
        onProgress
      );

      // Check if any skills failed
      const failedSkills = skillResults.filter((r) => !r.success);
      if (failedSkills.length > 0) {
        const errorMsg = `部分技能安装失败: ${failedSkills.map((s) => s.skillId).join(', ')}`;
        console.warn(errorMsg);

        // 如果所有技能都失败，返回失败状态
        if (failedSkills.length === skillDependencies.length) {
          return {
            success: false,
            appId,
            version: actualVersion,
            error: `所有技能依赖安装失败，应用无法正常使用`,
            skillResults,
          };
        }

        // 部分失败，继续安装但标记警告
        onProgress?.({
          stage: 'installing-app',
          percent: 85,
          message: `警告：${failedSkills.length}/${skillDependencies.length} 个技能安装失败`,
          appId,
        });
      }
    }

    // Copy app package to workspace
    onProgress?.({
      stage: 'installing-app',
      percent: 90,
      message: '正在复制应用到工作区...',
      appId,
    });

    const creatorFlowDir = join(workspaceRoot, '.creator-flow');
    if (!existsSync(creatorFlowDir)) {
      mkdirSync(creatorFlowDir, { recursive: true });
    }

    // Copy manifest.json to workspace (always overwrite)
    const sourceManifest = join(packageDir, 'manifest.json');
    const targetManifest = join(creatorFlowDir, 'app-manifest.json');
    copyFileSync(sourceManifest, targetManifest);

    // Merge app directories - always use mergeDir to preserve user customizations
    // Extended list includes labels, statuses, resources for merge install
    const appDirs = ['prompts', 'guides', 'sources', 'labels', 'statuses', 'resources', 'skills'];
    for (const dir of appDirs) {
      const sourceDir = join(packageDir, dir);
      if (existsSync(sourceDir)) {
        const targetDir = join(creatorFlowDir, dir);
        mergeDir(sourceDir, targetDir);
      }
    }

    onProgress?.({
      stage: 'complete',
      percent: 100,
      message: '应用安装完成',
      appId,
    });

    return {
      success: true,
      appId,
      version: actualVersion,
      appPath: creatorFlowDir,
      skillResults,
      backupPath,
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : '未知错误';

    onProgress?.({
      stage: 'error',
      percent: 0,
      message: errorMessage,
      appId,
      error: errorMessage,
    });

    return {
      success: false,
      appId,
      version,
      error: errorMessage,
    };
  }
}

/**
 * Install app from a local package path (for bundled apps)
 * @param workspaceRoot - Workspace root directory
 * @param packageDir - Path to the local app package
 * @param options - Install options (force, skipBackup, onProgress)
 */
export async function installAppFromLocal(
  workspaceRoot: string,
  packageDir: string,
  options: AppInstallOptions | AppInstallProgressCallback = {}
): Promise<AppInstallResult> {
  // Support legacy signature: installAppFromLocal(workspaceRoot, packageDir, onProgress)
  const opts: AppInstallOptions =
    typeof options === 'function' ? { onProgress: options } : options;
  const { force = false, merge = false, skipBackup = false, onProgress } = opts;

  try {
    // Read manifest to get app info and skill dependencies
    const manifestPath = join(packageDir, 'manifest.json');
    if (!existsSync(manifestPath)) {
      throw new Error('无效的应用包：缺少 manifest.json 文件');
    }

    const manifest = JSON.parse(readFileSync(manifestPath, 'utf-8'));
    const appId = manifest.id || 'local-app';

    // Check for existing installation
    const existingApp = checkInstalledApp(workspaceRoot);
    if (existingApp) {
      // merge mode allows installation without force flag
      if (!force && !merge) {
        return {
          success: false,
          appId,
          version: manifest.version || '1.0.0',
          error: `工作区已安装应用 "${existingApp.name}" (${existingApp.version})。如需覆盖安装，请使用 force 选项。`,
        };
      }

      onProgress?.({
        stage: 'installing-app',
        percent: 2,
        message: merge ? '检测到已有应用，正在准备合并安装...' : '检测到已有应用，正在准备覆盖安装...',
        appId,
      });
    }

    // Backup existing data if force installing (not merge) and not skipping backup
    let backupPath: string | undefined;
    if (existingApp && force && !merge && !skipBackup) {
      onProgress?.({
        stage: 'installing-app',
        percent: 5,
        message: '正在备份现有数据...',
        appId,
      });
      backupPath = backupAppData(workspaceRoot) || undefined;
    }

    onProgress?.({
      stage: 'installing-app',
      percent: 10,
      message: '正在解析应用配置...',
      appId,
    });

    // Parse skill dependencies from manifest
    const skillDependencies = manifest.skill_dependencies || [];
    let skillResults: Array<{ skillId: string; success: boolean; error?: string }> = [];

    if (skillDependencies.length > 0) {
      onProgress?.({
        stage: 'installing-app',
        percent: 20,
        message: `正在安装 ${skillDependencies.length} 个技能...`,
        appId,
        totalSkills: skillDependencies.length,
      });

      // Install skills from cloud
      const total = skillDependencies.length;
      for (let i = 0; i < skillDependencies.length; i++) {
        const skillId = skillDependencies[i];

        onProgress?.({
          stage: 'installing-app',
          percent: 20 + Math.round(((i + 0.5) / total) * 60), // 20-80%
          message: `正在安装技能: ${skillId} (${i + 1}/${total})`,
          appId,
          currentSkill: skillId,
          totalSkills: total,
          installedSkills: i,
        });

        try {
          const result = await installSkill(workspaceRoot, skillId, 'latest');
          skillResults.push({
            skillId,
            success: result.success,
            error: result.error,
          });
        } catch (error) {
          skillResults.push({
            skillId,
            success: false,
            error: error instanceof Error ? error.message : '未知错误',
          });
        }
      }
    }

    // Copy app package to workspace
    onProgress?.({
      stage: 'installing-app',
      percent: 85,
      message: '正在复制应用到工作区...',
      appId,
    });

    const creatorFlowDir = join(workspaceRoot, '.creator-flow');
    if (!existsSync(creatorFlowDir)) {
      mkdirSync(creatorFlowDir, { recursive: true });
    }

    // Copy manifest.json (always overwrite)
    const targetManifest = join(creatorFlowDir, 'app-manifest.json');
    copyFileSync(manifestPath, targetManifest);

    // Merge app directories - always use mergeDir to preserve user customizations
    // Extended list includes labels, statuses, resources for merge install
    const appDirs = ['prompts', 'guides', 'sources', 'labels', 'statuses', 'resources', 'skills'];
    for (const dir of appDirs) {
      const sourceDir = join(packageDir, dir);
      if (existsSync(sourceDir)) {
        const targetDir = join(creatorFlowDir, dir);
        mergeDir(sourceDir, targetDir);
      }
    }

    onProgress?.({
      stage: 'complete',
      percent: 100,
      message: '应用安装完成',
      appId,
    });

    return {
      success: true,
      appId,
      version: manifest.version || '1.0.0',
      appPath: creatorFlowDir,
      skillResults,
      backupPath,
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : '未知错误';

    onProgress?.({
      stage: 'error',
      percent: 0,
      message: errorMessage,
      error: errorMessage,
    });

    return {
      success: false,
      appId: 'unknown',
      version: 'unknown',
      error: errorMessage,
    };
  }
}

// ============================================================
// Utility Functions
// ============================================================

/**
 * Extract version from download URL
 */
function extractVersionFromUrl(url: string): string | null {
  const match = url.match(/\/(\d+\.\d+\.\d+)\.zip/);
  return match?.[1] ?? null;
}
