/**
 * 浏览器 Profile 管理器 — 管理平台账号的 Chromium Profile 和指纹配置
 *
 * Profile 目录结构：
 * {workspace}/.sprouty-ai/browser-profiles/{platformAccountId}/
 * ├── Default/          # Chromium userDataDir（供 Playwright 使用）
 * └── fingerprint.json  # 指纹配置
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  readdirSync,
  rmSync,
} from 'fs';
import { join } from 'path';
import { atomicWriteFileSync } from '../utils/files.ts';

/** 浏览器指纹配置 */
export interface BrowserFingerprint {
  userAgent: string;
  viewport: { width: number; height: number };
  locale: string;
  timezone: string;
  platform: string;
  webglVendor?: string;
  webglRenderer?: string;
  hardwareConcurrency?: number;
  deviceMemory?: number;
  extraHeaders?: Record<string, string>;
  createdAt: string;
  updatedAt: string;
}

/** BrowserProfileManager 配置 */
export interface BrowserProfileManagerOptions {
  /** 工作区根路径 */
  workspaceRootPath: string;
}

const WORKSPACE_DATA_DIR = '.sprouty-ai';
const BROWSER_PROFILES_DIR = 'browser-profiles';
const CHROMIUM_PROFILE_DIR = 'Default';
const FINGERPRINT_FILE = 'fingerprint.json';

/** 浏览器 Profile 管理器 */
export class BrowserProfileManager {
  private basePath: string;

  constructor(options: BrowserProfileManagerOptions) {
    this.basePath = join(options.workspaceRootPath, WORKSPACE_DATA_DIR, BROWSER_PROFILES_DIR);
  }

  /** 获取 Profile 根目录路径 */
  getProfilePath(platformAccountId: string): string {
    return join(this.basePath, platformAccountId, CHROMIUM_PROFILE_DIR);
  }

  /** 确保 Profile 目录存在 */
  ensureProfileDir(platformAccountId: string): string {
    const profilePath = this.getProfilePath(platformAccountId);
    if (!existsSync(profilePath)) {
      mkdirSync(profilePath, { recursive: true });
    }
    return profilePath;
  }

  /**
   * 加载 Profile 目录路径（供 Playwright 使用）
   * 返回 userDataDir 路径（Profile 根目录的父目录），如果不存在则自动创建
   */
  loadProfile(platformAccountId: string): string {
    this.ensureProfileDir(platformAccountId);
    // Playwright 的 userDataDir 是 Default/ 的父目录
    return join(this.basePath, platformAccountId);
  }

  /**
   * Profile 保存（Playwright 关闭浏览器时自动保存到 userDataDir）
   * 此方法仅确保目录存在，实际保存由 Playwright 完成
   */
  saveProfile(platformAccountId: string): void {
    this.ensureProfileDir(platformAccountId);
  }

  /** 删除 Profile 目录 */
  deleteProfile(platformAccountId: string): boolean {
    const profileRoot = join(this.basePath, platformAccountId);
    if (!existsSync(profileRoot)) return false;

    try {
      rmSync(profileRoot, { recursive: true });
      return true;
    } catch {
      return false;
    }
  }

  /** 检查 Profile 是否存在 */
  profileExists(platformAccountId: string): boolean {
    return existsSync(this.getProfilePath(platformAccountId));
  }

  /** 读取 fingerprint.json */
  getFingerprint(platformAccountId: string): BrowserFingerprint | null {
    const fingerprintPath = join(this.basePath, platformAccountId, FINGERPRINT_FILE);
    if (!existsSync(fingerprintPath)) return null;

    try {
      return JSON.parse(readFileSync(fingerprintPath, 'utf-8')) as BrowserFingerprint;
    } catch {
      return null;
    }
  }

  /** 保存 fingerprint.json */
  saveFingerprint(platformAccountId: string, config: BrowserFingerprint): void {
    const profileRoot = join(this.basePath, platformAccountId);
    if (!existsSync(profileRoot)) {
      mkdirSync(profileRoot, { recursive: true });
    }

    const fingerprintPath = join(profileRoot, FINGERPRINT_FILE);
    atomicWriteFileSync(fingerprintPath, JSON.stringify(config, null, 2));
  }

  /** 列出所有已存在的 Profile ID */
  listProfiles(): string[] {
    if (!existsSync(this.basePath)) return [];

    try {
      return readdirSync(this.basePath, { withFileTypes: true })
        .filter((d) => d.isDirectory())
        .map((d) => d.name);
    } catch {
      return [];
    }
  }
}

/** 工厂函数 */
export function createBrowserProfileManager(workspaceRootPath: string): BrowserProfileManager {
  return new BrowserProfileManager({ workspaceRootPath });
}
