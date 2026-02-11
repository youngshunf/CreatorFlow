/**
 * social-publisher — 浏览器生命周期管理
 *
 * 基于 Playwright 实现反检测浏览器的启动、Profile 持久化和关闭。
 * 支持有头模式（用户登录）和无头模式（自动化操作）。
 *
 * Profile 存储路径: {workspacePath}/.sprouty-ai/browser-profiles/{platform}/
 */

import { join } from 'path';
import { mkdir } from 'fs/promises';
import type { Browser, BrowserContext, Page } from 'playwright';
import type { BrowserFingerprint, SocialPlatform, PlatformAdapter } from '../types.ts';
import { FingerprintManager } from './fingerprint.ts';
import { getAntiDetectScripts } from './stealth.ts';

// ============================================================
// 常量
// ============================================================

/** 浏览器 Profile 根目录名 */
const PROFILES_DIR = '.sprouty-ai/browser-profiles';

/** 登录检测轮询间隔（毫秒） */
const LOGIN_CHECK_INTERVAL = 2000;

/** 登录超时时间（毫秒）— 5 分钟 */
const LOGIN_TIMEOUT = 5 * 60 * 1000;

// ============================================================
// BrowserManager
// ============================================================

export class BrowserManager {
  private browser: Browser | null = null;
  private context: BrowserContext | null = null;
  private readonly fingerprintManager = new FingerprintManager();
  private readonly workspacePath: string;

  constructor(workspacePath: string) {
    this.workspacePath = workspacePath;
  }

  /**
   * 获取平台的 Profile 存储路径
   */
  getProfilePath(platform: SocialPlatform): string {
    return join(this.workspacePath, PROFILES_DIR, platform);
  }

  /**
   * 启动有头浏览器用于用户登录
   *
   * 流程：
   * 1. 启动带 UI 的 Chromium 浏览器
   * 2. 导航到平台登录页
   * 3. 轮询检测登录状态
   * 4. 登录成功后保存 Profile 和指纹
   *
   * @param adapter - 平台适配器
   * @param fingerprint - 可选的指纹配置，不传则自动生成/加载
   * @returns 登录成功的页面实例
   */
  async launchForLogin(
    adapter: PlatformAdapter,
    fingerprint?: BrowserFingerprint,
  ): Promise<{ page: Page; context: BrowserContext }> {
    const platform = adapter.platform;
    const profilePath = this.getProfilePath(platform);
    await mkdir(profilePath, { recursive: true });

    // 解析指纹：优先使用传入的 > 已保存的 > 新生成的
    const fp = fingerprint
      ?? await this.fingerprintManager.loadFingerprint(profilePath)
      ?? this.fingerprintManager.generateFingerprint(platform);

    // 动态导入 playwright（避免硬依赖）
    const { chromium } = await import('playwright');

    // 启动有头浏览器，使用 userDataDir 持久化 Profile
    this.context = await chromium.launchPersistentContext(profilePath, {
      headless: false,
      viewport: fp.viewport,
      userAgent: fp.userAgent,
      locale: fp.locale,
      timezoneId: fp.timezone,
      args: [
        '--disable-blink-features=AutomationControlled',
        '--disable-features=IsolateOrigins,site-per-process',
        '--no-first-run',
        '--no-default-browser-check',
      ],
      ignoreDefaultArgs: ['--enable-automation'],
    });

    // 应用指纹和反检测脚本
    await this.fingerprintManager.applyFingerprint(this.context, fp);
    await this.context.addInitScript(getAntiDetectScripts());

    // 保存指纹配置
    await this.fingerprintManager.saveFingerprint(profilePath, fp);

    // 获取页面并导航到登录页
    const page = this.context.pages()[0] ?? await this.context.newPage();
    await page.goto(adapter.getLoginUrl(), { waitUntil: 'domcontentloaded' });

    // 轮询检测登录状态
    const loginSuccess = await this.waitForLogin(page, adapter);
    if (!loginSuccess) {
      throw new Error(`[social-publisher] ${platform} 登录超时，请在 ${LOGIN_TIMEOUT / 1000} 秒内完成登录`);
    }

    return { page, context: this.context };
  }

  /**
   * 启动无头浏览器用于自动化操作
   *
   * 复用已保存的 Profile（Cookie、LocalStorage 等），
   * 以无头模式运行，适用于发布、采集等自动化任务。
   *
   * @param platform - 目标平台
   * @param fingerprint - 可选的指纹配置
   * @returns 浏览器上下文
   */
  async launchHeadless(
    platform: SocialPlatform,
    fingerprint?: BrowserFingerprint,
  ): Promise<BrowserContext> {
    const profilePath = this.getProfilePath(platform);

    // 解析指纹
    const fp = fingerprint
      ?? await this.fingerprintManager.loadFingerprint(profilePath)
      ?? this.fingerprintManager.generateFingerprint(platform);

    const { chromium } = await import('playwright');

    // 启动无头浏览器，复用 Profile
    this.context = await chromium.launchPersistentContext(profilePath, {
      headless: true,
      viewport: fp.viewport,
      userAgent: fp.userAgent,
      locale: fp.locale,
      timezoneId: fp.timezone,
      args: [
        '--disable-blink-features=AutomationControlled',
        '--disable-features=IsolateOrigins,site-per-process',
        '--no-first-run',
        '--no-default-browser-check',
      ],
      ignoreDefaultArgs: ['--enable-automation'],
    });

    // 应用指纹和反检测脚本
    await this.fingerprintManager.applyFingerprint(this.context, fp);
    await this.context.addInitScript(getAntiDetectScripts());

    return this.context;
  }

  /**
   * 关闭浏览器和上下文，释放资源
   */
  async close(): Promise<void> {
    if (this.context) {
      await this.context.close().catch(() => {});
      this.context = null;
    }
    if (this.browser) {
      await this.browser.close().catch(() => {});
      this.browser = null;
    }
  }

  // ============================================================
  // 内部方法
  // ============================================================

  /**
   * 轮询等待用户完成登录
   *
   * 每隔 LOGIN_CHECK_INTERVAL 调用适配器的 detectLoginSuccess 方法，
   * 直到检测到登录成功或超时。
   */
  private async waitForLogin(page: Page, adapter: PlatformAdapter): Promise<boolean> {
    const startTime = Date.now();

    while (Date.now() - startTime < LOGIN_TIMEOUT) {
      try {
        const success = await adapter.detectLoginSuccess(page);
        if (success) return true;
      } catch {
        // 页面可能正在跳转，忽略检测错误
      }
      await sleep(LOGIN_CHECK_INTERVAL);
    }

    return false;
  }
}

// ============================================================
// 工具函数
// ============================================================

/** 异步等待指定毫秒数 */
function sleep(ms: number): Promise<void> {
  return new Promise(resolve => setTimeout(resolve, ms));
}
