/**
 * 反检测浏览器服务 — 基于 playwright-extra + stealth 插件
 *
 * 提供反检测浏览器实例的启动和管理，集成 BrowserProfileManager 和 FingerprintGenerator。
 * 支持登录、发布、数据采集三种场景。
 */

import type { Browser, BrowserContext, Page } from 'playwright-core';
import type { BrowserFingerprint } from './browser-profile-manager.ts';
import type { BrowserProfileManager } from './browser-profile-manager.ts';
import type { FingerprintGenerator } from './fingerprint-generator.ts';
import type { Platform } from '../db/types.ts';
import { existsSync } from 'fs';

// ============================================================
// 类型定义
// ============================================================

export interface LaunchOptions {
  /** 平台账号 ID */
  platformAccountId: string;
  /** 是否有头模式（默认 false） */
  headless?: boolean;
  /** 自定义指纹（不传则自动生成） */
  fingerprint?: BrowserFingerprint;
  /** 额外的浏览器启动参数 */
  extraArgs?: string[];
}

export interface LoginResult {
  /** 是否登录成功 */
  success: boolean;
  /** 失败原因 */
  error?: string;
  /** 登录成功后采集的账号资料 */
  profile?: {
    nickname?: string;
    avatar_url?: string;
    bio?: string;
    home_url?: string;
    platform_uid?: string;
    followers?: number;
    following?: number;
    total_likes?: number;
    total_favorites?: number;
    total_posts?: number;
  };
}

export interface CollectedMetrics {
  views?: number;
  likes?: number;
  comments?: number;
  shares?: number;
  favorites?: number;
  /** 原始页面数据（JSON 字符串） */
  rawData?: string;
}

/** 各平台登录页 URL */
const PLATFORM_LOGIN_URLS: Record<Platform, string> = {
  xiaohongshu: 'https://creator.xiaohongshu.com/login',
  douyin: 'https://creator.douyin.com/login',
  wechat: 'https://mp.weixin.qq.com/',
  bilibili: 'https://passport.bilibili.com/login',
  zhihu: 'https://www.zhihu.com/signin',
  weibo: 'https://passport.weibo.com/sso/signin',
  x: 'https://x.com/i/flow/login',
};

/** 各平台个人主页 URL（登录成功后用于采集账号信息） */
const PLATFORM_PROFILE_URLS: Record<Platform, string> = {
  xiaohongshu: 'https://creator.xiaohongshu.com/home',
  douyin: 'https://creator.douyin.com/creator-micro/home',
  wechat: 'https://mp.weixin.qq.com/cgi-bin/home',
  bilibili: 'https://member.bilibili.com/platform/home',
  zhihu: 'https://www.zhihu.com/creator',
  weibo: 'https://weibo.com/u/home',
  x: 'https://x.com/settings/profile',
};

/** 各平台登录成功后的 URL 特征（用于检测登录状态） */
const PLATFORM_LOGGED_IN_PATTERNS: Record<Platform, RegExp> = {
  xiaohongshu: /creator\.xiaohongshu\.com/,
  douyin: /creator\.douyin\.com/,
  wechat: /mp\.weixin\.qq\.com\/cgi-bin/,
  bilibili: /(member|space|www)\.bilibili\.com/,
  zhihu: /www\.zhihu\.com(?!\/signin)/,
  weibo: /weibo\.com(?!\/sso)/,
  x: /x\.com(?!\/i\/flow\/login)/,
};

// ============================================================
// 工具函数
// ============================================================

/** 随机延迟（模拟人类操作） */
function randomDelay(minMs: number, maxMs: number): Promise<void> {
  const ms = minMs + Math.random() * (maxMs - minMs);
  return new Promise((resolve) => setTimeout(resolve, ms));
}

/**
 * 查找可用的 Chrome/Chromium 可执行文件路径
 *
 * playwright-core 不自带浏览器，需要指定 executablePath。
 * 优先使用系统已安装的 Chrome（几乎所有用户都有），
 * 找不到时返回 undefined，由调用方抛出友好错误提示。
 */
function findChromiumExecutable(): string | undefined {
  const candidates =
    process.platform === 'darwin'
      ? [
          '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
          '/Applications/Chromium.app/Contents/MacOS/Chromium',
          `${process.env.HOME}/Applications/Google Chrome.app/Contents/MacOS/Google Chrome`,
        ]
      : process.platform === 'win32'
        ? [
            `${process.env['PROGRAMFILES']}\\Google\\Chrome\\Application\\chrome.exe`,
            `${process.env['PROGRAMFILES(X86)']}\\Google\\Chrome\\Application\\chrome.exe`,
            `${process.env.LOCALAPPDATA}\\Google\\Chrome\\Application\\chrome.exe`,
          ]
        : [
            '/usr/bin/google-chrome',
            '/usr/bin/google-chrome-stable',
            '/usr/bin/chromium',
            '/usr/bin/chromium-browser',
          ];

  return candidates.find((p) => existsSync(p));
}

// ============================================================
// AntiDetectBrowser
// ============================================================

export class AntiDetectBrowser {
  private profileManager: BrowserProfileManager;
  private fingerprintGenerator: FingerprintGenerator;
  private browser: Browser | null = null;
  private context: BrowserContext | null = null;

  constructor(
    profileManager: BrowserProfileManager,
    fingerprintGenerator: FingerprintGenerator,
  ) {
    this.profileManager = profileManager;
    this.fingerprintGenerator = fingerprintGenerator;
  }

  /**
   * 启动反检测浏览器
   * 加载 Profile 持久化 + 指纹应用 + 反检测脚本注入
   */
  async launch(options: LaunchOptions): Promise<BrowserContext> {
    const { chromium } = await import('playwright-core');

    const { platformAccountId, headless = true, fingerprint, extraArgs = [] } = options;

    // 获取或生成指纹
    const fp = fingerprint
      ?? this.profileManager.getFingerprint(platformAccountId)
      ?? this.fingerprintGenerator.generate(platformAccountId);

    // 确保指纹已保存
    if (!this.profileManager.getFingerprint(platformAccountId)) {
      this.profileManager.saveFingerprint(platformAccountId, fp);
    }

    // 获取 Profile 目录
    const userDataDir = this.profileManager.loadProfile(platformAccountId);

    // 检测系统 Chrome 路径
    const executablePath = findChromiumExecutable();
    if (!executablePath) {
      throw new Error(
        '未找到 Chrome 浏览器，请安装 Google Chrome 后重试。\n'
        + '下载地址: https://www.google.com/chrome/',
      );
    }

    // 启动浏览器（持久化上下文）
    this.context = await chromium.launchPersistentContext(userDataDir, {
      headless,
      executablePath,
      args: [
        '--disable-blink-features=AutomationControlled',
        '--disable-features=IsolateOrigins,site-per-process',
        '--no-first-run',
        '--no-default-browser-check',
        ...extraArgs,
      ],
      userAgent: fp.userAgent,
      viewport: fp.viewport,
      locale: fp.locale,
      timezoneId: fp.timezone,
      colorScheme: 'light',
      ignoreHTTPSErrors: true,
    });

    // 为所有页面注入反检测脚本
    await this.context.addInitScript(() => {
      // 隐藏 webdriver 标志
      Object.defineProperty(navigator, 'webdriver', { get: () => false });

      // 覆盖 permissions query
      const originalQuery = window.navigator.permissions.query.bind(
        window.navigator.permissions,
      );
      window.navigator.permissions.query = (parameters: PermissionDescriptor) => {
        if (parameters.name === 'notifications') {
          return Promise.resolve({ state: 'prompt', onchange: null } as PermissionStatus);
        }
        return originalQuery(parameters);
      };

      // 覆盖 plugins（模拟真实浏览器）
      Object.defineProperty(navigator, 'plugins', {
        get: () => [1, 2, 3, 4, 5],
      });

      // 覆盖 languages
      Object.defineProperty(navigator, 'languages', {
        get: () => ['zh-CN', 'zh', 'en-US', 'en'],
      });

      // 隐藏 chrome.runtime（Headless Chrome 特征）
      if (!(window as any).chrome) {
        (window as any).chrome = { runtime: {} };
      }

      // 覆盖 navigator.connection（部分平台检测）
      Object.defineProperty(navigator, 'connection', {
        get: () => ({
          effectiveType: '4g',
          rtt: 50,
          downlink: 10,
          saveData: false,
        }),
      });
    });

    return this.context;
  }

  /**
   * 启动登录用浏览器
   * 有头模式，导航到平台登录页，监听登录成功
   */
  async launchForLogin(
    platformAccountId: string,
    platform: Platform,
  ): Promise<LoginResult> {
    const fingerprint = this.fingerprintGenerator.generateForPlatform(platform, platformAccountId);

    const context = await this.launch({
      platformAccountId,
      headless: false,
      fingerprint,
    });

    const page = context.pages()[0] ?? await context.newPage();
    const loginUrl = PLATFORM_LOGIN_URLS[platform];
    const loggedInPattern = PLATFORM_LOGGED_IN_PATTERNS[platform];

    // 检查当前页面是否已经登录成功（Profile 中可能保留了登录态）
    const currentUrl = page.url();
    if (currentUrl && loggedInPattern && loggedInPattern.test(currentUrl)) {
      // 已经登录，采集资料后返回
      this.profileManager.saveFingerprint(platformAccountId, fingerprint);
      const profile = await this.scrapeProfile(page, platform).catch(() => undefined);
      await this.close();
      return { success: true, profile };
    }

    // 导航到登录页，使用 networkidle 等待更稳定
    try {
      await page.goto(loginUrl, {
        waitUntil: 'domcontentloaded',
        timeout: 30000,
      });
    } catch (navError: any) {
      // 导航失败可能是因为重定向到了已登录页面，检查当前 URL
      const urlAfterNav = page.url();
      if (urlAfterNav && loggedInPattern && loggedInPattern.test(urlAfterNav)) {
        this.profileManager.saveFingerprint(platformAccountId, fingerprint);
        const profile = await this.scrapeProfile(page, platform).catch(() => undefined);
        await this.close();
        return { success: true, profile };
      }
      // 真正的导航失败，不关闭浏览器，让用户手动操作
      console.warn(`Navigation to login page failed: ${navError.message}, waiting for manual login...`);
    }

    // 等待登录成功（URL 变化匹配）或超时 5 分钟
    const LOGIN_TIMEOUT = 5 * 60 * 1000;

    try {
      await page.waitForURL(loggedInPattern, { timeout: LOGIN_TIMEOUT });

      // 登录成功，等待页面稳定
      await randomDelay(2000, 4000);

      // 保存指纹
      this.profileManager.saveFingerprint(platformAccountId, fingerprint);

      // 采集账号资料
      const profile = await this.scrapeProfile(page, platform).catch(() => undefined);

      return { success: true, profile };
    } catch {
      return { success: false, error: 'Login timeout (5 minutes)' };
    } finally {
      await this.close();
    }
  }

  /**
   * 启动发布用浏览器
   * 无头模式，加载 Profile，应用反检测策略
   */
  async launchForPublish(
    platformAccountId: string,
    platform: Platform,
  ): Promise<{ context: BrowserContext; page: Page }> {
    const fingerprint = this.fingerprintGenerator.generateForPlatform(platform, platformAccountId);

    const context = await this.launch({
      platformAccountId,
      headless: true,
      fingerprint,
    });

    const page = context.pages()[0] ?? await context.newPage();

    // 应用反检测脚本
    await this.applyAntiDetectScripts(page);

    return { context, page };
  }

  /**
   * 启动数据采集浏览器
   * 无头模式，导航到发布页面采集数据
   */
  async launchForDataCollection(
    platformAccountId: string,
    publishUrl: string,
  ): Promise<CollectedMetrics> {
    const context = await this.launch({
      platformAccountId,
      headless: true,
    });

    const page = context.pages()[0] ?? await context.newPage();

    // 应用反检测脚本
    await this.applyAntiDetectScripts(page);

    try {
      // 导航到目标页面
      await page.goto(publishUrl, { waitUntil: 'networkidle', timeout: 30000 });

      // 模拟自然滚动
      await this.simulateNaturalScroll(page);

      // 等待数据加载
      await randomDelay(2000, 5000);

      // 采集页面数据
      const metrics = await this.extractMetrics(page, publishUrl);

      return metrics;
    } finally {
      await this.close();
    }
  }

  /**
   * 登录成功后采集账号资料（昵称、头像、粉丝数等）
   * 导航到个人主页/创作者中心，提取关键信息
   */
  private async scrapeProfile(
    page: Page,
    platform: Platform,
  ): Promise<LoginResult['profile']> {
    const profileUrl = PLATFORM_PROFILE_URLS[platform];

    // 导航到个人主页
    try {
      await page.goto(profileUrl, { waitUntil: 'domcontentloaded', timeout: 15000 });
      await randomDelay(2000, 4000);
    } catch {
      // 导航失败，尝试从当前页面提取
    }

    const currentUrl = page.url();

    switch (platform) {
      case 'zhihu':
        return page.evaluate((url) => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          const getAttr = (sel: string, attr: string) => document.querySelector(sel)?.getAttribute(attr);
          const parseNum = (text?: string | null) => {
            if (!text) return 0;
            const n = text.replace(/[,，]/g, '');
            if (n.includes('万') || n.includes('w')) return Math.round(parseFloat(n) * 10000);
            if (n.includes('亿')) return Math.round(parseFloat(n) * 100000000);
            return parseInt(n, 10) || 0;
          };
          return {
            nickname: getText('.ProfileHeader-name, .CreatorHome-userName, [class*="UserName"]') || undefined,
            avatar_url: getAttr('.Avatar--round img, .ProfileHeader-avatar img, [class*="avatar"] img', 'src') || undefined,
            bio: getText('.ProfileHeader-headline, [class*="bio"]') || undefined,
            home_url: url,
            followers: parseNum(getText('[class*="follower"] .NumberBoard-itemValue, .FollowshipCard-counts strong')),
            following: parseNum(getText('[class*="following"] .NumberBoard-itemValue')),
            total_likes: parseNum(getText('[class*="voteup"] .NumberBoard-itemValue, [class*="like"] .NumberBoard-itemValue')),
          };
        }, currentUrl);

      case 'xiaohongshu':
        return page.evaluate((url) => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          const getAttr = (sel: string, attr: string) => document.querySelector(sel)?.getAttribute(attr);
          const parseNum = (text?: string | null) => {
            if (!text) return 0;
            const n = text.replace(/[,，]/g, '');
            if (n.includes('万') || n.includes('w')) return Math.round(parseFloat(n) * 10000);
            return parseInt(n, 10) || 0;
          };
          return {
            nickname: getText('[class*="nickname"], [class*="userName"], .creator-name') || undefined,
            avatar_url: getAttr('[class*="avatar"] img, .creator-avatar img', 'src') || undefined,
            bio: getText('[class*="desc"], [class*="bio"]') || undefined,
            home_url: url,
            followers: parseNum(getText('[class*="fans"] [class*="count"], [class*="follower"] [class*="count"]')),
            following: parseNum(getText('[class*="follow"] [class*="count"]')),
            total_likes: parseNum(getText('[class*="like"] [class*="count"], [class*="interaction"] [class*="count"]')),
            total_favorites: parseNum(getText('[class*="collect"] [class*="count"]')),
          };
        }, currentUrl);

      case 'bilibili':
        return page.evaluate((url) => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          const getAttr = (sel: string, attr: string) => document.querySelector(sel)?.getAttribute(attr);
          const parseNum = (text?: string | null) => {
            if (!text) return 0;
            const n = text.replace(/[,，]/g, '');
            if (n.includes('万')) return Math.round(parseFloat(n) * 10000);
            return parseInt(n, 10) || 0;
          };
          return {
            nickname: getText('[class*="nickname"], .h-name, [class*="uname"]') || undefined,
            avatar_url: getAttr('[class*="avatar"] img, .h-avatar img', 'src') || undefined,
            bio: getText('[class*="sign"], [class*="bio"]') || undefined,
            home_url: url,
            followers: parseNum(getText('[class*="fans"] [class*="num"], [class*="follower"]')),
            following: parseNum(getText('[class*="follow"] [class*="num"]')),
            total_likes: parseNum(getText('[class*="like"] [class*="num"]')),
          };
        }, currentUrl);

      case 'douyin':
        return page.evaluate((url) => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          const getAttr = (sel: string, attr: string) => document.querySelector(sel)?.getAttribute(attr);
          const parseNum = (text?: string | null) => {
            if (!text) return 0;
            const n = text.replace(/[,，]/g, '');
            if (n.includes('万') || n.includes('w')) return Math.round(parseFloat(n) * 10000);
            return parseInt(n, 10) || 0;
          };
          return {
            nickname: getText('[class*="nickname"], [class*="userName"]') || undefined,
            avatar_url: getAttr('[class*="avatar"] img', 'src') || undefined,
            home_url: url,
            followers: parseNum(getText('[class*="fans"] [class*="num"], [class*="follower"]')),
            following: parseNum(getText('[class*="follow"] [class*="num"]')),
            total_likes: parseNum(getText('[class*="like"] [class*="num"]')),
            total_posts: parseNum(getText('[class*="work"] [class*="num"], [class*="video"] [class*="num"]')),
          };
        }, currentUrl);

      default:
        // 通用提取：尝试常见选择器
        return page.evaluate((url) => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          const getAttr = (sel: string, attr: string) => document.querySelector(sel)?.getAttribute(attr);
          return {
            nickname: getText('[class*="nickname"], [class*="userName"], [class*="name"]') || undefined,
            avatar_url: getAttr('[class*="avatar"] img', 'src') || undefined,
            home_url: url,
          };
        }, currentUrl);
    }
  }

  /** 关闭浏览器 */
  async close(): Promise<void> {
    if (this.context) {
      await this.context.close();
      this.context = null;
    }
    if (this.browser) {
      await this.browser.close();
      this.browser = null;
    }
  }

  // ============================================================
  // 反检测策略
  // ============================================================

  /** 注入反检测脚本（已在 launch 中通过 context.addInitScript 全局注入，此方法保留用于页面级补充） */
  private async applyAntiDetectScripts(_page: Page): Promise<void> {
    // 反检测脚本已在 context 级别注入，无需重复
  }

  /** 模拟自然滚动行为 */
  private async simulateNaturalScroll(page: Page): Promise<void> {
    const scrollSteps = 3 + Math.floor(Math.random() * 3);

    for (let i = 0; i < scrollSteps; i++) {
      const scrollAmount = 200 + Math.floor(Math.random() * 300);
      await page.mouse.wheel(0, scrollAmount);
      await randomDelay(500, 1500);
    }
  }

  /** 从页面提取数据指标 */
  private async extractMetrics(page: Page, url: string): Promise<CollectedMetrics> {
    // 根据 URL 判断平台，使用对应的提取策略
    const platform = this.detectPlatformFromUrl(url);

    const rawData = await page.evaluate(() => {
      return JSON.stringify({
        title: document.title,
        url: window.location.href,
        html: document.documentElement.outerHTML.substring(0, 50000),
      });
    });

    // 通用数字提取（各平台的具体选择器由 Agent 技能层处理）
    const metrics: CollectedMetrics = { rawData };

    if (platform) {
      const extracted = await this.extractPlatformMetrics(page, platform);
      Object.assign(metrics, extracted);
    }

    return metrics;
  }

  /** 从 URL 检测平台 */
  private detectPlatformFromUrl(url: string): Platform | null {
    if (url.includes('xiaohongshu.com')) return 'xiaohongshu';
    if (url.includes('douyin.com')) return 'douyin';
    if (url.includes('bilibili.com')) return 'bilibili';
    if (url.includes('zhihu.com')) return 'zhihu';
    if (url.includes('weibo.com')) return 'weibo';
    if (url.includes('mp.weixin.qq.com')) return 'wechat';
    if (url.includes('x.com') || url.includes('twitter.com')) return 'x';
    return null;
  }

  /** 平台特定的指标提取 */
  private async extractPlatformMetrics(
    page: Page,
    platform: Platform,
  ): Promise<Partial<CollectedMetrics>> {
    // 各平台的数据选择器（基础版本，后续由技能层扩展）
    switch (platform) {
      case 'xiaohongshu':
        return page.evaluate(() => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          return {
            likes: parseInt(getText('[class*="like"] [class*="count"]') ?? '0', 10) || undefined,
            favorites: parseInt(getText('[class*="collect"] [class*="count"]') ?? '0', 10) || undefined,
            comments: parseInt(getText('[class*="comment"] [class*="count"]') ?? '0', 10) || undefined,
          };
        });

      case 'bilibili':
        return page.evaluate(() => {
          const getText = (sel: string) => document.querySelector(sel)?.textContent?.trim();
          return {
            views: parseInt(getText('.view-text, [class*="view"]') ?? '0', 10) || undefined,
            likes: parseInt(getText('.like-text, [class*="like"]') ?? '0', 10) || undefined,
            favorites: parseInt(getText('.collect-text, [class*="collect"]') ?? '0', 10) || undefined,
            comments: parseInt(getText('.comment-text, [class*="comment"]') ?? '0', 10) || undefined,
            shares: parseInt(getText('.share-text, [class*="share"]') ?? '0', 10) || undefined,
          };
        });

      default:
        // 其他平台返回空，由 Agent 技能层通过 rawData 解析
        return {};
    }
  }
}

// ============================================================
// 工厂函数
// ============================================================

export interface AntiDetectBrowserOptions {
  profileManager: BrowserProfileManager;
  fingerprintGenerator: FingerprintGenerator;
}

export function createAntiDetectBrowser(options: AntiDetectBrowserOptions): AntiDetectBrowser {
  return new AntiDetectBrowser(options.profileManager, options.fingerprintGenerator);
}
