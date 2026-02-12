/**
 * social-publisher — 浏览器指纹管理
 *
 * 为每个平台生成合理的浏览器指纹配置，并提供持久化存储和应用能力。
 * 指纹一旦生成会保存到 Profile 目录，后续复用以保持一致性。
 */

import { readFile, writeFile, mkdir } from 'fs/promises';
import { join, dirname } from 'path';
import type { BrowserContext } from 'playwright';
import type { BrowserFingerprint, SocialPlatform } from '../types.ts';

// ============================================================
// 常量：各平台默认视口尺寸
// ============================================================

/** 各平台推荐的视口尺寸（模拟移动端或桌面端） */
const PLATFORM_VIEWPORTS: Record<SocialPlatform, { width: number; height: number }> = {
  xiaohongshu: { width: 390, height: 844 },   // iPhone 14 Pro
  douyin: { width: 375, height: 812 },         // iPhone X
  wechat: { width: 390, height: 844 },         // iPhone 14 Pro
  bilibili: { width: 1920, height: 1080 },     // 桌面端
  zhihu: { width: 1440, height: 900 },         // 桌面端
  weibo: { width: 1440, height: 900 },         // 桌面端
  x: { width: 1440, height: 900 },             // 桌面端
  toutiao: { width: 1440, height: 900 },       // 桌面端
  sina: { width: 1440, height: 900 },           // 桌面端
  sohu: { width: 1440, height: 900 },           // 桌面端
};

/** 常见 Chrome 版本号（保持更新） */
const CHROME_VERSIONS = [
  '120.0.6099.109',
  '121.0.6167.85',
  '122.0.6261.94',
  '123.0.6312.58',
  '124.0.6367.91',
  '125.0.6422.76',
];

/** 常见 WebGL 渲染器配置 */
const WEBGL_CONFIGS = [
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M1, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M2, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M3, OpenGL 4.1)' },
  { vendor: 'Google Inc. (NVIDIA)', renderer: 'ANGLE (NVIDIA, NVIDIA GeForce RTX 3060, OpenGL 4.5)' },
  { vendor: 'Google Inc. (NVIDIA)', renderer: 'ANGLE (NVIDIA, NVIDIA GeForce RTX 4070, OpenGL 4.5)' },
  { vendor: 'Google Inc. (Intel)', renderer: 'ANGLE (Intel, Intel(R) UHD Graphics 630, OpenGL 4.5)' },
];

/** 指纹配置文件名 */
const FINGERPRINT_FILENAME = 'fingerprint.json';

// ============================================================
// FingerprintManager
// ============================================================

export class FingerprintManager {
  /**
   * 为指定平台生成合理的浏览器指纹
   */
  generateFingerprint(platform: SocialPlatform): BrowserFingerprint {
    const viewport = PLATFORM_VIEWPORTS[platform];
    const chromeVersion = randomPick(CHROME_VERSIONS);
    const webgl = randomPick(WEBGL_CONFIGS);
    const isMobile = viewport.width < 500;

    const platformStr = isMobile ? 'Linux armv81' : 'MacIntel';
    const mobileUA = isMobile
      ? ' Mobile' : '';

    return {
      userAgent: `Mozilla/5.0 (${isMobile ? 'Linux; Android 14; Pixel 8 Pro' : 'Macintosh; Intel Mac OS X 10_15_7'}) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/${chromeVersion}${mobileUA} Safari/537.36`,
      viewport,
      locale: 'zh-CN',
      timezone: 'Asia/Shanghai',
      webglVendor: webgl.vendor,
      webglRenderer: webgl.renderer,
      platform: platformStr,
      hardwareConcurrency: randomInt(4, 16),
      deviceMemory: randomPick([4, 8, 8, 16]),
    };
  }

  /**
   * 从 Profile 目录加载已保存的指纹配置
   * @returns 指纹配置，文件不存在时返回 null
   */
  async loadFingerprint(profilePath: string): Promise<BrowserFingerprint | null> {
    try {
      const filePath = join(profilePath, FINGERPRINT_FILENAME);
      const raw = await readFile(filePath, 'utf-8');
      return JSON.parse(raw) as BrowserFingerprint;
    } catch {
      return null;
    }
  }

  /**
   * 保存指纹配置到 Profile 目录
   */
  async saveFingerprint(profilePath: string, fingerprint: BrowserFingerprint): Promise<void> {
    await mkdir(profilePath, { recursive: true });
    const filePath = join(profilePath, FINGERPRINT_FILENAME);
    await writeFile(filePath, JSON.stringify(fingerprint, null, 2), 'utf-8');
  }

  /**
   * 将指纹配置应用到浏览器上下文
   *
   * 注入 WebGL 参数覆盖、hardwareConcurrency、deviceMemory 等属性，
   * 使浏览器环境与指纹配置一致。
   */
  async applyFingerprint(context: BrowserContext, fingerprint: BrowserFingerprint): Promise<void> {
    // 注入指纹覆盖脚本（在每个新页面加载前执行）
    await context.addInitScript(`
      // 覆盖 WebGL 调试信息
      const getParameterProxyHandler = {
        apply(target, thisArg, args) {
          const param = args[0];
          // UNMASKED_VENDOR_WEBGL
          if (param === 0x9245) return '${fingerprint.webglVendor}';
          // UNMASKED_RENDERER_WEBGL
          if (param === 0x9246) return '${fingerprint.webglRenderer}';
          return Reflect.apply(target, thisArg, args);
        },
      };

      // 拦截 WebGL getParameter
      const origGetParam = WebGLRenderingContext.prototype.getParameter;
      WebGLRenderingContext.prototype.getParameter = new Proxy(origGetParam, getParameterProxyHandler);

      // 拦截 WebGL2 getParameter
      if (typeof WebGL2RenderingContext !== 'undefined') {
        const origGetParam2 = WebGL2RenderingContext.prototype.getParameter;
        WebGL2RenderingContext.prototype.getParameter = new Proxy(origGetParam2, getParameterProxyHandler);
      }

      // 覆盖 hardwareConcurrency
      Object.defineProperty(navigator, 'hardwareConcurrency', {
        get: () => ${fingerprint.hardwareConcurrency},
      });

      // 覆盖 deviceMemory
      Object.defineProperty(navigator, 'deviceMemory', {
        get: () => ${fingerprint.deviceMemory},
      });

      // 覆盖 platform
      Object.defineProperty(navigator, 'platform', {
        get: () => '${fingerprint.platform}',
      });
    `);
  }
}

// ============================================================
// 工具函数
// ============================================================

/** 从数组中随机选取一个元素 */
function randomPick<T>(arr: T[]): T {
  return arr[Math.floor(Math.random() * arr.length)]!;
}

/** 生成 [min, max] 范围内的随机整数 */
function randomInt(min: number, max: number): number {
  return Math.floor(Math.random() * (max - min + 1)) + min;
}
