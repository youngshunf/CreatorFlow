/**
 * 浏览器指纹生成器 — 生成随机但一致的浏览器指纹配置
 *
 * 基于 platformAccountId 的确定性随机，确保同一账号始终生成相同指纹。
 * 与 BrowserProfileManager.saveFingerprint() 集成使用。
 */

import type { BrowserFingerprint } from './browser-profile-manager.ts';

// ============================================================
// 指纹数据池
// ============================================================

/** 常见 macOS Chrome User-Agent */
const USER_AGENTS = [
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36',
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36',
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36',
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36',
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36',
];

/** 常见屏幕分辨率 */
const VIEWPORTS = [
  { width: 1920, height: 1080 },
  { width: 2560, height: 1440 },
  { width: 1680, height: 1050 },
  { width: 1440, height: 900 },
  { width: 1536, height: 864 },
  { width: 1366, height: 768 },
];

/** WebGL 渲染器 */
const WEBGL_RENDERERS = [
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M1, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M1 Pro, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M2, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M2 Pro, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Apple)', renderer: 'ANGLE (Apple, Apple M3, OpenGL 4.1)' },
  { vendor: 'Google Inc. (Intel)', renderer: 'ANGLE (Intel, Intel(R) UHD Graphics 630, OpenGL 4.1)' },
];

/** 硬件并发数 */
const HARDWARE_CONCURRENCIES = [4, 8, 10, 12, 16];

/** 设备内存 (GB) */
const DEVICE_MEMORIES = [4, 8, 16, 32];

/** 时区 */
const TIMEZONES = [
  'Asia/Shanghai',
  'Asia/Shanghai',  // 权重加倍
  'Asia/Shanghai',  // 权重加倍
  'Asia/Chongqing',
  'Asia/Hong_Kong',
];

// ============================================================
// 确定性随机数生成器（基于 seed 的简单哈希）
// ============================================================

/** 基于字符串 seed 的确定性伪随机数生成器 */
class SeededRandom {
  private state: number;

  constructor(seed: string) {
    // 简单的字符串哈希
    this.state = 0;
    for (let i = 0; i < seed.length; i++) {
      this.state = ((this.state << 5) - this.state + seed.charCodeAt(i)) | 0;
    }
    // 确保非零
    if (this.state === 0) this.state = 1;
  }

  /** 返回 [0, 1) 的伪随机数 */
  next(): number {
    // xorshift32
    this.state ^= this.state << 13;
    this.state ^= this.state >> 17;
    this.state ^= this.state << 5;
    return (this.state >>> 0) / 4294967296;
  }

  /** 从数组中随机选择一个元素 */
  pick<T>(arr: T[]): T {
    return arr[Math.floor(this.next() * arr.length)];
  }
}

// ============================================================
// FingerprintGenerator
// ============================================================

/** 指纹生成器 */
export class FingerprintGenerator {
  /**
   * 生成随机但一致的浏览器指纹配置
   * @param platformAccountId - 平台账号 ID，作为确定性随机的 seed
   */
  generate(platformAccountId: string): BrowserFingerprint {
    const rng = new SeededRandom(platformAccountId);
    const now = new Date().toISOString();

    const webgl = rng.pick(WEBGL_RENDERERS);

    return {
      userAgent: rng.pick(USER_AGENTS),
      viewport: rng.pick(VIEWPORTS),
      locale: 'zh-CN',
      timezone: rng.pick(TIMEZONES),
      platform: 'MacIntel',
      webglVendor: webgl.vendor,
      webglRenderer: webgl.renderer,
      hardwareConcurrency: rng.pick(HARDWARE_CONCURRENCIES),
      deviceMemory: rng.pick(DEVICE_MEMORIES),
      createdAt: now,
      updatedAt: now,
    };
  }

  /**
   * 根据平台特点生成指纹
   * 不同平台可能需要不同的浏览器特征（如移动端 UA）
   * @param platform - 平台标识
   * @param platformAccountId - 平台账号 ID
   */
  generateForPlatform(platform: string, platformAccountId: string): BrowserFingerprint {
    const fingerprint = this.generate(platformAccountId);

    // 根据平台调整指纹
    switch (platform) {
      case 'douyin':
      case 'xiaohongshu':
        // 抖音和小红书桌面端使用较新的 Chrome 版本
        fingerprint.userAgent = fingerprint.userAgent.replace(
          /Chrome\/\d+/,
          'Chrome/125'
        );
        break;

      case 'bilibili':
        // B站对分辨率不敏感，保持默认
        break;

      case 'weibo':
      case 'zhihu':
        // 微博和知乎使用标准配置
        break;

      case 'wechat':
        // 微信公众号后台使用标准桌面配置
        break;

      default:
        break;
    }

    return fingerprint;
  }
}

/** 工厂函数 */
export function createFingerprintGenerator(): FingerprintGenerator {
  return new FingerprintGenerator();
}
