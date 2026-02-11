/**
 * social-publisher 模块导出
 */

// 类型
export type {
  SocialPlatform,
  BrowserFingerprint,
  PublishContent,
  PublishResult,
  MetricsResult,
  CommentResult,
  PlatformAdapter,
} from './types.ts';

// 浏览器引擎
export { BrowserManager } from './browser/browser-manager.ts';
export { FingerprintManager } from './browser/fingerprint.ts';
export { HumanBehavior } from './browser/human-behavior.ts';
export { getStealthPlugin, getAntiDetectScripts } from './browser/stealth.ts';

// 平台适配器
export { BaseAdapter } from './adapters/base-adapter.ts';
export { XiaohongshuAdapter } from './adapters/xiaohongshu.ts';
export { getAdapter, getSupportedPlatforms } from './adapters/registry.ts';
