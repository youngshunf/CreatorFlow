/**
 * 平台适配器注册表
 *
 * 根据平台名称获取对应的适配器实例
 */

import type { PlatformAdapter, SocialPlatform } from '../types.ts';
import { XiaohongshuAdapter } from './xiaohongshu.ts';

/** 适配器工厂函数映射 */
const adapterFactories: Record<string, () => PlatformAdapter> = {
  xiaohongshu: () => new XiaohongshuAdapter(),
  // P1: douyin, wechat, bilibili
  // P2: zhihu, weibo, toutiao, sina, sohu
};

/**
 * 获取指定平台的适配器
 * @returns 适配器实例，不支持的平台返回 null
 */
export function getAdapter(platform: string): PlatformAdapter | null {
  const factory = adapterFactories[platform];
  return factory ? factory() : null;
}

/**
 * 获取所有已支持的平台列表
 */
export function getSupportedPlatforms(): SocialPlatform[] {
  return Object.keys(adapterFactories) as SocialPlatform[];
}
