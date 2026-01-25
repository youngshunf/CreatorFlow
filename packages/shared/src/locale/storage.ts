/**
 * 语言偏好存储服务
 */

/* eslint-disable @typescript-eslint/no-explicit-any */
declare const window: any;
declare const navigator: any;

const LOCALE_STORAGE_KEY = 'creator-flow-locale';

/**
 * 获取存储的语言偏好
 */
export function getStoredLocale(): string | null {
  if (typeof window === 'undefined' || !window.localStorage) {
    return null;
  }
  
  try {
    return window.localStorage.getItem(LOCALE_STORAGE_KEY);
  } catch {
    return null;
  }
}

/**
 * 设置存储的语言偏好
 */
export function setStoredLocale(locale: string): void {
  if (typeof window === 'undefined' || !window.localStorage) {
    return;
  }
  
  try {
    window.localStorage.setItem(LOCALE_STORAGE_KEY, locale);
  } catch {
    // 存储失败时静默处理
  }
}

/**
 * 清除存储的语言偏好
 */
export function clearStoredLocale(): void {
  if (typeof window === 'undefined' || !window.localStorage) {
    return;
  }
  
  try {
    window.localStorage.removeItem(LOCALE_STORAGE_KEY);
  } catch {
    // 清除失败时静默处理
  }
}

/**
 * 获取系统语言
 */
export function getSystemLocale(): string {
  if (typeof navigator === 'undefined') {
    return 'en';
  }
  
  // 获取浏览器语言
  const browserLang = navigator.language || (navigator as any).userLanguage || 'en';
  
  // 标准化语言代码
  return normalizeLocale(browserLang);
}

/**
 * 标准化语言代码
 * 例如: zh-CN -> zh-cn, en-US -> en
 */
export function normalizeLocale(locale: string): string {
  const normalized = locale.toLowerCase().replace('_', '-');
  
  // 处理常见的语言代码映射
  const mapping: Record<string, string> = {
    'zh': 'zh-cn',
    'zh-hans': 'zh-cn',
    'zh-hant': 'zh-tw',
    'en-us': 'en',
    'en-gb': 'en',
  };
  
  // 尝试直接匹配
  if (mapping[normalized]) {
    return mapping[normalized];
  }
  
  // 尝试匹配语言前缀
  const prefix = normalized.split('-')[0];
  if (prefix && mapping[prefix]) {
    return mapping[prefix];
  }
  
  return normalized;
}
