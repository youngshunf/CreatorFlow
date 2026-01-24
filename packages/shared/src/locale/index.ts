/**
 * 多语言国际化核心服务
 * 
 * 使用方式:
 * ```tsx
 * import { t, $t, initLocale, setLocale } from '@creator-flow/shared/locale';
 * 
 * // 初始化（默认使用中文）
 * initLocale('zh-cn');
 * 
 * // 基础翻译
 * t('设置')  // -> '设置' 或 'Settings'
 * 
 * // 带参数的翻译
 * $t('欢迎{name}', { name: '用户' })  // -> '欢迎用户' 或 'Welcome 用户'
 * 
 * // 切换语言
 * setLocale('en');
 * ```
 */

import type { LocaleMessages } from './types.ts';
import { getStoredLocale, setStoredLocale, getSystemLocale, normalizeLocale } from './storage.ts';

// 重新导出类型
export * from './types.ts';
export { getStoredLocale, setStoredLocale, getSystemLocale, normalizeLocale } from './storage.ts';

// 语言包映射对象
const messages: Record<string, LocaleMessages> = {};

// 当前语言
let currentLocale = 'zh-cn';

// 默认语言
let defaultLocale = 'zh-cn';

// 语言变化监听器
type LocaleChangeListener = (locale: string) => void;
const listeners: LocaleChangeListener[] = [];

/**
 * 获取当前语言
 */
export function getLocale(): string {
  return currentLocale;
}

/**
 * 设置当前语言
 */
export function setLocale(locale: string): void {
  const normalized = normalizeLocale(locale);
  
  if (currentLocale !== normalized) {
    currentLocale = normalized;
    setStoredLocale(normalized);
    
    // 通知所有监听器
    listeners.forEach(listener => listener(normalized));
  }
}

/**
 * 获取默认语言
 */
export function getDefaultLocale(): string {
  return defaultLocale;
}

/**
 * 添加语言变化监听器
 */
export function onLocaleChange(listener: LocaleChangeListener): () => void {
  listeners.push(listener);
  
  // 返回取消监听的函数
  return () => {
    const index = listeners.indexOf(listener);
    if (index > -1) {
      listeners.splice(index, 1);
    }
  };
}

/**
 * 追加语言包数据
 * 
 * @param name - 语言代码 (例如: 'en', 'zh-cn')
 * @param data - 语言包数据
 */
export function appendLocale(name: string, data: LocaleMessages): void {
  const normalized = normalizeLocale(name);
  
  if (!messages[normalized]) {
    messages[normalized] = [];
  }
  
  // 将新数据追加到前面，优先使用后添加的翻译
  messages[normalized].unshift(...data);
}

/**
 * 加载语言包
 * 
 * @param name - 语言代码
 * @param data - 语言包数据
 */
export function loadLocale(name: string, data: LocaleMessages): void {
  const normalized = normalizeLocale(name);
  messages[normalized] = data;
}

/**
 * 获取所有已加载的语言
 */
export function getLoadedLocales(): string[] {
  return Object.keys(messages);
}

/**
 * 检查语言包是否已加载
 */
export function isLocaleLoaded(locale: string): boolean {
  return normalizeLocale(locale) in messages;
}

/**
 * 基础翻译方法
 * 
 * @param text - 源文本（中文）
 * @returns 翻译后的文本
 * 
 * @example
 * t('设置')  // -> 'Settings' (当语言为英文时)
 */
export function t(text: string): string {
  const localeMessages = messages[currentLocale];
  
  if (!localeMessages) {
    return text;
  }
  
  // 查找翻译
  const entry = localeMessages.find(([source]) => source === text);
  
  if (entry && entry[1] && entry[1] !== '') {
    return entry[1];
  }
  
  return text;
}

/**
 * 带参数的翻译方法
 * 
 * @param text - 源文本（包含 {参数名} 占位符）
 * @param params - 参数对象
 * @returns 翻译并替换参数后的文本
 * 
 * @example
 * $t('欢迎{name}', { name: '张三' })  // -> 'Welcome 张三'
 * $t('共{count}条记录', { count: 10 })  // -> '共10条记录'
 */
export function $t(text: string, params: Record<string, string | number>): string {
  let result = t(text);
  
  // 替换参数
  if (params) {
    Object.entries(params).forEach(([key, value]) => {
      const strValue = typeof value === 'number' ? value.toString() : value;
      result = result.replaceAll(`{${key}}`, strValue);
    });
  }
  
  return result;
}

/**
 * 初始化多语言设置
 * 
 * @param locale - 默认语言代码，'system' 表示跟随系统
 * 
 * @example
 * // 默认使用中文
 * initLocale('zh-cn');
 * 
 * // 跟随系统语言
 * initLocale('system');
 */
export function initLocale(locale: string = 'zh-cn'): void {
  // 设置默认语言
  if (locale === 'system') {
    defaultLocale = getSystemLocale();
  } else {
    defaultLocale = normalizeLocale(locale);
  }
  
  // 优先使用存储的语言偏好
  const storedLocale = getStoredLocale();
  
  if (storedLocale) {
    currentLocale = normalizeLocale(storedLocale);
  } else {
    currentLocale = defaultLocale;
  }
}

/**
 * 重置为默认语言
 */
export function resetLocale(): void {
  setLocale(defaultLocale);
}
