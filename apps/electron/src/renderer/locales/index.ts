/**
 * 语言包加载器
 * 
 * 自动加载所有语言包并注册到 locale 服务
 */

import { loadLocale, initLocale } from '@creator-flow/shared/locale';
import type { LocaleMessages } from '@creator-flow/shared/locale';

// 导入语言包
import zhCN from './zh-cn.json';
import en from './en.json';

/**
 * 语言包配置
 */
export const localeConfig = {
  // 默认语言
  defaultLocale: 'zh-cn',
  
  // 支持的语言列表
  languages: ['zh-cn', 'en'],
};

/**
 * 初始化语言包
 * 在应用启动时调用此函数
 */
export function setupLocales(defaultLocale: string = 'zh-cn'): void {
  // 加载语言包
  loadLocale('zh-cn', zhCN as LocaleMessages);
  loadLocale('en', en as LocaleMessages);
  
  // 初始化 locale 服务
  initLocale(defaultLocale);
}

// 导出语言包数据（供其他模块使用）
export const locales = {
  'zh-cn': zhCN as LocaleMessages,
  'en': en as LocaleMessages,
};
