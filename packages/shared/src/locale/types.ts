/**
 * 多语言国际化类型定义
 */

/**
 * 语言消息条目
 * [源文本, 翻译文本]
 */
export type LocaleMessage = [string, string];

/**
 * 语言消息列表
 */
export type LocaleMessages = LocaleMessage[];

/**
 * 语言包配置
 */
export interface LocaleConfig {
  /**
   * 支持的语言列表
   */
  languages: string[];
  
  /**
   * 默认语言
   */
  defaultLocale: string;
}

/**
 * 语言信息
 */
export interface LanguageInfo {
  code: string;
  name: string;
  nativeName: string;
}

/**
 * 支持的语言列表
 */
export const SUPPORTED_LANGUAGES: LanguageInfo[] = [
  { code: 'zh-cn', name: 'Chinese (Simplified)', nativeName: '简体中文' },
  { code: 'zh-tw', name: 'Chinese (Traditional)', nativeName: '繁體中文' },
  { code: 'en', name: 'English', nativeName: 'English' },
  { code: 'ja', name: 'Japanese', nativeName: '日本語' },
  { code: 'ko', name: 'Korean', nativeName: '한국어' },
  { code: 'es', name: 'Spanish', nativeName: 'Español' },
  { code: 'fr', name: 'French', nativeName: 'Français' },
  { code: 'de', name: 'German', nativeName: 'Deutsch' },
  { code: 'it', name: 'Italian', nativeName: 'Italiano' },
  { code: 'pt', name: 'Portuguese', nativeName: 'Português' },
  { code: 'ru', name: 'Russian', nativeName: 'Русский' },
  { code: 'ar', name: 'Arabic', nativeName: 'العربية' },
  { code: 'th', name: 'Thai', nativeName: 'ไทย' },
  { code: 'vi', name: 'Vietnamese', nativeName: 'Tiếng Việt' },
];

/**
 * 获取语言信息
 */
export function getLanguageInfo(code: string): LanguageInfo | undefined {
  return SUPPORTED_LANGUAGES.find(lang => lang.code === code);
}
