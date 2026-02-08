/**
 * 多语言 React Context
 * 
 * 提供 React 组件中的多语言支持
 * 
 * @example
 * ```tsx
 * import { useLocale, useT } from '@/context/LocaleContext';
 * 
 * function MyComponent() {
 *   const { locale, setLocale } = useLocale();
 *   const t = useT();
 *   
 *   return (
 *     <div>
 *       <p>{t('设置')}</p>
 *       <button onClick={() => setLocale('en')}>English</button>
 *     </div>
 *   );
 * }
 * ```
 */

import React, { createContext, useContext, useState, useEffect, useCallback, useMemo } from 'react';
import {
  getLocale,
  setLocale as setLocaleService,
  onLocaleChange,
  t as tService,
  $t as $tService,
  SUPPORTED_LANGUAGES,
  type LanguageInfo,
} from '@sprouty-ai/shared/locale';
import { setupLocales, localeConfig } from '../locales';

/**
 * Locale Context 类型
 */
interface LocaleContextType {
  /**
   * 当前语言
   */
  locale: string;
  
  /**
   * 设置当前语言
   */
  setLocale: (locale: string) => void;
  
  /**
   * 支持的语言列表
   */
  languages: LanguageInfo[];
  
  /**
   * 基础翻译方法
   */
  t: (text: string) => string;
  
  /**
   * 带参数的翻译方法
   */
  $t: (text: string, params: Record<string, string | number>) => string;
}

// 创建 Context
const LocaleContext = createContext<LocaleContextType | null>(null);

/**
 * LocaleProvider 组件
 * 
 * 包装应用根组件，提供多语言支持
 */
export function LocaleProvider({ children }: { children: React.ReactNode }) {
  // 初始化语言包
  const [initialized, setInitialized] = useState(false);
  const [locale, setLocaleState] = useState(localeConfig.defaultLocale);
  
  // 初始化
  useEffect(() => {
    try {
      setupLocales(localeConfig.defaultLocale);
      setLocaleState(getLocale());
      setInitialized(true);
      
      // 监听语言变化
      const unsubscribe = onLocaleChange((newLocale) => {
        setLocaleState(newLocale);
      });
      
      return unsubscribe;
    } catch (error) {
      console.error('[LocaleProvider] Initialization failed:', error);
      // Fallback: Initialize anyway to prevent white screen
      setInitialized(true);
    }
  }, []);
  
  // 设置语言
  const setLocale = useCallback((newLocale: string) => {
    setLocaleService(newLocale);
  }, []);
  
  // 获取支持的语言列表（过滤出配置中的语言）
  const languages = useMemo(() => {
    return SUPPORTED_LANGUAGES.filter(lang => 
      localeConfig.languages.includes(lang.code)
    );
  }, []);
  
  // 创建翻译函数（带有状态依赖，触发重新渲染）
  const t = useCallback((text: string) => {
    return tService(text);
  }, [locale]); // locale 变化时重新创建函数
  
  const $t = useCallback((text: string, params: Record<string, string | number>) => {
    return $tService(text, params);
  }, [locale]);
  
  // Context 值
  const contextValue = useMemo<LocaleContextType>(() => ({
    locale,
    setLocale,
    languages,
    t,
    $t,
  }), [locale, setLocale, languages, t, $t]);
  
  // 等待初始化完成
  if (!initialized) {
    return null;
  }
  
  return (
    <LocaleContext.Provider value={contextValue}>
      {children}
    </LocaleContext.Provider>
  );
}

/**
 * useLocale Hook
 * 
 * 获取当前语言和切换语言的方法
 */
export function useLocale() {
  const context = useContext(LocaleContext);
  
  if (!context) {
    throw new Error('useLocale must be used within a LocaleProvider');
  }
  
  return {
    locale: context.locale,
    setLocale: context.setLocale,
    languages: context.languages,
  };
}

/**
 * useT Hook
 * 
 * 获取翻译函数
 * 
 * @example
 * ```tsx
 * const t = useT();
 * return <div>{t('设置')}</div>;
 * ```
 */
export function useT() {
  const context = useContext(LocaleContext);
  
  if (!context) {
    throw new Error('useT must be used within a LocaleProvider');
  }
  
  return context.t;
}

/**
 * useParamT Hook
 * 
 * 获取带参数的翻译函数
 * 
 * @example
 * ```tsx
 * const $t = useParamT();
 * return <div>{$t('欢迎{name}', { name: '用户' })}</div>;
 * ```
 */
export function useParamT() {
  const context = useContext(LocaleContext);
  
  if (!context) {
    throw new Error('useParamT must be used within a LocaleProvider');
  }
  
  return context.$t;
}

/**
 * @deprecated Use useParamT instead
 */
export const use$T = useParamT;

/**
 * useTranslation Hook
 * 
 * 获取 t 和 $t 函数
 * 
 * @example
 * ```tsx
 * const { t, $t } = useTranslation();
 * return (
 *   <div>
 *     <p>{t('设置')}</p>
 *     <p>{$t('欢迎{name}', { name: '用户' })}</p>
 *   </div>
 * );
 * ```
 */
export function useTranslation() {
  const context = useContext(LocaleContext);
  
  if (!context) {
    throw new Error('useTranslation must be used within a LocaleProvider');
  }
  
  return {
    t: context.t,
    $t: context.$t,
  };
}

// 导出 Context（供高级用例使用）
export { LocaleContext };
