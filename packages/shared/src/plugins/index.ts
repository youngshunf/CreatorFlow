/**
 * Plugins Module
 *
 * 全局插件管理（存储在 ~/.sprouty-ai/plugins/）
 */

export {
  getGlobalPluginsDir,
  initializeGlobalPlugins,
  isGlobalPluginsInitialized,
  getGlobalPluginsMeta,
} from './global-plugins.ts';
