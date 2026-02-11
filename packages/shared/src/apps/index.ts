/**
 * Application Framework
 *
 * Public API for the application framework.
 * Applications provide pre-configured bundles for specific use cases.
 */

// Types
export type {
  AppManifest,
  AppDependencies,
  AppCapabilities,
  AppPricing,
  AppAuthor,
  WorkspaceInitConfig,
  DirectoryStructure,
  PresetData,
  MarketplaceInfo,
  PluginCompatibility,
  AppCompatibility,
  LoadedApp,
  AppSummary,
  InstallAppOptions,
  SkillReference,
  OSType,
} from './types.ts';

// Storage
export {
  getAppsDir,
  ensureAppsDir,
  getAppPath,
  loadApp,
  loadAppById,
  listAllApps,
  listMainApps,
  listPluginApps,
  getAppSummaries,
  saveAppManifest,
  createAppDirectory,
  deleteApp,
  appExists,
  getCompatiblePlugins,
} from './storage.ts';

// Initializer
export type {
  InitializeWorkspaceResult,
  InitializeWorkspaceOptions,
} from './initializer.ts';

export {
  initializeWorkspaceFromApp,
  installPluginToWorkspace,
  uninstallPluginFromWorkspace,
  migrateWorkspaceToApp,
  installSkillsFromCloud,
  createDirectoryStructure,
} from './initializer.ts';
