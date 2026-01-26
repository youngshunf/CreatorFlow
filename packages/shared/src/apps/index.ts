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
  getBundledAppsDir,
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
  isAppBundled,
  getCompatiblePlugins,
  registerBundledApp,
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
} from './initializer.ts';

// Bundled Apps
export {
  GENERAL_APP,
  CREATOR_MEDIA_APP,
  registerBundledApps,
  getBundledAppManifests,
} from './bundled-apps.ts';
