/**
 * Marketplace Module
 *
 * Cloud marketplace integration for skills and apps.
 * Provides:
 * - API client for marketplace browsing and download
 * - Local storage and caching
 * - Skill installation and update management
 */

// Types
export type {
  // Market items
  MarketplaceSkill,
  MarketplaceSkillVersion,
  MarketplaceApp,
  MarketplaceAppVersion,
  MarketplaceCategory,
  // API responses
  PaginatedResponse,
  DownloadResponse,
  AppDownloadResponse,
  SkillDependencyDownload,
  SyncRequest,
  SyncResponse,
  InstalledItem,
  UpdateItem,
  SearchResponse,
  // Local storage
  SkillMeta,
  MarketplaceCacheData,
  InstalledSkillInfo,
  // Installation
  InstallProgress,
  InstallProgressCallback,
  InstallResult,
} from './types.ts';

// API Client
export {
  listSkills,
  getSkill,
  getSkillVersions,
  getSkillDownload,
  listApps,
  getApp,
  getAppVersions,
  getAppDownload,
  listCategories,
  search,
  syncInstalled,
  type ListSkillsParams,
  type ListAppsParams,
  type SearchParams,
} from './api.ts';

// Local Storage
export {
  // Path utilities
  getConfigRoot,
  getMarketplaceDir,
  getMarketplaceCacheDir,
  getMarketplaceSkillsDir,
  getMarketplaceAppsDir,
  // Cache operations
  getMarketplaceCache,
  saveMarketplaceCache,
  clearMarketplaceCache,
  // Skill meta operations
  getSkillMetaPath,
  readSkillMeta,
  writeSkillMeta,
  createSkillMeta,
  markSkillModified,
  // Installed skills
  getInstalledSkills,
  isSkillInstalled,
  getInstalledSkillVersion,
  // Package cache
  getCachedSkillPackagePath,
  isSkillPackageCached,
  getCachedAppPackagePath,
  isAppPackageCached,
} from './storage.ts';

// Skill Installer
export {
  downloadSkillPackage,
  installSkill,
  installSkills,
  updateSkill,
} from './skill-installer.ts';

// Update Checker
export {
  checkForUpdates,
  getSkillsWithUpdateStatus,
  hasSkillUpdate,
  startBackgroundSync,
  stopBackgroundSync,
  forceSync,
  compareVersions,
  isVersionOlder,
  isLocalVersion,
  type UpdateCallback,
} from './update-checker.ts';
