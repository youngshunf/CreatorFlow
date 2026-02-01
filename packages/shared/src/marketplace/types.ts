/**
 * Marketplace Types
 *
 * Type definitions for skill/app marketplace integration.
 * These types align with the cloud backend API schema.
 */

// ============================================================
// Market Item Types
// ============================================================

/**
 * Skill information from marketplace
 */
export interface MarketplaceSkill {
  skill_id: string;
  name: string;
  description: string | null;
  icon_url: string | null;
  author_id: number | null;
  author_name: string | null;
  category: string | null;
  tags: string | null;
  pricing_type: 'free' | 'paid';
  price: number;
  is_private: boolean;
  is_official: boolean;
  download_count: number;
  latest_version?: string | null;
  created_at?: string;
  updated_at?: string;
}

/**
 * Skill version information
 */
export interface MarketplaceSkillVersion {
  id: number;
  skill_id: string;
  version: string;
  changelog: string | null;
  package_url: string | null;
  file_hash: string | null;
  file_size: number | null;
  is_latest: boolean;
  published_at: string;
}

/**
 * App information from marketplace
 */
export interface MarketplaceApp {
  app_id: string;
  name: string;
  description: string | null;
  icon_url: string | null;
  author_id: number | null;
  author_name: string | null;
  pricing_type: 'free' | 'paid' | 'subscription';
  price: number;
  is_private: boolean;
  is_official: boolean;
  download_count: number;
  skill_dependencies: string | null;
  latest_version?: string | null;
  created_at?: string;
  updated_at?: string;
}

/**
 * App version information
 */
export interface MarketplaceAppVersion {
  id: number;
  app_id: string;
  version: string;
  changelog: string | null;
  skill_dependencies_versioned: Record<string, string> | null;
  package_url: string | null;
  file_hash: string | null;
  file_size: number | null;
  is_latest: boolean;
  published_at: string;
}

/**
 * Category information
 */
export interface MarketplaceCategory {
  id: number;
  slug: string;
  name: string;
  icon: string | null;
  parent_slug: string | null;
  sort_order: number;
}

// ============================================================
// API Response Types
// ============================================================

/**
 * Paginated list response
 */
export interface PaginatedResponse<T> {
  items: T[];
  total: number;
  page?: number;
  size?: number;
}

/**
 * Download response with presigned URL
 */
export interface DownloadResponse {
  download_url: string;
  file_hash: string | null;
  file_size: number | null;
}

/**
 * App download response with skill dependencies
 */
export interface AppDownloadResponse extends DownloadResponse {
  skill_dependencies: SkillDependencyDownload[] | null;
}

/**
 * Skill dependency download info
 */
export interface SkillDependencyDownload {
  id: string;
  version: string;
  download_url: string;
  file_hash: string | null;
}

/**
 * Sync request - list of installed items
 */
export interface SyncRequest {
  installed: InstalledItem[];
}

/**
 * Installed item info for sync
 */
export interface InstalledItem {
  id: string;
  version: string;
  type: 'skill' | 'app';
}

/**
 * Sync response - list of available updates
 */
export interface SyncResponse {
  updates: UpdateItem[];
}

/**
 * Update item info
 */
export interface UpdateItem {
  id: string;
  latest_version: string;
  changelog: string | null;
  type: 'skill' | 'app';
}

/**
 * Search response
 */
export interface SearchResponse {
  skills: MarketplaceSkill[];
  apps: MarketplaceApp[];
  total_skills: number;
  total_apps: number;
}

// ============================================================
// Local Storage Types
// ============================================================

/**
 * Local skill metadata for version tracking
 * Stored in {workspace}/skills/{skill_id}/.skill-meta.json
 */
export interface SkillMeta {
  /** Original skill ID from marketplace */
  sourceId: string;
  /** Base version from marketplace */
  baseVersion: string;
  /** Local version (e.g., "1.0.0-local.2") */
  localVersion: string;
  /** Whether the skill has been modified locally */
  isModified: boolean;
  /** Last sync timestamp */
  lastSynced: string;
  /** Whether this is a private/local-only version */
  isPrivate?: boolean;
}

/**
 * Cached marketplace metadata
 * Stored in ~/.creator-flow/marketplace/cache/meta.json
 */
export interface MarketplaceCacheData {
  /** Cache timestamp */
  cachedAt: string;
  /** Skills list */
  skills: MarketplaceSkill[];
  /** Apps list */
  apps: MarketplaceApp[];
  /** Categories list */
  categories: MarketplaceCategory[];
}

/**
 * Installed skill info (combination of local + meta)
 */
export interface InstalledSkillInfo {
  /** Skill ID */
  skillId: string;
  /** Display name */
  name: string;
  /** Description */
  description: string;
  /** Installed version */
  version: string;
  /** Whether has local modifications */
  isModified: boolean;
  /** Whether update is available */
  hasUpdate: boolean;
  /** Latest version if update available */
  latestVersion?: string;
  /** Path to skill directory */
  path: string;
}

// ============================================================
// Installation Types
// ============================================================

/**
 * Install progress callback
 */
export type InstallProgressCallback = (progress: InstallProgress) => void;

/**
 * Install progress info
 */
export interface InstallProgress {
  stage: 'downloading' | 'extracting' | 'installing' | 'complete' | 'error';
  percent: number;
  message: string;
  skillId?: string;
  error?: string;
}

/**
 * Install result
 */
export interface InstallResult {
  success: boolean;
  skillId: string;
  version: string;
  path?: string;
  error?: string;
}
