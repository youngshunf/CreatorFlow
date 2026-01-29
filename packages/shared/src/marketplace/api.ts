/**
 * Marketplace API Client
 *
 * HTTP client for communicating with the cloud marketplace API.
 * Handles skill/app listing, downloading, and sync operations.
 */

import type {
  MarketplaceSkill,
  MarketplaceSkillVersion,
  MarketplaceApp,
  MarketplaceAppVersion,
  MarketplaceCategory,
  PaginatedResponse,
  DownloadResponse,
  AppDownloadResponse,
  SyncRequest,
  SyncResponse,
  SearchResponse,
} from './types.ts';
import { getCloudConfig } from '../cloud/cloud-config.ts';
import { getEnvConfig } from '../cloud/env-config.ts';

// ============================================================
// Configuration
// ============================================================

/**
 * Get the marketplace API base URL
 * Uses /client prefix for public desktop API endpoints
 * 
 * Note: cloudApiUrl already contains /api/v1, so we only append /marketplace/client
 */
function getMarketplaceApiUrl(): string {
  const config = getCloudConfig();
  if (config?.apiBaseUrl) {
    // apiBaseUrl might or might not contain /api/v1, normalize it
    const base = config.apiBaseUrl.replace(/\/$/, '').replace(/\/api\/v1$/, '');
    return `${base}/api/v1/marketplace/client`;
  }
  const envConfig = getEnvConfig();
  // cloudApiUrl already contains /api/v1, so just append /marketplace/client
  return `${envConfig.cloudApiUrl.replace(/\/$/, '')}/marketplace/client`;
}

/**
 * Get auth headers for API calls
 */
function getAuthHeaders(): Record<string, string> {
  const config = getCloudConfig();
  const headers: Record<string, string> = {
    'Content-Type': 'application/json',
  };
  if (config?.accessToken) {
    headers['Authorization'] = `Bearer ${config.accessToken}`;
  }
  return headers;
}

/**
 * Wrapper for fetch with error handling
 */
async function apiFetch<T>(
  path: string,
  options: RequestInit = {}
): Promise<T> {
  const url = `${getMarketplaceApiUrl()}${path}`;
  const response = await fetch(url, {
    ...options,
    headers: {
      ...getAuthHeaders(),
      ...(options.headers || {}),
    },
  });

  if (!response.ok) {
    const errorText = await response.text();
    throw new Error(`Marketplace API error: ${response.status} - ${errorText}`);
  }

  const json = await response.json();
  // Handle wrapped response format: { code: 200, data: {...} }
  if (json && typeof json === 'object' && 'data' in json) {
    return json.data as T;
  }
  return json as T;
}

// ============================================================
// Skills API
// ============================================================

export interface ListSkillsParams {
  page?: number;
  size?: number;
  category?: string;
  tags?: string;
  pricing_type?: string;
  is_official?: boolean;
}

/**
 * List skills from marketplace
 */
export async function listSkills(
  params: ListSkillsParams = {}
): Promise<PaginatedResponse<MarketplaceSkill>> {
  const searchParams = new URLSearchParams();
  if (params.page) searchParams.set('page', params.page.toString());
  if (params.size) searchParams.set('size', params.size.toString());
  if (params.category) searchParams.set('category', params.category);
  if (params.tags) searchParams.set('tags', params.tags);
  if (params.pricing_type) searchParams.set('pricing_type', params.pricing_type);
  if (params.is_official !== undefined) {
    searchParams.set('is_official', params.is_official.toString());
  }

  const query = searchParams.toString();
  return apiFetch<PaginatedResponse<MarketplaceSkill>>(
    `/skills${query ? `?${query}` : ''}`
  );
}

/**
 * Get skill details by ID
 */
export async function getSkill(skillId: string): Promise<MarketplaceSkill> {
  return apiFetch<MarketplaceSkill>(`/skills/${skillId}`);
}

/**
 * Get skill versions
 */
export async function getSkillVersions(
  skillId: string
): Promise<PaginatedResponse<MarketplaceSkillVersion>> {
  return apiFetch<PaginatedResponse<MarketplaceSkillVersion>>(
    `/skills/versions?skill_id=${skillId}`
  );
}

/**
 * Get download URL for a skill
 */
export async function getSkillDownload(
  skillId: string,
  version: string = 'latest'
): Promise<DownloadResponse> {
  return apiFetch<DownloadResponse>(`/download/skill/${skillId}/${version}`);
}

// ============================================================
// Apps API
// ============================================================

export interface ListAppsParams {
  page?: number;
  size?: number;
  pricing_type?: string;
  is_official?: boolean;
}

/**
 * List apps from marketplace
 */
export async function listApps(
  params: ListAppsParams = {}
): Promise<PaginatedResponse<MarketplaceApp>> {
  const searchParams = new URLSearchParams();
  if (params.page) searchParams.set('page', params.page.toString());
  if (params.size) searchParams.set('size', params.size.toString());
  if (params.pricing_type) searchParams.set('pricing_type', params.pricing_type);
  if (params.is_official !== undefined) {
    searchParams.set('is_official', params.is_official.toString());
  }

  const query = searchParams.toString();
  return apiFetch<PaginatedResponse<MarketplaceApp>>(
    `/apps${query ? `?${query}` : ''}`
  );
}

/**
 * Get app details by ID
 */
export async function getApp(appId: string): Promise<MarketplaceApp> {
  return apiFetch<MarketplaceApp>(`/apps/${appId}`);
}

/**
 * Get app versions
 */
export async function getAppVersions(
  appId: string
): Promise<PaginatedResponse<MarketplaceAppVersion>> {
  return apiFetch<PaginatedResponse<MarketplaceAppVersion>>(
    `/apps/versions?app_id=${appId}`
  );
}

/**
 * Get download URL for an app (includes skill dependencies)
 */
export async function getAppDownload(
  appId: string,
  version: string = 'latest'
): Promise<AppDownloadResponse> {
  return apiFetch<AppDownloadResponse>(`/download/app/${appId}/${version}`);
}

// ============================================================
// Categories API
// ============================================================

/**
 * List all categories
 */
export async function listCategories(): Promise<PaginatedResponse<MarketplaceCategory>> {
  return apiFetch<PaginatedResponse<MarketplaceCategory>>('/categories');
}

// ============================================================
// Search API
// ============================================================

export interface SearchParams {
  q: string;
  type?: 'skill' | 'app' | 'all';
  category?: string;
  limit?: number;
}

/**
 * Search skills and apps
 */
export async function search(params: SearchParams): Promise<SearchResponse> {
  const searchParams = new URLSearchParams();
  searchParams.set('q', params.q);
  if (params.type) searchParams.set('type', params.type);
  if (params.category) searchParams.set('category', params.category);
  if (params.limit) searchParams.set('limit', params.limit.toString());

  return apiFetch<SearchResponse>(`/search?${searchParams.toString()}`);
}

// ============================================================
// Sync API
// ============================================================

/**
 * Sync installed items and check for updates
 */
export async function syncInstalled(
  request: SyncRequest
): Promise<SyncResponse> {
  return apiFetch<SyncResponse>('/sync', {
    method: 'POST',
    body: JSON.stringify(request),
  });
}
