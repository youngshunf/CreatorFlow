/**
 * Cloud LLM Configuration (Browser-safe)
 *
 * Manages cloud gateway configuration for LLM calls.
 * This module is designed to work in browser/renderer context.
 *
 * Authentication:
 * - The cloud gateway uses x-api-key header for authentication
 * - LLM Token (sk-cf-xxx) serves as both authentication and rate limiting key
 * - No JWT required - simplifies desktop integration
 *
 * Environment Configuration:
 * - Configuration is defined in packages/shared/src/config/environments.ts
 * - Environment is selected at build time via APP_ENV
 * - Supports development, staging, and production environments
 *
 * Note: This module does NOT import Node.js-specific modules (crypto, fs, etc.)
 * to ensure browser compatibility. Credential storage is handled via IPC in Electron.
 */

import {
  getEnvironmentConfig,
  getLlmGatewayUrl,
  getCloudApiUrl,
  getCurrentEnv,
} from '../config/environments';

// Cloud configuration interface
export interface CloudConfig {
  /** Cloud API base URL (from env or override) */
  apiBaseUrl: string;
  /** JWT access token for cloud API calls (not needed for LLM gateway) */
  accessToken: string;
  /** LLM token for API key authentication (used as x-api-key) */
  llmToken: string;
  /** Token expiry time */
  expiresAt?: number;
}

// Local storage keys
const CLOUD_CONFIG_KEY = 'cloud_config';

/**
 * Get cloud configuration from storage
 */
export function getCloudConfig(): CloudConfig | null {
  try {
    const stored = localStorage.getItem(CLOUD_CONFIG_KEY);
    if (!stored) return null;
    return JSON.parse(stored) as CloudConfig;
  } catch {
    return null;
  }
}

/**
 * Save cloud configuration to storage
 */
export function saveCloudConfig(config: CloudConfig): void {
  localStorage.setItem(CLOUD_CONFIG_KEY, JSON.stringify(config));
}

/**
 * Clear cloud configuration from storage
 */
export function clearCloudConfig(): void {
  localStorage.removeItem(CLOUD_CONFIG_KEY);
}

/**
 * Check if cloud mode is enabled
 */
export function isCloudModeEnabled(): boolean {
  const config = getCloudConfig();
  return config !== null;
}

/**
 * Get the cloud LLM gateway URL
 * Returns the full URL to the cloud LLM proxy endpoint
 * Uses environment configuration as default
 */
export function getCloudGatewayUrl(): string {
  // First try to get from saved cloud config
  const config = getCloudConfig();
  if (config?.apiBaseUrl) {
    const envConfig = getEnvironmentConfig();
    return `${config.apiBaseUrl.replace(/\/$/, '')}${envConfig.llmGatewayPath}`;
  }

  // Fall back to environment configuration
  return getLlmGatewayUrl();
}

/**
 * Configure SDK to use cloud gateway (browser-safe)
 * 
 * This saves cloud config to localStorage. The main process
 * will read this config and configure the SDK accordingly.
 * 
 * @param config Cloud configuration from login
 */
export function enableCloudMode(config: CloudConfig): void {
  const envConfig = getEnvironmentConfig();
  console.log(`[CloudConfig] Enabling cloud mode (env: ${envConfig.name})`);

  // Use API URL from env config if not provided
  if (!config.apiBaseUrl) {
    config.apiBaseUrl = envConfig.cloudApiUrl;
  }

  // Save cloud config to localStorage
  // Main process will read this and configure SDK
  saveCloudConfig(config);

  // Also save gateway URL for easy access
  const gatewayUrl = `${config.apiBaseUrl.replace(/\/$/, '')}${envConfig.llmGatewayPath}`;
  localStorage.setItem('cloud_gateway_url', gatewayUrl);
  localStorage.setItem('cloud_llm_token', config.llmToken);

  console.log(`[CloudConfig] Cloud mode enabled, gateway: ${gatewayUrl}`);
}

/**
 * Disable cloud mode and clear configuration
 */
export function disableCloudMode(): void {
  console.log('[CloudConfig] Disabling cloud mode');
  
  // Clear cloud config from localStorage
  clearCloudConfig();
  localStorage.removeItem('cloud_gateway_url');
  localStorage.removeItem('cloud_llm_token');
  
  console.log('[CloudConfig] Cloud mode disabled');
}

/**
 * Get access token for cloud API calls
 */
export function getCloudAccessToken(): string | null {
  const config = getCloudConfig();
  return config?.accessToken || null;
}

/**
 * Get LLM token for cloud gateway authentication
 */
export function getCloudLlmToken(): string | null {
  const config = getCloudConfig();
  return config?.llmToken || null;
}

/**
 * Update tokens after refresh
 */
export function updateCloudTokens(accessToken: string, llmToken?: string): void {
  const config = getCloudConfig();
  if (!config) return;
  
  config.accessToken = accessToken;
  if (llmToken) {
    config.llmToken = llmToken;
  }
  saveCloudConfig(config);
}

/**
 * Check if cloud tokens need refresh
 */
export function needsTokenRefresh(): boolean {
  const config = getCloudConfig();
  if (!config || !config.expiresAt) return false;
  
  // Refresh if within 5 minutes of expiry
  const bufferMs = 5 * 60 * 1000;
  return Date.now() > config.expiresAt - bufferMs;
}
