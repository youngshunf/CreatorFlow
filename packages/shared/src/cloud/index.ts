/**
 * Cloud Integration Module
 * 
 * Provides cloud gateway configuration and subscription management.
 * 
 * Environment Configuration:
 * - Configuration is loaded from .env files (VITE_* variables)
 * - Supports: development, staging, production environments
 * - Set VITE_APP_ENV to switch environments
 */

// Cloud configuration
export {
  type CloudConfig,
  getCloudConfig,
  saveCloudConfig,
  clearCloudConfig,
  isCloudModeEnabled,
  getCloudGatewayUrl,
  enableCloudMode,
  disableCloudMode,
  getCloudAccessToken,
  getCloudLlmToken,
  updateCloudTokens,
  needsTokenRefresh,
} from './cloud-config.ts';

// Environment configuration
export {
  type AppEnvironment,
  type EnvConfig,
  getEnvConfig,
  loadEnvConfig,
  resetEnvConfig,
  getCurrentEnv,
  isDevelopment,
  isStaging,
  isProduction,
  isCloudModeEnabled as isEnvCloudModeEnabled,
  isDebugMode,
  getCloudApiUrl,
  getLlmGatewayUrl,
} from './env-config.ts';
