/**
 * Cloud Integration Module
 *
 * Provides cloud gateway configuration and subscription management.
 *
 * Environment Configuration:
 * - Configuration is defined in packages/shared/src/config/environments.ts
 * - Environment is selected at build time via APP_ENV
 * - Supports: development, staging, production environments
 *
 * Usage:
 *   - Development: bun run electron:dev (default)
 *   - Staging: APP_ENV=staging bun run electron:build
 *   - Production: APP_ENV=production bun run electron:build
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
