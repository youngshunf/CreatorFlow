/**
 * Environment Configuration Manager
 *
 * This module re-exports from the centralized environment configuration.
 * It provides backward compatibility for existing code that imports from here.
 *
 * The actual configuration is defined in @creator-flow/shared/config/environments.ts
 * and is determined at build time via the __APP_ENV__ variable.
 *
 * Usage:
 *   - Development: bun run electron:dev (default)
 *   - Staging: APP_ENV=staging bun run electron:build
 *   - Production: APP_ENV=production bun run electron:build
 */

import {
  type AppEnvironment,
  type EnvironmentConfig,
  getCurrentEnv,
  getEnvironmentConfig,
  isDevelopment,
  isStaging,
  isProduction,
  isDebugMode,
  getCloudApiUrl,
  getLlmGatewayUrl,
} from '../config/environments';

// Re-export types
export type { AppEnvironment };

// Legacy EnvConfig interface for backward compatibility
export interface EnvConfig {
  env: AppEnvironment;
  cloudApiUrl: string;
  llmGatewayPath: string;
  llmGatewayUrl: string;
  cloudModeEnabled: boolean;
  debugMode: boolean;
}

/**
 * Convert EnvironmentConfig to legacy EnvConfig format
 */
function toEnvConfig(config: EnvironmentConfig): EnvConfig {
  return {
    env: config.name,
    cloudApiUrl: config.cloudApiUrl,
    llmGatewayPath: config.llmGatewayPath,
    llmGatewayUrl: config.llmGatewayUrl,
    cloudModeEnabled: true, // Cloud mode is always enabled
    debugMode: config.debugMode,
  };
}

// Cached configuration
let cachedConfig: EnvConfig | null = null;

/**
 * Get environment configuration (cached)
 */
export function getEnvConfig(): EnvConfig {
  if (!cachedConfig) {
    cachedConfig = toEnvConfig(getEnvironmentConfig());
  }
  return cachedConfig;
}

/**
 * Load environment configuration (returns fresh config without caching)
 */
export function loadEnvConfig(): EnvConfig {
  return toEnvConfig(getEnvironmentConfig());
}

/**
 * Reset cached configuration (for testing)
 */
export function resetEnvConfig(): void {
  cachedConfig = null;
}

/**
 * Check if cloud mode is enabled (always true in this implementation)
 */
export function isCloudModeEnabled(): boolean {
  return true;
}

// Re-export functions from environments.ts
export {
  getCurrentEnv,
  isDevelopment,
  isStaging,
  isProduction,
  isDebugMode,
  getCloudApiUrl,
  getLlmGatewayUrl,
};
