/**
 * Environment Configuration
 *
 * Static environment configurations for different deployment targets.
 * The active environment is determined at build time via __APP_ENV__.
 *
 * Usage:
 *   - Development: bun run electron:dev (default)
 *   - Staging: APP_ENV=staging bun run electron:build
 *   - Production: APP_ENV=production bun run electron:build
 */

// Environment types
export type AppEnvironment = 'development' | 'staging' | 'production';

// Environment configuration interface
export interface EnvironmentConfig {
  /** Environment name */
  name: AppEnvironment;
  /** Cloud API base URL */
  cloudApiUrl: string;
  /** LLM gateway path (relative to API URL) */
  llmGatewayPath: string;
  /** Full LLM gateway URL */
  llmGatewayUrl: string;
  /** Whether debug mode is enabled */
  debugMode: boolean;
  /** Whether to show DevTools by default */
  showDevTools: boolean;
}

// Environment configurations
const environments: Record<AppEnvironment, EnvironmentConfig> = {
  development: {
    name: 'development',
    cloudApiUrl: 'http://localhost:8020/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'http://localhost:8020/api/v1/llm/proxy',
    debugMode: true,
    showDevTools: true,
  },
  staging: {
    name: 'staging',
    cloudApiUrl: 'http://192.168.1.92:8020/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'http://192.168.1.92:8020/api/v1/llm/proxy',
    debugMode: true,
    showDevTools: true,
  },
  production: {
    name: 'production',
    cloudApiUrl: 'http://api.ai.dcfuture.cn/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'http://api.ai.dcfuture.cn/api/v1/llm/proxy',
    debugMode: false,
    showDevTools: false,
  },
};

// Build-time injected environment (set by esbuild/vite --define)
declare const __APP_ENV__: AppEnvironment | undefined;

/**
 * Get the current environment name
 */
export function getCurrentEnv(): AppEnvironment {
  // Check build-time injected value
  if (typeof __APP_ENV__ !== 'undefined' && __APP_ENV__ in environments) {
    return __APP_ENV__;
  }
  // Default to development
  return 'development';
}

/**
 * Get the configuration for the current environment
 */
export function getEnvironmentConfig(): EnvironmentConfig {
  return environments[getCurrentEnv()];
}

/**
 * Get the configuration for a specific environment
 */
export function getEnvironmentConfigFor(env: AppEnvironment): EnvironmentConfig {
  return environments[env];
}

// Convenience exports
export const currentEnv = getCurrentEnv();
export const config = getEnvironmentConfig();

// Helper functions for common checks
export function isDevelopment(): boolean {
  return getCurrentEnv() === 'development';
}

export function isStaging(): boolean {
  return getCurrentEnv() === 'staging';
}

export function isProduction(): boolean {
  return getCurrentEnv() === 'production';
}

export function isDebugMode(): boolean {
  return getEnvironmentConfig().debugMode;
}

// Cloud API helpers
export function getCloudApiUrl(): string {
  return getEnvironmentConfig().cloudApiUrl;
}

export function getLlmGatewayUrl(): string {
  return getEnvironmentConfig().llmGatewayUrl;
}
