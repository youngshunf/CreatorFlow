/**
 * Environment Configuration Manager
 * 
 * Manages environment-specific configuration for cloud services.
 * Automatically loads configuration based on VITE_APP_ENV.
 * 
 * Supported environments:
 * - development: Local development (localhost)
 * - staging: Test/staging environment
 * - production: Production environment
 */

// Environment types
export type AppEnvironment = 'development' | 'staging' | 'production';

// Environment configuration interface
export interface EnvConfig {
  /** Current environment */
  env: AppEnvironment;
  /** Cloud API base URL */
  cloudApiUrl: string;
  /** LLM gateway path (relative to API URL) */
  llmGatewayPath: string;
  /** Full LLM gateway URL */
  llmGatewayUrl: string;
  /** Whether cloud mode is enabled */
  cloudModeEnabled: boolean;
  /** Whether debug mode is enabled */
  debugMode: boolean;
}

// Default configuration values
const DEFAULT_CONFIG: Record<AppEnvironment, Omit<EnvConfig, 'env'>> = {
  development: {
    cloudApiUrl: 'http://localhost:8020/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'http://localhost:8020/api/v1/llm/proxy',
    cloudModeEnabled: true,
    debugMode: true,
  },
  staging: {
    // Default to localhost for local development/staging
    // For remote staging server, set VITE_CLOUD_API_URL in .env.staging
    cloudApiUrl: 'http://localhost:8020/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'http://localhost:8020/api/v1/llm/proxy',
    cloudModeEnabled: true,
    debugMode: true,
  },
  production: {
    cloudApiUrl: 'https://api.creator-flow.com/api/v1',
    llmGatewayPath: '/llm/proxy',
    llmGatewayUrl: 'https://api.creator-flow.com/api/v1/llm/proxy',
    cloudModeEnabled: true,
    debugMode: false,
  },
};

/**
 * Get Vite environment variables safely
 * 
 * This function is designed to work around esbuild warnings about import.meta
 * in CJS output. We use Function constructor to dynamically access import.meta
 * which prevents esbuild from statically analyzing and warning about it.
 */
function getViteEnv(): Record<string, string> | undefined {
  // Only attempt in browser/renderer context
  if (typeof window === 'undefined') {
    return undefined;
  }
  
  try {
    // Use indirect eval to access import.meta without triggering esbuild warning
    // This is safe because it only runs in the browser context where ESM is used
    // eslint-disable-next-line no-new-func
    const getImportMeta = new Function('return typeof import.meta !== "undefined" ? import.meta : undefined');
    const meta = getImportMeta();
    return meta?.env;
  } catch {
    return undefined;
  }
}

/**
 * Get environment variable value (works in both Vite and Node.js)
 * 
 * Priority:
 * 1. process.env (Node.js / Electron main process)
 * 2. import.meta.env (Vite / browser renderer)
 */
function getEnvVar(key: string): string | undefined {
  // Node.js / Electron main process environment
  if (typeof process !== 'undefined' && process.env) {
    return process.env[key];
  }
  
  // Vite environment (browser/renderer)
  const viteEnv = getViteEnv();
  if (viteEnv) {
    return viteEnv[key];
  }
  
  return undefined;
}

/**
 * Parse boolean environment variable
 */
function parseBool(value: string | undefined, defaultValue: boolean): boolean {
  if (value === undefined) return defaultValue;
  return value.toLowerCase() === 'true' || value === '1';
}

/**
 * Get current environment from VITE_APP_ENV
 */
export function getCurrentEnv(): AppEnvironment {
  const env = getEnvVar('VITE_APP_ENV');
  if (env === 'staging' || env === 'production') {
    return env;
  }
  return 'development';
}

/**
 * Load environment configuration
 * 
 * Priority:
 * 1. Environment variables (VITE_*)
 * 2. Default values for current environment
 */
export function loadEnvConfig(): EnvConfig {
  const env = getCurrentEnv();
  const defaults = DEFAULT_CONFIG[env];
  
  // Get values from environment variables or use defaults
  const cloudApiUrl = getEnvVar('VITE_CLOUD_API_URL') || defaults.cloudApiUrl;
  const llmGatewayPath = getEnvVar('VITE_LLM_GATEWAY_PATH') || defaults.llmGatewayPath;
  const cloudModeEnabled = parseBool(getEnvVar('VITE_CLOUD_MODE_ENABLED'), defaults.cloudModeEnabled);
  const debugMode = parseBool(getEnvVar('VITE_DEBUG_MODE'), defaults.debugMode);
  
  // Construct full gateway URL
  const baseUrl = cloudApiUrl.replace(/\/$/, ''); // Remove trailing slash
  const gatewayPath = llmGatewayPath.startsWith('/') ? llmGatewayPath : `/${llmGatewayPath}`;
  const llmGatewayUrl = `${baseUrl}${gatewayPath}`;
  
  return {
    env,
    cloudApiUrl,
    llmGatewayPath,
    llmGatewayUrl,
    cloudModeEnabled,
    debugMode,
  };
}

// Cached configuration (loaded once)
let cachedConfig: EnvConfig | null = null;

/**
 * Get environment configuration (cached)
 */
export function getEnvConfig(): EnvConfig {
  if (!cachedConfig) {
    cachedConfig = loadEnvConfig();
  }
  return cachedConfig;
}

/**
 * Reset cached configuration (for testing)
 */
export function resetEnvConfig(): void {
  cachedConfig = null;
}

/**
 * Check if current environment is development
 */
export function isDevelopment(): boolean {
  return getEnvConfig().env === 'development';
}

/**
 * Check if current environment is staging
 */
export function isStaging(): boolean {
  return getEnvConfig().env === 'staging';
}

/**
 * Check if current environment is production
 */
export function isProduction(): boolean {
  return getEnvConfig().env === 'production';
}

/**
 * Check if cloud mode is enabled
 */
export function isCloudModeEnabled(): boolean {
  return getEnvConfig().cloudModeEnabled;
}

/**
 * Check if debug mode is enabled
 */
export function isDebugMode(): boolean {
  return getEnvConfig().debugMode;
}

/**
 * Get cloud API URL for current environment
 */
export function getCloudApiUrl(): string {
  return getEnvConfig().cloudApiUrl;
}

/**
 * Get LLM gateway URL for current environment
 */
export function getLlmGatewayUrl(): string {
  return getEnvConfig().llmGatewayUrl;
}
