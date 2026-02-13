/**
 * SourceServerBuilder
 *
 * Builds MCP and API server configurations from LoadedSource objects.
 * This module handles URL normalization and server config creation,
 * but does NOT fetch credentials - credentials are passed in.
 *
 * This replaces SourceService's server building logic with a cleaner
 * separation of concerns:
 * - SourceCredentialManager: handles credentials
 * - SourceServerBuilder: handles server configuration
 */

import { execSync } from 'node:child_process';
import { existsSync } from 'node:fs';
import { homedir } from 'node:os';
import { join } from 'node:path';
import type { LoadedSource, ApiConfig } from './types.ts';
import type { ApiCredential } from './credential-manager.ts';
import { isSourceUsable } from './storage.ts';
import { createApiServer, type SummarizeCallback } from './api-tools.ts';
import { createSdkMcpServer } from '@anthropic-ai/claude-agent-sdk';
import { debug } from '../utils/debug.ts';

/**
 * 命令解析器类型 — 将命令名（如 'bun'）解析为可执行文件的绝对路径。
 * 由宿主环境（如 Electron 主进程）注入，以覆盖默认的 which 查找逻辑。
 */
export type CommandResolver = (command: string) => string;

/**
 * 工作目录解析器类型 — 将 cwd 配置（可能含 'app:' 前缀）解析为实际路径。
 * 由宿主环境注入，以正确处理 Electron 特有的路径前缀。
 */
export type CwdResolver = (cwd: string | undefined, workspaceRootPath: string) => string | undefined;

/**
 * Well-known install locations for common runtimes.
 * Used as fallback when `which` fails (e.g. Electron subprocess with limited PATH).
 */
const WELL_KNOWN_PATHS: Record<string, string[]> = {
  bun: [
    join(homedir(), '.bun', 'bin', 'bun'),
    '/usr/local/bin/bun',
    '/opt/homebrew/bin/bun',
  ],
  node: [
    '/usr/local/bin/node',
    '/opt/homebrew/bin/node',
    join(homedir(), '.nvm', 'versions', 'node', 'current', 'bin', 'node'),
  ],
  npx: [
    '/usr/local/bin/npx',
    '/opt/homebrew/bin/npx',
  ],
  deno: [
    join(homedir(), '.deno', 'bin', 'deno'),
    '/usr/local/bin/deno',
    '/opt/homebrew/bin/deno',
  ],
};

/**
 * Resolve a command name to its absolute path.
 * Returns the original command if already absolute or resolution fails.
 *
 * Strategy:
 * 1. If already an absolute path, return as-is
 * 2. Try `which` to find it on PATH
 * 3. Check well-known install locations
 */
export function resolveCommand(command: string): string {
  // Already absolute
  if (command.startsWith('/')) return command;

  // Try `which`
  try {
    const resolved = execSync(`which ${command}`, {
      encoding: 'utf-8',
      timeout: 3000,
      env: {
        ...process.env,
        // Extend PATH with common locations Electron might miss
        PATH: [
          process.env.PATH,
          join(homedir(), '.bun', 'bin'),
          join(homedir(), '.deno', 'bin'),
          '/usr/local/bin',
          '/opt/homebrew/bin',
        ].filter(Boolean).join(':'),
      },
    }).trim();
    if (resolved) {
      debug(`[resolveCommand] Resolved '${command}' → '${resolved}' via which`);
      return resolved;
    }
  } catch {
    // which failed, try well-known paths
  }

  // Check well-known paths
  const candidates = WELL_KNOWN_PATHS[command];
  if (candidates) {
    for (const candidate of candidates) {
      if (existsSync(candidate)) {
        debug(`[resolveCommand] Resolved '${command}' → '${candidate}' via well-known path`);
        return candidate;
      }
    }
  }

  debug(`[resolveCommand] Could not resolve '${command}', using as-is`);
  return command;
}

/**
 * Standard error messages for server build failures.
 * Use these constants instead of string literals to ensure consistent matching.
 */
export const SERVER_BUILD_ERRORS = {
  AUTH_REQUIRED: 'Authentication required',
  CREDENTIALS_NEEDED: 'Credentials needed',
} as const;

/**
 * MCP server configuration compatible with Claude Agent SDK
 * Supports HTTP/SSE (remote) and stdio (local subprocess) transports.
 */
export type McpServerConfig =
  | { type: 'http' | 'sse'; url: string; headers?: Record<string, string> }
  | { type: 'stdio'; command: string; args?: string[]; env?: Record<string, string> };

/**
 * Source with its credential pre-loaded
 */
export interface SourceWithCredential {
  source: LoadedSource;
  /** Token for MCP sources, or ApiCredential for API sources */
  token?: string | null;
  credential?: ApiCredential | null;
}

/**
 * Result of building servers from sources
 */
export interface BuiltServers {
  /** MCP server configs keyed by source slug */
  mcpServers: Record<string, McpServerConfig>;
  /** In-process API servers keyed by source slug */
  apiServers: Record<string, ReturnType<typeof createSdkMcpServer>>;
  /** Sources that failed to build (missing auth, etc.) */
  errors: Array<{ sourceSlug: string; error: string }>;
}

/**
 * Resolve cwd path for stdio MCP servers.
 * Supports special prefixes:
 * - 'app:' - Relative to app installation directory
 * - Absolute path - Used as-is
 * - Relative path - Relative to workspace root
 */
export function resolveCwd(cwd: string | undefined, workspaceRootPath: string): string | undefined {
  if (!cwd) return undefined;

  // app: prefix - resolve relative to app installation directory
  if (cwd.startsWith('app:')) {
    const relativePath = cwd.slice(4); // Remove 'app:' prefix
    // In Electron, __dirname points to app.asar or app/ in development
    // We need to get the actual app root
    const appRoot = process.resourcesPath || join(__dirname, '..', '..');
    const resolved = join(appRoot, relativePath);
    debug(`[resolveCwd] Resolved 'app:${relativePath}' → '${resolved}'`);
    return resolved;
  }

  // Absolute path - use as-is
  if (cwd.startsWith('/')) {
    return cwd;
  }

  // Relative path - resolve relative to workspace root
  const resolved = join(workspaceRootPath, cwd);
  debug(`[resolveCwd] Resolved relative path '${cwd}' → '${resolved}'`);
  return resolved;
}

/**
 * SourceServerBuilder - builds server configs from sources
 *
 * Usage:
 * ```typescript
 * const builder = new SourceServerBuilder();
 *
 * // Build MCP server config
 * const mcpConfig = builder.buildMcpServer(source, token);
 *
 * // Build all servers from sources with credentials
 * const { mcpServers, apiServers, errors } = await builder.buildAll([
 *   { source, token: 'abc123' },
 *   { source: apiSource, credential: 'api-key' },
 * ]);
 * ```
 */
export class SourceServerBuilder {
  private commandResolver: CommandResolver;
  private cwdResolver: CwdResolver;

  constructor(options?: {
    commandResolver?: CommandResolver;
    cwdResolver?: CwdResolver;
  }) {
    this.commandResolver = options?.commandResolver ?? resolveCommand;
    this.cwdResolver = options?.cwdResolver ?? resolveCwd;
  }

  /** Resolve a command name using the injected resolver (e.g. 'bun' → absolute path). */
  resolveCommand(command: string): string {
    return this.commandResolver(command);
  }

  /** Resolve a cwd path using the injected resolver (handles 'app:' prefix etc.). */
  resolveCwd(cwd: string | undefined, workspaceRootPath: string): string | undefined {
    return this.cwdResolver(cwd, workspaceRootPath);
  }

  /**
   * Build MCP server config from a source
   *
   * @param source - The source configuration
   * @param token - Authentication token (null for public/stdio sources)
   */
  buildMcpServer(source: LoadedSource, token: string | null): McpServerConfig | null {
    if (source.config.type !== 'mcp' || !source.config.mcp) {
      return null;
    }

    const mcp = source.config.mcp;

    // Handle stdio transport (local subprocess servers)
    if (mcp.transport === 'stdio') {
      if (!mcp.command) {
        debug(`[SourceServerBuilder] Stdio source ${source.config.slug} missing command`);
        return null;
      }
      const resolvedCommand = this.commandResolver(mcp.command);
      const resolvedCwd = this.cwdResolver(mcp.cwd, source.workspaceRootPath);

      // Claude CLI's Zod schema for stdio MCP servers does NOT include 'cwd',
      // so it gets stripped during validation. To work around this, resolve
      // relative args to absolute paths using the resolved cwd.
      let resolvedArgs = mcp.args;
      if (resolvedCwd && resolvedArgs?.length) {
        resolvedArgs = resolvedArgs.map((arg) =>
          arg.startsWith('/') || arg.startsWith('-') ? arg : join(resolvedCwd, arg)
        );
        debug(`[SourceServerBuilder] Resolved args with cwd: ${JSON.stringify(resolvedArgs)}`);
      }

      return {
        type: 'stdio',
        command: resolvedCommand,
        args: resolvedArgs,
        env: mcp.env,
      };
    }

    // Handle HTTP/SSE transport (remote servers)
    if (!mcp.url) {
      debug(`[SourceServerBuilder] HTTP/SSE source ${source.config.slug} missing URL`);
      return null;
    }

    const url = normalizeMcpUrl(mcp.url);

    const config: McpServerConfig = {
      type: url.includes('/sse') ? 'sse' : 'http',
      url,
    };

    // Handle authentication for HTTP/SSE
    // Note: Direct isAuthenticated check is safe here because we're inside authType !== 'none' block
    if (mcp.authType !== 'none') {
      if (token) {
        (config as { headers?: Record<string, string> }).headers = { Authorization: `Bearer ${token}` };
      } else if (source.config.isAuthenticated) {
        // Source claims to be authenticated but token is missing - needs re-auth
        debug(`[SourceServerBuilder] Source ${source.config.slug} needs re-authentication`);
        return null;
      }
    }

    return config;
  }

  /**
   * Build API server from a source
   *
   * @param source - The source configuration
   * @param credential - API credential (null for public APIs)
   * @param getToken - Token getter for OAuth APIs (Google, etc.) - supports auto-refresh
   * @param sessionPath - Optional path to session folder for saving large responses
   */
  async buildApiServer(
    source: LoadedSource,
    credential: ApiCredential | null,
    getToken?: () => Promise<string>,
    sessionPath?: string,
    summarize?: SummarizeCallback
  ): Promise<ReturnType<typeof createSdkMcpServer> | null> {
    if (source.config.type !== 'api') return null;
    if (!source.config.api) {
      debug(`[SourceServerBuilder] API source ${source.config.slug} missing api config`);
      return null;
    }

    const apiConfig = source.config.api;
    const authType = apiConfig.authType;
    const provider = source.config.provider;

    // Google APIs - use token getter with auto-refresh
    // Note: Direct isAuthenticated check is safe - Google OAuth always requires auth
    if (provider === 'google') {
      if (!source.config.isAuthenticated || !getToken) {
        debug(`[SourceServerBuilder] Google API source ${source.config.slug} not authenticated`);
        return null;
      }
      debug(`[SourceServerBuilder] Building Google API server for ${source.config.slug}`);
      const config = this.buildApiConfig(source);
      // Pass the token getter function - it will be called before each request
      // to get a fresh token (with auto-refresh if expired)
      return createApiServer(config, getToken, sessionPath, summarize);
    }

    // Slack APIs - use token getter with auto-refresh
    // Note: Direct isAuthenticated check is safe - Slack OAuth always requires auth
    if (provider === 'slack') {
      if (!source.config.isAuthenticated || !getToken) {
        debug(`[SourceServerBuilder] Slack API source ${source.config.slug} not authenticated`);
        return null;
      }
      debug(`[SourceServerBuilder] Building Slack API server for ${source.config.slug}`);
      const config = this.buildApiConfig(source);
      // Pass the token getter function - it will be called before each request
      // to get a fresh token (with auto-refresh if expired)
      return createApiServer(config, getToken, sessionPath, summarize);
    }

    // Public APIs (no auth) can be used immediately
    if (authType === 'none') {
      debug(`[SourceServerBuilder] Building public API server for ${source.config.slug}`);
      const config = this.buildApiConfig(source);
      return createApiServer(config, '', sessionPath, summarize);
    }

    // API key/bearer/header/query/basic auth - use static credential
    if (!credential) {
      debug(`[SourceServerBuilder] API source ${source.config.slug} needs credentials`);
      return null;
    }

    debug(`[SourceServerBuilder] Building API server for ${source.config.slug} (auth: ${authType})`);
    const config = this.buildApiConfig(source);
    return createApiServer(config, credential, sessionPath, summarize);
  }

  /**
   * Build ApiConfig from a LoadedSource
   */
  buildApiConfig(source: LoadedSource): ApiConfig {
    const api = source.config.api!;

    const config: ApiConfig = {
      name: source.config.slug,
      baseUrl: api.baseUrl,
      documentation: source.guide?.raw || '',
      defaultHeaders: api.defaultHeaders,
    };

    // Map auth type
    switch (api.authType) {
      case 'bearer':
        config.auth = { type: 'bearer', authScheme: api.authScheme ?? 'Bearer' };
        break;
      case 'header':
        config.auth = { type: 'header', headerName: api.headerName || 'x-api-key' };
        break;
      case 'query':
        config.auth = { type: 'query', queryParam: api.queryParam || 'api_key' };
        break;
      case 'basic':
        config.auth = { type: 'basic' };
        break;
      case 'none':
      default:
        config.auth = { type: 'none' };
    }

    return config;
  }

  /**
   * Build all MCP and API servers for enabled sources
   *
   * @param sourcesWithCredentials - Sources with their pre-loaded credentials
   * @param getTokenForSource - Function to get token getter for OAuth sources
   * @param sessionPath - Optional path to session folder for saving large API responses
   */
  async buildAll(
    sourcesWithCredentials: SourceWithCredential[],
    getTokenForSource?: (source: LoadedSource) => (() => Promise<string>) | undefined,
    sessionPath?: string,
    summarize?: SummarizeCallback
  ): Promise<BuiltServers> {
    const mcpServers: Record<string, McpServerConfig> = {};
    const apiServers: Record<string, ReturnType<typeof createSdkMcpServer>> = {};
    const errors: BuiltServers['errors'] = [];

    for (const { source, token, credential } of sourcesWithCredentials) {
      if (!isSourceUsable(source)) continue;

      try {
        if (source.config.type === 'mcp') {
          const config = this.buildMcpServer(source, token ?? null);
          if (config) {
            debug(`[SourceServerBuilder] Built MCP server for ${source.config.slug}`);
            mcpServers[source.config.slug] = config;
          } else if (source.config.mcp?.transport !== 'stdio' && source.config.mcp?.authType !== 'none') {
            // Only report auth error for HTTP/SSE sources that need auth
            // Stdio sources don't need auth
            debug(`[SourceServerBuilder] MCP server ${source.config.slug} needs auth`);
            errors.push({
              sourceSlug: source.config.slug,
              error: SERVER_BUILD_ERRORS.AUTH_REQUIRED,
            });
          }
        } else if (source.config.type === 'api') {
          const getToken = getTokenForSource?.(source);
          const server = await this.buildApiServer(source, credential ?? null, getToken, sessionPath, summarize);
          if (server) {
            apiServers[source.config.slug] = server;
          }
        }
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        debug(`[SourceServerBuilder] Failed to build server for ${source.config.slug}: ${message}`);
        errors.push({ sourceSlug: source.config.slug, error: message });
      }
    }

    return { mcpServers, apiServers, errors };
  }
}

/**
 * Normalize MCP URL to standard format
 * - Removes trailing slashes
 * - Preserves /sse suffix for SSE type detection
 * - Ensures /mcp suffix for HTTP type
 */
export function normalizeMcpUrl(url: string): string {
  url = url.replace(/\/+$/, '');

  // If URL ends with /sse, keep it for SSE type detection
  if (url.endsWith('/sse')) {
    return url;
  }

  // Ensure /mcp suffix for HTTP type
  if (!url.endsWith('/mcp')) {
    url = url + '/mcp';
  }

  return url;
}

// Singleton instance
let instance: SourceServerBuilder | null = null;

/**
 * 获取共享的 SourceServerBuilder 实例。
 * 首次调用时可传入 options 注入自定义解析器，后续调用返回同一实例。
 */
export function getSourceServerBuilder(options?: {
  commandResolver?: CommandResolver;
  cwdResolver?: CwdResolver;
}): SourceServerBuilder {
  if (!instance) {
    instance = new SourceServerBuilder(options);
  }
  return instance;
}
