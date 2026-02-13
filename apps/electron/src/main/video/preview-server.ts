/**
 * Preview Server
 *
 * Manages Remotion Studio preview server for real-time video preview.
 * Supports starting, stopping, and hot-reloading of video compositions.
 *
 * 使用内置 Bun 运行时启动 Remotion Studio，而非依赖外部 npx 命令。
 *
 * @requirements 3.1, 3.2, 3.3, 3.4, 3.5, 3.6
 */

import { join, dirname } from 'path';
import { existsSync } from 'fs';
import { ChildProcess, spawn } from 'child_process';
import log from '../logger';
import type { BunPathResolver } from '../bun-path';
import { getDefaultBunPathResolver } from '../bun-path';

const previewLog = log.scope('video:preview');

// ============================================================================
// Types and Interfaces
// ============================================================================

/**
 * Preview server configuration
 */
export interface PreviewServerConfig {
  /** Port to run the preview server on */
  port?: number;
  /** Host to bind the server to */
  host?: string;
  /** Startup timeout in milliseconds */
  startupTimeout?: number;
}

/**
 * Preview server result
 */
export interface PreviewServerResult {
  /** URL to access the preview */
  url: string;
  /** Port the server is running on */
  port: number;
  /** Process ID of the server */
  pid: number;
}

/**
 * Interface for Preview Server
 */
export interface IPreviewServer {
  start(projectPath: string, config?: PreviewServerConfig): Promise<PreviewServerResult>;
  stop(projectPath: string): Promise<void>;
  stopAll(): Promise<void>;
  isRunning(projectPath: string): boolean;
  getUrl(projectPath: string): string | null;
}

/**
 * Preview server error types
 */
export enum PreviewServerError {
  /** Project path not found */
  INVALID_PROJECT_PATH = 'INVALID_PROJECT_PATH',
  /** Failed to start preview server */
  PREVIEW_START_FAILED = 'PREVIEW_START_FAILED',
  /** Port is already in use */
  PREVIEW_PORT_IN_USE = 'PREVIEW_PORT_IN_USE',
  /** Server startup timeout */
  PREVIEW_TIMEOUT = 'PREVIEW_TIMEOUT',
  /** Bun runtime not found */
  BUN_NOT_FOUND = 'BUN_NOT_FOUND',
}

/**
 * Preview server error class
 */
export class PreviewServerException extends Error {
  constructor(
    message: string,
    public readonly code: PreviewServerError,
    public readonly projectPath?: string
  ) {
    super(message);
    this.name = 'PreviewServerException';
  }
}

// ============================================================================
// Preview Server Implementation
// ============================================================================

/**
 * Preview Server class for managing Remotion Studio instances
 *
 * 使用 Bun 子进程启动 Remotion Studio，支持多项目并发预览。
 *
 * @requirements 3.1, 3.2, 3.3, 3.4, 3.5, 3.6
 */
export class PreviewServer implements IPreviewServer {
  /** Map of project paths to their preview server processes */
  private servers: Map<string, {
    process: ChildProcess;
    port: number;
    url: string;
    pid: number;
  }> = new Map();

  /** Default port for preview servers */
  private defaultPort: number = 3100;

  /** Port counter for multiple servers */
  private portCounter: number = 0;

  /** Bun path resolver for getting Bun executable path */
  private bunResolver: BunPathResolver;

  /**
   * Create a new PreviewServer instance
   *
   * @param bunResolver - Optional BunPathResolver instance for dependency injection
   */
  constructor(bunResolver?: BunPathResolver) {
    this.bunResolver = bunResolver ?? getDefaultBunPathResolver();
    previewLog.info('PreviewServer initialized with Bun runtime support');

    // Clean up servers on process exit
    process.on('exit', () => {
      this.stopAll().catch(() => {});
    });
  }

  /**
   * Start a preview server for a video project
   *
   * @param projectPath - Path to the video project
   * @param config - Optional server configuration
   * @returns Preview server result with URL and PID
   *
   * @requirements 3.1, 3.2, 3.3
   */
  async start(
    projectPath: string,
    config?: PreviewServerConfig
  ): Promise<PreviewServerResult> {
    // Check if server is already running for this project
    const existing = this.servers.get(projectPath);
    if (existing) {
      previewLog.info(`Preview server already running for ${projectPath} at ${existing.url}`);
      return {
        url: existing.url,
        port: existing.port,
        pid: existing.pid,
      };
    }

    // Validate project path
    if (!existsSync(projectPath)) {
      throw new PreviewServerException(
        `Project path not found: ${projectPath}`,
        PreviewServerError.INVALID_PROJECT_PATH,
        projectPath
      );
    }

    // Determine port
    const port = config?.port || this.getNextPort();
    const host = config?.host || 'localhost';
    const startupTimeout = config?.startupTimeout || 30000;

    previewLog.info(`Starting preview server for ${projectPath} on port ${port}`);

    try {
      // Find the entry point for the project
      const entryPoint = this.findEntryPoint(projectPath);
      previewLog.info(`Using entry point: ${entryPoint}`);

      // Start Remotion Studio using Bun
      const serverProcess = await this.startRemotionStudio(entryPoint, port, host, startupTimeout);

      // Wait for server to be ready
      const url = `http://${host}:${port}`;
      await this.waitForServerReady(url, startupTimeout);

      const pid = serverProcess.pid!;

      // Store server info
      this.servers.set(projectPath, {
        process: serverProcess,
        port,
        url,
        pid,
      });

      previewLog.info(`Preview server started at ${url} (PID: ${pid})`);

      return { url, port, pid };
    } catch (error) {
      previewLog.error(`Failed to start preview server:`, error);

      // Re-throw PreviewServerException as-is
      if (error instanceof PreviewServerException) {
        throw error;
      }

      // Wrap other errors
      const message = error instanceof Error ? error.message : String(error);
      throw new PreviewServerException(
        `Failed to start preview server: ${message}`,
        PreviewServerError.PREVIEW_START_FAILED,
        projectPath
      );
    }
  }

  /**
   * Stop a preview server for a specific project
   *
   * @param projectPath - Path to the video project
   *
   * @requirements 3.6
   */
  async stop(projectPath: string): Promise<void> {
    const server = this.servers.get(projectPath);
    if (!server) {
      previewLog.warn(`No preview server running for ${projectPath}`);
      return;
    }

    previewLog.info(`Stopping preview server for ${projectPath} (PID: ${server.pid})`);

    try {
      // Kill the server process
      await this.killProcess(server.process);

      // Remove from map
      this.servers.delete(projectPath);

      previewLog.info(`Preview server stopped for ${projectPath}`);
    } catch (error) {
      previewLog.error(`Error stopping preview server:`, error);
      // Still remove from map even if kill fails
      this.servers.delete(projectPath);
    }
  }

  /**
   * Stop all running preview servers
   *
   * @requirements 3.6
   */
  async stopAll(): Promise<void> {
    previewLog.info(`Stopping all preview servers (${this.servers.size} running)`);

    const stopPromises: Promise<void>[] = [];
    for (const projectPath of this.servers.keys()) {
      stopPromises.push(this.stop(projectPath));
    }

    await Promise.allSettled(stopPromises);
    previewLog.info('All preview servers stopped');
  }

  /**
   * Check if a preview server is running for a project
   *
   * @requirements 3.5
   */
  isRunning(projectPath: string): boolean {
    return this.servers.has(projectPath);
  }

  /**
   * Get the URL for a running preview server
   */
  getUrl(projectPath: string): string | null {
    const server = this.servers.get(projectPath);
    return server?.url || null;
  }

  // ============================================================================
  // Private Helper Methods
  // ============================================================================

  /**
   * Get the next available port
   */
  private getNextPort(): number {
    const port = this.defaultPort + this.portCounter;
    this.portCounter++;
    return port;
  }

  /**
   * Find the entry point for a video project
   */
  private findEntryPoint(projectPath: string): string {
    // Try common entry point locations
    const possibleEntryPoints = [
      join(projectPath, 'src', 'Root.tsx'),
      join(projectPath, 'src', 'index.tsx'),
      join(projectPath, 'Root.tsx'),
      join(projectPath, 'index.tsx'),
    ];

    for (const entryPoint of possibleEntryPoints) {
      if (existsSync(entryPoint)) {
        return entryPoint;
      }
    }

    // If no entry point found in project, use the video package's Root
    try {
      const videoPackageRoot = require.resolve('@sprouty-ai/video');
      const videoPackageDir = dirname(videoPackageRoot);
      const defaultEntryPoint = join(videoPackageDir, 'Root.tsx');
      if (existsSync(defaultEntryPoint)) {
        return defaultEntryPoint;
      }
    } catch {
      // Ignore if video package not found
    }

    // Fallback to project path itself
    return projectPath;
  }

  /**
   * Find the Remotion CLI path
   *
   * 查找 @remotion/cli 的路径，优先使用项目本地安装的版本
   */
  private findRemotionCliPath(projectPath: string): string {
    // Try to find @remotion/cli in the project's node_modules
    const possiblePaths = [
      // Project local installation
      join(projectPath, 'node_modules', '@remotion', 'cli', 'dist', 'cli.js'),
      // Monorepo root installation (relative to project)
      join(projectPath, '..', '..', 'node_modules', '@remotion', 'cli', 'dist', 'cli.js'),
      // Try to resolve from current process
      join(process.cwd(), 'node_modules', '@remotion', 'cli', 'dist', 'cli.js'),
    ];

    for (const cliPath of possiblePaths) {
      if (existsSync(cliPath)) {
        return cliPath;
      }
    }

    // Try to resolve using require.resolve
    try {
      const resolved = require.resolve('@remotion/cli/dist/cli.js');
      if (existsSync(resolved)) {
        return resolved;
      }
    } catch {
      // Ignore resolution errors
    }

    // Fallback: assume it's in the project's node_modules
    return join(projectPath, 'node_modules', '@remotion', 'cli', 'dist', 'cli.js');
  }

  /**
   * Start Remotion Studio process using Bun
   *
   * @requirements 3.1, 3.2 - 使用内置 Bun 运行时启动 Remotion Studio
   */
  private startRemotionStudio(
    entryPoint: string,
    port: number,
    host: string,
    startupTimeout: number
  ): Promise<ChildProcess> {
    return new Promise((resolve, reject) => {
      // Get Bun executable path
      let bunPath: string;
      try {
        bunPath = this.bunResolver.getBunPath();
        previewLog.debug(`Using Bun at: ${bunPath}`);
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        reject(new PreviewServerException(
          `Bun runtime not found: ${message}`,
          PreviewServerError.BUN_NOT_FOUND
        ));
        return;
      }

      // Find Remotion CLI path
      const projectDir = dirname(entryPoint);
      const remotionCliPath = this.findRemotionCliPath(projectDir);
      previewLog.debug(`Using Remotion CLI at: ${remotionCliPath}`);

      // Build arguments for Bun to run Remotion CLI
      // Command: bun run <remotion-cli-path> studio <entry-point> --port <port> --log error
      const args = [
        'run',
        remotionCliPath,
        'studio',
        entryPoint,
        '--port',
        port.toString(),
        '--log',
        'error', // Reduce log verbosity
      ];

      previewLog.debug(`Spawning: ${bunPath} ${args.join(' ')}`);

      const serverProcess = spawn(bunPath, args, {
        cwd: projectDir,
        stdio: ['ignore', 'pipe', 'pipe'],
        detached: false,
        env: {
          ...process.env,
          // Ensure Bun uses the correct Node.js compatibility mode
          BUN_RUNTIME_TRANSPILER_CACHE_PATH: join(projectDir, '.bun-cache'),
        },
      });

      let started = false;
      let errorOutput = '';

      // Handle stdout
      serverProcess.stdout?.on('data', (data: Buffer) => {
        const output = data.toString();
        previewLog.debug(`[studio stdout] ${output}`);

        // Check if server has started
        if (!started && (output.includes('Server ready') || output.includes('http://') || output.includes('Listening'))) {
          started = true;
          resolve(serverProcess);
        }
      });

      // Handle stderr
      serverProcess.stderr?.on('data', (data: Buffer) => {
        const output = data.toString();
        errorOutput += output;
        previewLog.debug(`[studio stderr] ${output}`);

        // Some startup messages go to stderr
        if (!started && (output.includes('Server ready') || output.includes('http://') || output.includes('Listening'))) {
          started = true;
          resolve(serverProcess);
        }

        // Check for port in use error
        if (output.includes('EADDRINUSE') || output.includes('address already in use')) {
          reject(new PreviewServerException(
            `Port ${port} is already in use`,
            PreviewServerError.PREVIEW_PORT_IN_USE
          ));
        }
      });

      // Handle process error
      serverProcess.on('error', (error) => {
        previewLog.error('Failed to start Remotion Studio:', error);
        reject(new PreviewServerException(
          `Failed to start Remotion Studio: ${error.message}`,
          PreviewServerError.PREVIEW_START_FAILED
        ));
      });

      // Handle process exit
      serverProcess.on('exit', (code) => {
        if (!started) {
          reject(new PreviewServerException(
            `Remotion Studio exited with code ${code}: ${errorOutput}`,
            PreviewServerError.PREVIEW_START_FAILED
          ));
        }
      });

      // Timeout for startup
      const timeoutMs = Math.min(startupTimeout, 15000); // Use shorter timeout for initial startup detection
      setTimeout(() => {
        if (!started) {
          // Assume server started if process is still running
          if (serverProcess.exitCode === null && serverProcess.pid) {
            started = true;
            previewLog.info(`Assuming server started (process still running, PID: ${serverProcess.pid})`);
            resolve(serverProcess);
          } else {
            reject(new PreviewServerException(
              'Remotion Studio startup timeout',
              PreviewServerError.PREVIEW_TIMEOUT
            ));
          }
        }
      }, timeoutMs);
    });
  }

  /**
   * Wait for the server to be ready by polling the URL
   *
   * @requirements 3.3 - Return preview URL when server is ready
   */
  private async waitForServerReady(url: string, timeout: number): Promise<void> {
    const startTime = Date.now();
    const pollInterval = 500; // 500ms between polls

    while (Date.now() - startTime < timeout) {
      try {
        // Try to fetch the URL
        const response = await fetch(url, {
          method: 'HEAD',
          signal: AbortSignal.timeout(2000),
        });

        if (response.ok || response.status === 200) {
          previewLog.debug(`Server ready at ${url}`);
          return;
        }
      } catch {
        // Server not ready yet, continue polling
      }

      // Wait before next poll
      await new Promise((resolve) => setTimeout(resolve, pollInterval));
    }

    // If we get here, assume server is ready (some servers don't respond to HEAD)
    previewLog.warn(`Server readiness check timed out, assuming ready at ${url}`);
  }

  /**
   * Kill a server process
   *
   * @requirements 3.6 - Properly clean up resources
   */
  private async killProcess(childProcess: ChildProcess): Promise<void> {
    return new Promise((resolve) => {
      try {
        if (!childProcess.pid) {
          resolve();
          return;
        }

        const pid = childProcess.pid;

        // On Windows, use taskkill
        if (process.platform === 'win32') {
          const taskkill = spawn('taskkill', ['/pid', pid.toString(), '/f', '/t'], {
            stdio: 'ignore',
          });
          taskkill.on('close', () => resolve());
          taskkill.on('error', () => resolve());
        } else {
          // On Unix, send SIGTERM first
          childProcess.kill('SIGTERM');

          // Set up a timeout for SIGKILL
          const killTimeout = setTimeout(() => {
            try {
              childProcess.kill('SIGKILL');
            } catch {
              // Process may already be dead
            }
            resolve();
          }, 5000);

          // Listen for process exit
          childProcess.on('exit', () => {
            clearTimeout(killTimeout);
            resolve();
          });

          // Also resolve after a short delay if no exit event
          setTimeout(() => {
            clearTimeout(killTimeout);
            resolve();
          }, 6000);
        }
      } catch (error) {
        previewLog.warn('Error killing process:', error);
        resolve();
      }
    });
  }
}

/**
 * Create a new PreviewServer instance
 *
 * @param bunResolver - Optional BunPathResolver for dependency injection
 * @returns PreviewServer instance
 *
 * @example
 * ```typescript
 * import { createPreviewServer } from './preview-server';
 * import { createBunPathResolver } from '../bun-path';
 *
 * // Use default Bun resolver
 * const server = createPreviewServer();
 *
 * // Or inject custom resolver
 * const customResolver = createBunPathResolver({ isPackaged: false });
 * const server = createPreviewServer(customResolver);
 *
 * // Start preview
 * const result = await server.start('/path/to/project');
 * console.log(`Preview at ${result.url} (PID: ${result.pid})`);
 * ```
 */
export function createPreviewServer(bunResolver?: BunPathResolver): PreviewServer {
  return new PreviewServer(bunResolver);
}
