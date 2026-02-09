/**
 * Video Service Manager
 *
 * 管理 MCP Video Server 和其他视频相关子进程的生命周期。
 * 使用内置 Bun 运行时启动服务，支持启动、停止、重启和状态跟踪。
 *
 * @requirements 2.1, 2.2, 2.3, 2.4, 2.5, 2.6
 */

import { ChildProcess, spawn } from 'child_process';
import { EventEmitter } from 'events';
import log from '../logger';
import type { BunPathResolver } from '../bun-path';

const serviceLog = log.scope('video:service');

// ============================================================================
// Types and Interfaces
// ============================================================================

/**
 * 视频服务配置
 */
export interface VideoServiceConfig {
  /** MCP Video Server 入口文件路径 */
  mcpServerEntry: string;
  /** 最大重试次数（默认 3） */
  maxRetries?: number;
  /** 重试间隔（毫秒，默认 2000） */
  retryInterval?: number;
  /** 启动超时（毫秒，默认 30000） */
  startupTimeout?: number;
  /** 停止超时（毫秒，默认 5000） */
  stopTimeout?: number;
}

/**
 * 服务状态
 */
export interface ServiceStatus {
  /** 服务是否运行中 */
  running: boolean;
  /** 进程 ID */
  pid?: number;
  /** 启动时间 */
  startedAt?: number;
  /** 重启次数 */
  restartCount: number;
  /** 最后错误 */
  lastError?: string;
}

/**
 * 视频服务管理器接口
 */
export interface VideoServiceManager {
  /**
   * 启动 MCP Video Server
   * @throws Error 如果启动失败
   */
  startMcpServer(): Promise<void>;

  /**
   * 停止 MCP Video Server
   */
  stopMcpServer(): Promise<void>;

  /**
   * 重启 MCP Video Server
   */
  restartMcpServer(): Promise<void>;

  /**
   * 获取服务状态
   */
  getStatus(): ServiceStatus;

  /**
   * 获取 MCP 服务器进程（用于 stdio 通信）
   */
  getMcpProcess(): ChildProcess | null;

  /**
   * 停止所有服务
   */
  stopAll(): Promise<void>;

  /**
   * 添加事件监听器
   */
  on(event: 'status-changed', listener: (status: ServiceStatus) => void): void;
  on(event: 'error', listener: (error: Error) => void): void;
  on(event: 'log', listener: (level: 'info' | 'error', message: string) => void): void;
}

/**
 * 视频服务错误类型
 */
export enum VideoServiceError {
  // Bun 相关错误
  BUN_NOT_FOUND = 'BUN_NOT_FOUND',
  BUN_EXECUTION_FAILED = 'BUN_EXECUTION_FAILED',

  // MCP 服务器错误
  MCP_SERVER_START_FAILED = 'MCP_SERVER_START_FAILED',
  MCP_SERVER_CRASHED = 'MCP_SERVER_CRASHED',
  MCP_SERVER_TIMEOUT = 'MCP_SERVER_TIMEOUT',

  // 通用错误
  MAX_RETRIES_EXCEEDED = 'MAX_RETRIES_EXCEEDED',
  PROCESS_COMMUNICATION_FAILED = 'PROCESS_COMMUNICATION_FAILED',
}

/**
 * 视频服务错误
 */
export class VideoServiceManagerError extends Error {
  constructor(
    message: string,
    public readonly code: VideoServiceError,
    public readonly details?: string
  ) {
    super(message);
    this.name = 'VideoServiceManagerError';
  }
}

// ============================================================================
// Default Configuration
// ============================================================================

const DEFAULT_CONFIG: Required<Omit<VideoServiceConfig, 'mcpServerEntry'>> = {
  maxRetries: 3,
  retryInterval: 2000,
  startupTimeout: 30000,
  stopTimeout: 5000,
};

// ============================================================================
// Video Service Manager Implementation
// ============================================================================

/**
 * 创建视频服务管理器
 *
 * @param bunResolver - Bun 路径解析器
 * @param config - 服务配置
 * @returns VideoServiceManager 实例
 *
 * @example
 * ```typescript
 * import { createBunPathResolver } from '../bun-path';
 * import { createVideoServiceManager } from './service-manager';
 *
 * const bunResolver = createBunPathResolver();
 * const manager = createVideoServiceManager(bunResolver, {
 *   mcpServerEntry: '/path/to/mcp-video/src/index.ts',
 * });
 *
 * await manager.startMcpServer();
 * const status = manager.getStatus();
 * console.log('MCP Server running:', status.running);
 * ```
 */
export function createVideoServiceManager(
  bunResolver: BunPathResolver,
  config: VideoServiceConfig
): VideoServiceManager {
  // 合并默认配置
  const fullConfig: Required<VideoServiceConfig> = {
    ...DEFAULT_CONFIG,
    ...config,
  };

  // 事件发射器
  const emitter = new EventEmitter();

  // MCP 服务器进程
  let mcpProcess: ChildProcess | null = null;

  // 服务状态
  let status: ServiceStatus = {
    running: false,
    restartCount: 0,
  };

  // 是否正在停止（防止重复停止）
  let isStopping = false;

  // 是否正在启动（防止重复启动）
  let isStarting = false;

  /**
   * 更新服务状态并发出事件
   */
  function updateStatus(updates: Partial<ServiceStatus>): void {
    status = { ...status, ...updates };
    emitter.emit('status-changed', status);
    serviceLog.debug('Status updated:', status);
  }

  /**
   * 发出日志事件
   */
  function emitLog(level: 'info' | 'error', message: string): void {
    emitter.emit('log', level, message);
    if (level === 'error') {
      serviceLog.error(message);
    } else {
      serviceLog.info(message);
    }
  }

  /**
   * 发出错误事件
   */
  function emitError(error: Error): void {
    emitter.emit('error', error);
    serviceLog.error('Service error:', error);
  }

  /**
   * 等待指定时间
   */
  function delay(ms: number): Promise<void> {
    return new Promise((resolve) => setTimeout(resolve, ms));
  }

  /**
   * 启动 MCP Video Server 进程
   * @requirements 2.1 - 使用内置 Bun 运行时启动 MCP Video Server
   */
  async function spawnMcpServer(): Promise<ChildProcess> {
    // 获取 Bun 路径
    const bunPath = bunResolver.getBunPath();
    serviceLog.info(`Using Bun at: ${bunPath}`);

    // 启动参数
    const args = ['run', fullConfig.mcpServerEntry];

    serviceLog.info(`Spawning MCP Video Server: ${bunPath} ${args.join(' ')}`);

    // 创建子进程
    // @requirements 6.1 - 通过 stdio 管道与 MCP Video Server 通信
    const childProcess = spawn(bunPath, args, {
      stdio: ['pipe', 'pipe', 'pipe'],
      detached: false,
      env: {
        ...process.env,
        // 确保子进程使用正确的 NODE_ENV
        NODE_ENV: process.env.NODE_ENV || 'production',
      },
    });

    return childProcess;
  }

  /**
   * 设置进程事件处理器
   * @requirements 6.2, 6.3 - 转发日志和捕获错误
   */
  function setupProcessHandlers(childProcess: ChildProcess): void {
    // 处理 stdout
    // @requirements 6.2 - 将 stdout 转发到主进程日志系统
    childProcess.stdout?.on('data', (data: Buffer) => {
      const output = data.toString().trim();
      if (output) {
        emitLog('info', `[MCP stdout] ${output}`);
      }
    });

    // 处理 stderr
    // @requirements 6.3 - 捕获 stderr 并记录错误
    childProcess.stderr?.on('data', (data: Buffer) => {
      const output = data.toString().trim();
      if (output) {
        emitLog('error', `[MCP stderr] ${output}`);
      }
    });

    // 处理进程错误
    childProcess.on('error', (error) => {
      serviceLog.error('MCP Server process error:', error);
      updateStatus({
        running: false,
        lastError: error.message,
      });
      emitError(
        new VideoServiceManagerError(
          `MCP Server process error: ${error.message}`,
          VideoServiceError.MCP_SERVER_CRASHED,
          error.stack
        )
      );
    });

    // 处理进程退出
    // @requirements 2.5 - 记录错误并支持重启
    childProcess.on('exit', (code, signal) => {
      serviceLog.info(`MCP Server exited with code ${code}, signal ${signal}`);

      // 如果不是主动停止，则记录为意外退出
      if (!isStopping && status.running) {
        const errorMessage = `MCP Server unexpectedly exited with code ${code}`;
        updateStatus({
          running: false,
          pid: undefined,
          lastError: errorMessage,
        });
        emitError(
          new VideoServiceManagerError(
            errorMessage,
            VideoServiceError.MCP_SERVER_CRASHED
          )
        );
      } else {
        updateStatus({
          running: false,
          pid: undefined,
        });
      }

      mcpProcess = null;
    });

    // 处理进程关闭
    childProcess.on('close', (code, signal) => {
      serviceLog.debug(`MCP Server closed with code ${code}, signal ${signal}`);
    });
  }

  /**
   * 等待服务器就绪
   */
  async function waitForServerReady(
    childProcess: ChildProcess,
    timeout: number
  ): Promise<void> {
    return new Promise((resolve, reject) => {
      const timeoutId = setTimeout(() => {
        reject(
          new VideoServiceManagerError(
            'MCP Server startup timeout',
            VideoServiceError.MCP_SERVER_TIMEOUT
          )
        );
      }, timeout);

      let resolved = false;

      // 监听 stdout 检测服务器就绪
      const onData = (data: Buffer) => {
        const output = data.toString();
        // 检测服务器就绪信号（根据 MCP Video Server 的输出格式调整）
        if (
          output.includes('Server started') ||
          output.includes('ready') ||
          output.includes('listening') ||
          output.includes('MCP server running')
        ) {
          if (!resolved) {
            resolved = true;
            clearTimeout(timeoutId);
            childProcess.stdout?.off('data', onData);
            resolve();
          }
        }
      };

      childProcess.stdout?.on('data', onData);

      // 如果进程提前退出，则拒绝
      childProcess.once('exit', (code) => {
        if (!resolved) {
          resolved = true;
          clearTimeout(timeoutId);
          childProcess.stdout?.off('data', onData);
          reject(
            new VideoServiceManagerError(
              `MCP Server exited during startup with code ${code}`,
              VideoServiceError.MCP_SERVER_START_FAILED
            )
          );
        }
      });

      // 如果进程出错，则拒绝
      childProcess.once('error', (error) => {
        if (!resolved) {
          resolved = true;
          clearTimeout(timeoutId);
          childProcess.stdout?.off('data', onData);
          reject(
            new VideoServiceManagerError(
              `MCP Server error during startup: ${error.message}`,
              VideoServiceError.MCP_SERVER_START_FAILED,
              error.stack
            )
          );
        }
      });

      // 如果 5 秒内没有检测到就绪信号但进程仍在运行，假设已就绪
      setTimeout(() => {
        if (!resolved && childProcess.exitCode === null) {
          resolved = true;
          clearTimeout(timeoutId);
          childProcess.stdout?.off('data', onData);
          serviceLog.warn('No ready signal detected, assuming server is ready');
          resolve();
        }
      }, 5000);
    });
  }

  /**
   * 终止进程
   * @requirements 6.4, 6.5 - 支持 SIGTERM 和 SIGKILL 信号
   */
  async function killProcess(
    childProcess: ChildProcess,
    timeout: number
  ): Promise<void> {
    return new Promise((resolve) => {
      if (!childProcess.pid) {
        resolve();
        return;
      }

      const pid = childProcess.pid;
      let killed = false;

      // 设置超时强制终止
      // @requirements 6.5 - 在 SIGTERM 超时后发送 SIGKILL
      const forceKillTimeout = setTimeout(() => {
        if (!killed) {
          serviceLog.warn(`Process ${pid} did not respond to SIGTERM, sending SIGKILL`);
          try {
            if (process.platform === 'win32') {
              // Windows: 使用 taskkill 强制终止进程树
              spawn('taskkill', ['/pid', pid.toString(), '/f', '/t'], {
                stdio: 'ignore',
              });
            } else {
              // Unix: 发送 SIGKILL
              childProcess.kill('SIGKILL');
            }
          } catch (error) {
            serviceLog.warn('Error sending SIGKILL:', error);
          }
        }
      }, timeout);

      // 监听进程退出
      childProcess.once('exit', () => {
        killed = true;
        clearTimeout(forceKillTimeout);
        resolve();
      });

      // 发送 SIGTERM
      // @requirements 6.4 - 支持向子进程发送信号
      try {
        serviceLog.info(`Sending SIGTERM to process ${pid}`);
        if (process.platform === 'win32') {
          // Windows: 使用 taskkill
          spawn('taskkill', ['/pid', pid.toString(), '/t'], {
            stdio: 'ignore',
          });
        } else {
          // Unix: 发送 SIGTERM
          childProcess.kill('SIGTERM');
        }
      } catch (error) {
        serviceLog.warn('Error sending SIGTERM:', error);
        clearTimeout(forceKillTimeout);
        resolve();
      }
    });
  }

  // ============================================================================
  // Public API
  // ============================================================================

  return {
    /**
     * 启动 MCP Video Server
     * @requirements 2.1, 2.2, 2.3
     */
    async startMcpServer(): Promise<void> {
      // 防止重复启动
      if (isStarting) {
        serviceLog.warn('MCP Server is already starting');
        return;
      }

      // 如果已经运行，直接返回
      if (status.running && mcpProcess) {
        serviceLog.info('MCP Server is already running');
        return;
      }

      isStarting = true;
      let lastError: Error | null = null;

      try {
        // 重试逻辑
        // @requirements 7.4 - 支持配置最大重试次数
        for (let attempt = 0; attempt <= fullConfig.maxRetries; attempt++) {
          if (attempt > 0) {
            serviceLog.info(
              `Retrying MCP Server start (attempt ${attempt}/${fullConfig.maxRetries})...`
            );
            await delay(fullConfig.retryInterval);
          }

          try {
            // 启动进程
            const childProcess = await spawnMcpServer();
            mcpProcess = childProcess;

            // 设置事件处理器
            setupProcessHandlers(childProcess);

            // 等待服务器就绪
            await waitForServerReady(childProcess, fullConfig.startupTimeout);

            // 更新状态
            // @requirements 2.2 - 管理生命周期
            updateStatus({
              running: true,
              pid: childProcess.pid,
              startedAt: Date.now(),
              lastError: undefined,
            });

            emitLog('info', `MCP Video Server started successfully (PID: ${childProcess.pid})`);
            return;
          } catch (error) {
            lastError = error instanceof Error ? error : new Error(String(error));
            serviceLog.error(`Start attempt ${attempt + 1} failed:`, lastError);

            // 清理失败的进程
            if (mcpProcess) {
              try {
                await killProcess(mcpProcess, 2000);
              } catch {
                // 忽略清理错误
              }
              mcpProcess = null;
            }

            updateStatus({
              running: false,
              lastError: lastError.message,
              restartCount: status.restartCount + 1,
            });
          }
        }

        // 所有重试都失败
        // @requirements 7.1 - 记录详细错误并通知用户
        throw new VideoServiceManagerError(
          `Failed to start MCP Server after ${fullConfig.maxRetries + 1} attempts`,
          VideoServiceError.MAX_RETRIES_EXCEEDED,
          lastError?.message
        );
      } finally {
        isStarting = false;
      }
    },

    /**
     * 停止 MCP Video Server
     * @requirements 2.2, 2.4
     */
    async stopMcpServer(): Promise<void> {
      // 防止重复停止
      if (isStopping) {
        serviceLog.warn('MCP Server is already stopping');
        return;
      }

      // 如果没有运行，直接返回
      if (!mcpProcess) {
        serviceLog.info('MCP Server is not running');
        updateStatus({ running: false });
        return;
      }

      isStopping = true;

      try {
        serviceLog.info('Stopping MCP Video Server...');

        // 终止进程
        // @requirements 2.4 - 优雅地停止所有视频服务子进程
        await killProcess(mcpProcess, fullConfig.stopTimeout);

        mcpProcess = null;
        updateStatus({
          running: false,
          pid: undefined,
        });

        emitLog('info', 'MCP Video Server stopped successfully');
      } catch (error) {
        const errorMessage = error instanceof Error ? error.message : String(error);
        serviceLog.error('Error stopping MCP Server:', error);
        updateStatus({
          running: false,
          lastError: errorMessage,
        });
      } finally {
        isStopping = false;
      }
    },

    /**
     * 重启 MCP Video Server
     * @requirements 2.2, 2.5, 7.5
     */
    async restartMcpServer(): Promise<void> {
      serviceLog.info('Restarting MCP Video Server...');

      // 先停止
      await this.stopMcpServer();

      // 等待一小段时间确保资源释放
      await delay(500);

      // 再启动
      // @requirements 7.5 - 重启时保持之前的配置状态
      await this.startMcpServer();

      emitLog('info', 'MCP Video Server restarted successfully');
    },

    /**
     * 获取服务状态
     * @requirements 2.6 - 提供获取连接信息的方法
     */
    getStatus(): ServiceStatus {
      return { ...status };
    },

    /**
     * 获取 MCP 服务器进程
     * @requirements 2.6 - 用于 stdio 通信
     */
    getMcpProcess(): ChildProcess | null {
      return mcpProcess;
    },

    /**
     * 停止所有服务
     * @requirements 2.4 - 优雅地停止所有视频服务子进程
     */
    async stopAll(): Promise<void> {
      serviceLog.info('Stopping all video services...');
      await this.stopMcpServer();
      emitLog('info', 'All video services stopped');
    },

    /**
     * 添加事件监听器
     */
    on(event: string, listener: (...args: unknown[]) => void): void {
      emitter.on(event, listener);
    },
  };
}
