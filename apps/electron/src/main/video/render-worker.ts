/**
 * Render Worker
 *
 * Handles video rendering using Bun subprocess to call Remotion's bundler and renderer APIs.
 * Supports multiple quality presets and output formats with progress reporting.
 *
 * 重构说明：
 * - 使用 Bun 子进程调用独立的渲染脚本 (apps/mcp-video/src/render-script.ts)
 * - 通过解析子进程 stdout 获取 JSON 格式的进度信息
 * - 支持通过 SIGTERM 信号取消渲染
 *
 * @requirements 4.1, 4.2, 4.3, 4.4, 4.5
 */

import { spawn, type ChildProcess } from 'child_process';
import { join, dirname } from 'path';
import { existsSync, mkdirSync } from 'fs';
import { app } from 'electron';
import type {
  RenderProgress,
  OutputFormat,
  QualityPreset,
} from '@sprouty-ai/video';
import type { BunPathResolver } from '../bun-path';
import { createBunPathResolver } from '../bun-path';
import log from '../logger';

const renderLog = log.scope('video:render');

// ============================================================================
// Types and Interfaces
// ============================================================================

/**
 * 渲染脚本支持的质量预设
 * 渲染脚本支持: draft, preview, production, high
 */
export type RenderScriptQualityPreset = 'draft' | 'preview' | 'production' | 'high';

/**
 * 将 @sprouty-ai/video 的 QualityPreset 映射到渲染脚本的质量预设
 */
function mapQualityPreset(quality: QualityPreset): RenderScriptQualityPreset {
  switch (quality) {
    case 'draft':
      return 'draft';
    case 'standard':
      return 'production';
    case 'high':
      return 'high';
    default:
      return 'production';
  }
}

/**
 * Render error types for categorizing failures
 * @requirements 4.5 - 返回分类的错误信息
 */
export enum RenderErrorType {
  BUNDLE_FAILED = 'BUNDLE_FAILED',
  COMPOSITION_NOT_FOUND = 'COMPOSITION_NOT_FOUND',
  RENDER_FAILED = 'RENDER_FAILED',
  OUTPUT_WRITE_FAILED = 'OUTPUT_WRITE_FAILED',
  CANCELLED = 'CANCELLED',
  INVALID_PROJECT = 'INVALID_PROJECT',
  BUN_NOT_FOUND = 'BUN_NOT_FOUND',
  SCRIPT_NOT_FOUND = 'SCRIPT_NOT_FOUND',
}

/**
 * Render error structure
 */
export interface RenderError {
  type: RenderErrorType;
  message: string;
  details?: string;
  stack?: string;
}

/**
 * Options for the render worker
 */
export interface RenderWorkerOptions {
  /** Path to the video project directory */
  projectPath: string;
  /** ID of the composition to render */
  compositionId: string;
  /** Output file path */
  outputPath: string;
  /** Quality preset (draft, standard, high) - will be mapped to render script presets */
  quality: QualityPreset;
  /** Output format */
  outputFormat?: OutputFormat;
  /** Progress callback */
  onProgress?: (progress: RenderProgress) => void;
}

/**
 * Interface for Render Worker
 */
export interface IRenderWorker {
  render(options: RenderWorkerOptions): Promise<string>;
  cancel(): void;
}

/**
 * 渲染脚本退出码映射
 */
const EXIT_CODE_TO_ERROR_TYPE: Record<number, RenderErrorType> = {
  1: RenderErrorType.INVALID_PROJECT, // INVALID_ARGS
  2: RenderErrorType.INVALID_PROJECT, // PROJECT_NOT_FOUND
  3: RenderErrorType.COMPOSITION_NOT_FOUND,
  4: RenderErrorType.BUNDLE_FAILED,
  5: RenderErrorType.RENDER_FAILED,
  6: RenderErrorType.OUTPUT_WRITE_FAILED,
};

// ============================================================================
// Render Worker Implementation
// ============================================================================

/**
 * Render Worker class for handling video rendering via Bun subprocess
 *
 * @requirements 4.1 - 使用内置 Bun 运行时执行渲染任务
 * @requirements 4.2 - 通过 Bun 子进程调用 Remotion 渲染 API
 */
export class RenderWorker implements IRenderWorker {
  /** Bun 路径解析器 */
  private bunResolver: BunPathResolver;

  /** 渲染脚本路径 */
  private renderScriptPath: string;

  /** 当前渲染子进程 */
  private currentProcess: ChildProcess | null = null;

  /** 是否已取消 */
  private isCancelled: boolean = false;

  /**
   * 创建 RenderWorker 实例
   *
   * @param bunResolver - Bun 路径解析器（可选，默认使用 createBunPathResolver()）
   * @param renderScriptPath - 渲染脚本路径（可选，默认自动检测）
   */
  constructor(bunResolver?: BunPathResolver, renderScriptPath?: string) {
    this.bunResolver = bunResolver ?? createBunPathResolver();
    this.renderScriptPath = renderScriptPath ?? this.detectRenderScriptPath();
    renderLog.info('RenderWorker initialized with Bun subprocess mode');
    renderLog.info(`Render script path: ${this.renderScriptPath}`);
  }

  /**
   * 检测渲染脚本路径
   */
  private detectRenderScriptPath(): string {
    // 在打包模式下，渲染脚本位于 app.asar 或 resources 目录
    if (app.isPackaged) {
      // 尝试多个可能的位置
      const possiblePaths = [
        // 打包后的路径（在 resources 目录下）
        join(process.resourcesPath, 'app.asar', 'apps', 'mcp-video', 'src', 'render-script.ts'),
        join(process.resourcesPath, 'apps', 'mcp-video', 'src', 'render-script.ts'),
        // 解压后的路径
        join(app.getAppPath(), '..', 'apps', 'mcp-video', 'src', 'render-script.ts'),
      ];

      for (const scriptPath of possiblePaths) {
        if (existsSync(scriptPath)) {
          return scriptPath;
        }
      }

      // 如果都找不到，返回默认路径（后续会报错）
      return join(process.resourcesPath, 'apps', 'mcp-video', 'src', 'render-script.ts');
    }

    // 开发模式：使用相对于项目根目录的路径
    // 从 apps/electron/src/main/video/ 向上找到项目根目录
    const devPath = join(__dirname, '..', '..', '..', '..', 'mcp-video', 'src', 'render-script.ts');
    if (existsSync(devPath)) {
      return devPath;
    }

    // 备用：从 app.getAppPath() 计算
    return join(app.getAppPath(), 'apps', 'mcp-video', 'src', 'render-script.ts');
  }

  /**
   * Render a video composition using Bun subprocess
   *
   * @param options - Render options
   * @returns Path to the rendered output file
   * @throws RenderError if rendering fails
   *
   * @requirements 4.1, 4.2, 4.3, 4.4, 4.5
   */
  async render(options: RenderWorkerOptions): Promise<string> {
    const {
      projectPath,
      compositionId,
      outputPath,
      quality,
      outputFormat = 'mp4',
      onProgress,
    } = options;

    // 重置状态
    this.isCancelled = false;
    this.currentProcess = null;

    renderLog.info(`Starting render: composition=${compositionId}, quality=${quality}, format=${outputFormat}`);

    // 报告初始状态
    this.reportProgress(onProgress, {
      status: 'bundling',
      progress: 0,
    });

    // 验证项目路径
    if (!existsSync(projectPath)) {
      const error = this.createError(
        RenderErrorType.INVALID_PROJECT,
        `Project path not found: ${projectPath}`
      );
      this.reportProgress(onProgress, {
        status: 'failed',
        progress: 0,
        error: error.message,
      });
      throw error;
    }

    // 确保输出目录存在
    const outputDir = dirname(outputPath);
    if (!existsSync(outputDir)) {
      mkdirSync(outputDir, { recursive: true });
    }

    // 获取 Bun 路径
    let bunPath: string;
    try {
      bunPath = this.bunResolver.getBunPath();
      renderLog.info(`Using Bun at: ${bunPath}`);
    } catch (error) {
      const renderError = this.createError(
        RenderErrorType.BUN_NOT_FOUND,
        `Bun runtime not found: ${error instanceof Error ? error.message : String(error)}`
      );
      this.reportProgress(onProgress, {
        status: 'failed',
        progress: 0,
        error: renderError.message,
      });
      throw renderError;
    }

    // 验证渲染脚本存在
    if (!existsSync(this.renderScriptPath)) {
      const error = this.createError(
        RenderErrorType.SCRIPT_NOT_FOUND,
        `Render script not found: ${this.renderScriptPath}`
      );
      this.reportProgress(onProgress, {
        status: 'failed',
        progress: 0,
        error: error.message,
      });
      throw error;
    }

    // 映射质量预设到渲染脚本支持的格式
    const scriptQuality = mapQualityPreset(quality);

    // 构建命令行参数
    const args = [
      'run',
      this.renderScriptPath,
      '--project', projectPath,
      '--composition', compositionId,
      '--output', outputPath,
      '--quality', scriptQuality,
      '--format', outputFormat,
    ];

    renderLog.info(`Spawning render process: ${bunPath} ${args.join(' ')}`);

    return new Promise<string>((resolve, reject) => {
      // 启动子进程
      // @requirements 4.2 - 通过 Bun 子进程调用 Remotion 渲染 API
      const childProcess = spawn(bunPath, args, {
        stdio: ['ignore', 'pipe', 'pipe'],
        detached: false,
        // 设置工作目录为项目路径
        cwd: projectPath,
      });

      this.currentProcess = childProcess;

      let lastError: string | undefined;
      let stdoutBuffer = '';

      // 处理 stdout - 解析 JSON 进度信息
      // @requirements 4.3 - 支持进度回调，报告渲染状态
      childProcess.stdout?.on('data', (data: Buffer) => {
        stdoutBuffer += data.toString();

        // 按行解析 JSON
        const lines = stdoutBuffer.split('\n');
        // 保留最后一个不完整的行
        stdoutBuffer = lines.pop() || '';

        for (const line of lines) {
          const trimmedLine = line.trim();
          if (!trimmedLine) continue;

          try {
            const progress = JSON.parse(trimmedLine) as RenderProgress;
            this.reportProgress(onProgress, progress);

            // 记录最后的错误信息
            if (progress.error) {
              lastError = progress.error;
            }
          } catch {
            // 非 JSON 输出，记录到日志
            renderLog.debug(`Non-JSON stdout: ${trimmedLine}`);
          }
        }
      });

      // 处理 stderr - 记录错误日志
      childProcess.stderr?.on('data', (data: Buffer) => {
        const message = data.toString().trim();
        if (message) {
          renderLog.warn(`Render stderr: ${message}`);
        }
      });

      // 处理进程退出
      childProcess.on('close', (code: number | null, signal: string | null) => {
        this.currentProcess = null;

        // 处理剩余的 stdout 缓冲
        if (stdoutBuffer.trim()) {
          try {
            const progress = JSON.parse(stdoutBuffer.trim()) as RenderProgress;
            this.reportProgress(onProgress, progress);
            if (progress.error) {
              lastError = progress.error;
            }
          } catch {
            // 忽略
          }
        }

        // 检查是否被取消
        // @requirements 4.4 - 支持取消正在进行的渲染任务
        if (this.isCancelled || signal === 'SIGTERM') {
          const error = this.createError(
            RenderErrorType.CANCELLED,
            'Render cancelled by user'
          );
          this.reportProgress(onProgress, {
            status: 'failed',
            progress: 0,
            error: error.message,
          });
          reject(error);
          return;
        }

        // 检查退出码
        if (code === 0) {
          // 成功
          this.reportProgress(onProgress, {
            status: 'completed',
            progress: 100,
          });
          renderLog.info(`Render completed: ${outputPath}`);
          resolve(outputPath);
        } else {
          // 失败 - 根据退出码分类错误
          // @requirements 4.5 - 返回分类的错误信息
          const errorType = EXIT_CODE_TO_ERROR_TYPE[code ?? -1] ?? RenderErrorType.RENDER_FAILED;
          const error = this.createError(
            errorType,
            lastError || `Render process exited with code ${code}`,
            `Exit code: ${code}, Signal: ${signal}`
          );
          this.reportProgress(onProgress, {
            status: 'failed',
            progress: 0,
            error: error.message,
          });
          renderLog.error(`Render failed:`, error);
          reject(error);
        }
      });

      // 处理进程错误
      childProcess.on('error', (err: Error) => {
        this.currentProcess = null;
        const error = this.createError(
          RenderErrorType.RENDER_FAILED,
          `Failed to spawn render process: ${err.message}`,
          err.stack
        );
        this.reportProgress(onProgress, {
          status: 'failed',
          progress: 0,
          error: error.message,
        });
        renderLog.error(`Render process error:`, error);
        reject(error);
      });
    });
  }

  /**
   * Cancel the current render operation
   *
   * @requirements 4.4 - 支持取消正在进行的渲染任务
   */
  cancel(): void {
    renderLog.info('Cancelling render...');
    this.isCancelled = true;

    if (this.currentProcess) {
      // 发送 SIGTERM 信号
      const killed = this.currentProcess.kill('SIGTERM');
      renderLog.info(`SIGTERM sent to render process, killed=${killed}`);

      // 如果 SIGTERM 没有效果，5 秒后发送 SIGKILL
      const process = this.currentProcess;
      setTimeout(() => {
        if (process && !process.killed) {
          renderLog.warn('Render process did not respond to SIGTERM, sending SIGKILL');
          process.kill('SIGKILL');
        }
      }, 5000);
    }
  }

  // ============================================================================
  // Private Helper Methods
  // ============================================================================

  /**
   * Report progress to the callback
   */
  private reportProgress(
    callback: ((progress: RenderProgress) => void) | undefined,
    progress: RenderProgress
  ): void {
    if (callback) {
      callback(progress);
    }
    renderLog.debug(`Progress: ${progress.status} - ${progress.progress}%`);
  }

  /**
   * Create a render error
   */
  private createError(type: RenderErrorType, message: string, details?: string): RenderError {
    return {
      type,
      message,
      details,
    };
  }
}
