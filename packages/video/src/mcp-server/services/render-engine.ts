/**
 * 渲染引擎服务
 *
 * 基于 Remotion 的视频渲染核心
 * 支持多种输出格式（MP4、WebM、GIF）和质量预设（draft、standard、high）
 *
 * @requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6
 */

import { existsSync, mkdirSync, statSync } from 'fs';
import { join, dirname } from 'path';
import { randomUUID } from 'crypto';
import {
  type RenderProgress,
  type RenderResult,
  type QualityPreset,
  type OutputFormat,
  type QualityPresetConfig,
  QUALITY_PRESETS,
} from '../types';
import {
  createRenderFailedError,
  createFileNotFoundError,
  createCompositionNotFoundError,
  MCPError,
} from '../types/errors';
import { getOutputPath, getOutputFilePath } from '../utils/paths';

// ============================================================================
// 类型定义
// ============================================================================

/**
 * 渲染配置
 */
export interface RenderConfig {
  /** 项目路径（包含 Remotion 入口文件的目录） */
  projectPath: string;
  /** 组合 ID */
  compositionId: string;
  /** 输出格式 */
  outputFormat: OutputFormat;
  /** 质量预设 */
  quality: QualityPreset;
  /** 输出文件路径（可选，默认自动生成） */
  outputPath?: string;
}

/**
 * 渲染错误类型
 */
export enum RenderErrorType {
  BUNDLE_FAILED = 'BUNDLE_FAILED',
  COMPOSITION_NOT_FOUND = 'COMPOSITION_NOT_FOUND',
  RENDER_FAILED = 'RENDER_FAILED',
  OUTPUT_WRITE_FAILED = 'OUTPUT_WRITE_FAILED',
  CANCELLED = 'CANCELLED',
  INVALID_PROJECT = 'INVALID_PROJECT',
}

/**
 * 渲染错误结构
 */
export interface RenderError {
  type: RenderErrorType;
  message: string;
  details?: string;
  stack?: string;
}

/**
 * 进度回调函数类型
 */
export type ProgressCallback = (progress: RenderProgress) => void;

/**
 * Remotion 取消信号类型
 */
interface CancelSignalHandle {
  cancelSignal: (callback: () => void) => void;
  cancel: () => void;
}

// ============================================================================
// RenderEngine 类
// ============================================================================

/**
 * 渲染引擎服务
 *
 * 提供视频渲染功能，支持：
 * - Remotion bundle 构建
 * - renderMedia 调用
 * - 质量预设配置
 * - 进度回调
 */
export class RenderEngine {
  /** 是否已取消渲染 */
  private isCancelled: boolean = false;

  /** 当前渲染的取消句柄 */
  private cancelHandle: CancelSignalHandle | null = null;

  constructor() {
    // 初始化
  }

  // ==========================================================================
  // 公共方法
  // ==========================================================================

  /**
   * 渲染视频
   *
   * @param config 渲染配置
   * @param onProgress 进度回调函数（可选）
   * @returns 渲染结果
   * @throws MCPError 如果渲染失败
   *
   * @requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6
   */
  async render(
    config: RenderConfig,
    onProgress?: ProgressCallback
  ): Promise<RenderResult> {
    const { projectPath, compositionId, outputFormat, quality, outputPath } = config;

    // 重置取消状态
    this.isCancelled = false;
    this.cancelHandle = null;

    // 记录开始时间
    const startTime = Date.now();

    console.error(`[RenderEngine] Starting render: composition=${compositionId}, quality=${quality}, format=${outputFormat}`);

    try {
      // 报告 bundling 状态
      this.reportProgress(onProgress, {
        status: 'bundling',
        progress: 0,
      });

      // 验证项目路径
      if (!existsSync(projectPath)) {
        throw createFileNotFoundError(projectPath);
      }

      // 确定输出文件路径
      const finalOutputPath = outputPath ?? this.generateOutputPath(projectPath, compositionId, outputFormat);

      // 确保输出目录存在
      const outputDir = dirname(finalOutputPath);
      if (!existsSync(outputDir)) {
        mkdirSync(outputDir, { recursive: true });
      }

      // 获取质量预设配置
      const presetConfig = this.getQualityConfig(quality);
      console.error(`[RenderEngine] Using quality preset: crf=${presetConfig.crf}, scale=${presetConfig.scale}, fps=${presetConfig.fps}`);

      // 检查是否已取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 构建 Remotion bundle
      const bundleLocation = await this.bundleProject(projectPath, onProgress);

      // 检查是否已取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 报告 preparing 状态
      this.reportProgress(onProgress, {
        status: 'preparing',
        progress: 35,
      });

      // 动态导入 Remotion renderer
      const { renderMedia, selectComposition, makeCancelSignal } = await import('@remotion/renderer');

      // 创建取消信号
      const { cancelSignal, cancel } = makeCancelSignal();
      this.cancelHandle = { cancelSignal, cancel };

      // 选择组合
      const composition = await selectComposition({
        serveUrl: bundleLocation,
        id: compositionId,
      });

      if (!composition) {
        throw createCompositionNotFoundError(compositionId);
      }

      console.error(`[RenderEngine] Selected composition: ${composition.id} (${composition.width}x${composition.height}, ${composition.fps}fps, ${composition.durationInFrames} frames)`);

      // 报告 preparing 完成
      this.reportProgress(onProgress, {
        status: 'preparing',
        progress: 40,
      });

      // 检查是否已取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 获取编解码器
      const codec = this.getCodecForFormat(outputFormat);

      // 计算缩放后的尺寸
      const scaledWidth = Math.round(composition.width * presetConfig.scale);
      const scaledHeight = Math.round(composition.height * presetConfig.scale);

      // 渲染视频
      // @requirements 5.1 - 启动渲染指定的组合
      // @requirements 5.2 - 支持 MP4、WebM、GIF 格式
      // @requirements 5.3 - 支持质量预设
      // @requirements 5.4 - 提供进度更新
      await renderMedia({
        composition: {
          ...composition,
          width: scaledWidth,
          height: scaledHeight,
          fps: presetConfig.fps,
        },
        serveUrl: bundleLocation,
        codec,
        outputLocation: finalOutputPath,
        crf: presetConfig.crf,
        onProgress: ({ progress: renderProgress }: { progress: number }) => {
          // 将渲染进度 (0-1) 映射到我们的进度范围 (40-100)
          const progress = 40 + Math.round(renderProgress * 60);
          this.reportProgress(onProgress, {
            status: 'rendering',
            progress,
          });
        },
        cancelSignal,
      });

      // 检查是否在渲染过程中被取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 报告完成
      this.reportProgress(onProgress, {
        status: 'completed',
        progress: 100,
      });

      // 计算渲染耗时
      const duration = Date.now() - startTime;

      // 获取文件大小
      const fileSize = this.getFileSize(finalOutputPath);

      console.error(`[RenderEngine] Render completed: ${finalOutputPath} (${fileSize} bytes, ${duration}ms)`);

      // @requirements 5.5 - 返回输出文件路径和渲染统计
      return {
        success: true,
        outputPath: finalOutputPath,
        duration,
        fileSize,
      };
    } catch (error) {
      // 处理取消
      if (this.isCancelled || (error instanceof Error && error.message.includes('cancelled'))) {
        const cancelError = this.createCancelledError();
        this.reportProgress(onProgress, {
          status: 'failed',
          progress: 0,
          error: cancelError.message,
        });
        throw createRenderFailedError(cancelError.message);
      }

      // @requirements 5.6 - 返回错误信息和部分进度
      const renderError = this.handleRenderError(error);
      this.reportProgress(onProgress, {
        status: 'failed',
        progress: 0,
        error: renderError.message,
      });

      console.error(`[RenderEngine] Render failed:`, renderError);
      throw createRenderFailedError(renderError.message, {
        path: config.projectPath,
      });
    } finally {
      this.cancelHandle = null;
    }
  }

  /**
   * 取消当前渲染操作
   */
  cancel(): void {
    console.error('[RenderEngine] Cancelling render...');
    this.isCancelled = true;
    if (this.cancelHandle) {
      this.cancelHandle.cancel();
    }
  }

  // ==========================================================================
  // 私有方法
  // ==========================================================================

  /**
   * 获取质量预设配置
   *
   * @param quality 质量预设
   * @returns 质量预设配置
   *
   * @requirements 5.3 - 支持 draft、standard、high 质量预设
   */
  private getQualityConfig(quality: QualityPreset): QualityPresetConfig {
    return QUALITY_PRESETS[quality];
  }

  /**
   * 构建 Remotion bundle
   *
   * @param projectPath 项目路径
   * @param onProgress 进度回调
   * @returns bundle 位置
   *
   * @requirements 5.1 - 使用 @remotion/bundler 构建项目
   */
  private async bundleProject(
    projectPath: string,
    onProgress?: ProgressCallback
  ): Promise<string> {
    // 动态导入 Remotion bundler
    const { bundle } = await import('@remotion/bundler');

    // 报告 bundling 进度
    this.reportProgress(onProgress, {
      status: 'bundling',
      progress: 10,
    });

    // 查找入口文件
    const entryPoint = this.findEntryPoint(projectPath);
    console.error(`[RenderEngine] Using entry point: ${entryPoint}`);

    // 构建 bundle
    const bundleLocation = await bundle({
      entryPoint,
      onProgress: (bundleProgress: number) => {
        // 将 bundle 进度 (0-1) 映射到我们的进度范围 (10-30)
        const progress = 10 + Math.round(bundleProgress * 20);
        this.reportProgress(onProgress, {
          status: 'bundling',
          progress,
        });
      },
    });

    console.error(`[RenderEngine] Bundle created at: ${bundleLocation}`);
    return bundleLocation;
  }

  /**
   * 查找 Remotion 入口文件
   *
   * @param projectPath 项目路径
   * @returns 入口文件路径
   */
  private findEntryPoint(projectPath: string): string {
    // 尝试常见的入口文件位置
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

    // 如果项目中没有找到入口文件，尝试使用 @sprouty-ai/video 包的 Root
    // 这是使用模板的项目的默认情况
    try {
      const videoPackageRoot = require.resolve('@sprouty-ai/video');
      const videoPackageDir = dirname(videoPackageRoot);
      const defaultEntryPoint = join(videoPackageDir, 'Root.tsx');
      if (existsSync(defaultEntryPoint)) {
        return defaultEntryPoint;
      }
    } catch {
      // 忽略 video 包未找到的情况
    }

    // 回退到项目路径本身
    return projectPath;
  }

  /**
   * 根据输出格式获取编解码器
   *
   * @param format 输出格式
   * @returns 编解码器
   *
   * @requirements 5.2 - 支持 mp4 (H.264), webm (VP8), gif
   */
  private getCodecForFormat(format: OutputFormat): 'h264' | 'vp8' | 'gif' {
    switch (format) {
      case 'mp4':
        return 'h264';
      case 'webm':
        return 'vp8'; // VP8 比 VP9 有更广泛的支持
      case 'gif':
        return 'gif';
      default:
        return 'h264';
    }
  }

  /**
   * 生成输出文件路径
   *
   * @param projectPath 项目路径
   * @param compositionId 组合 ID
   * @param format 输出格式
   * @returns 输出文件路径
   */
  private generateOutputPath(
    projectPath: string,
    compositionId: string,
    format: OutputFormat
  ): string {
    // 生成时间戳
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);

    // 生成文件名
    const fileName = `${compositionId}_${timestamp}.${format}`;

    // 查找输出目录
    // 首先尝试项目的 "输出" 目录
    const outputDir = join(projectPath, '输出');
    if (existsSync(outputDir) || existsSync(dirname(outputDir))) {
      return join(outputDir, fileName);
    }

    // 回退到项目目录下的 output 目录
    return join(projectPath, 'output', fileName);
  }

  /**
   * 获取文件大小
   *
   * @param filePath 文件路径
   * @returns 文件大小（字节）
   */
  private getFileSize(filePath: string): number {
    try {
      const stats = statSync(filePath);
      return stats.size;
    } catch {
      return 0;
    }
  }

  /**
   * 报告进度
   *
   * @param callback 进度回调
   * @param progress 进度信息
   */
  private reportProgress(
    callback: ProgressCallback | undefined,
    progress: RenderProgress
  ): void {
    if (callback) {
      callback(progress);
    }
    console.error(`[RenderEngine] Progress: ${progress.status} - ${progress.progress}%`);
  }

  /**
   * 创建取消错误
   */
  private createCancelledError(): RenderError {
    return {
      type: RenderErrorType.CANCELLED,
      message: '渲染已被用户取消',
    };
  }

  /**
   * 处理渲染错误
   *
   * @param error 原始错误
   * @returns 渲染错误
   *
   * @requirements 5.6 - 返回详细的错误信息
   */
  private handleRenderError(error: unknown): RenderError {
    if (error instanceof MCPError) {
      return {
        type: RenderErrorType.RENDER_FAILED,
        message: error.message,
        details: JSON.stringify(error.details),
      };
    }

    if (error instanceof Error) {
      const message = error.message.toLowerCase();

      // Bundle 错误
      if (message.includes('bundle') || message.includes('webpack') || message.includes('esbuild')) {
        return {
          type: RenderErrorType.BUNDLE_FAILED,
          message: '视频项目打包失败',
          details: error.message,
          stack: error.stack,
        };
      }

      // 组合未找到
      if (message.includes('composition') && message.includes('not found')) {
        return {
          type: RenderErrorType.COMPOSITION_NOT_FOUND,
          message: '找不到指定的视频组件',
          details: error.message,
          stack: error.stack,
        };
      }

      // 输出写入错误
      if (message.includes('write') || message.includes('permission') || message.includes('eacces')) {
        return {
          type: RenderErrorType.OUTPUT_WRITE_FAILED,
          message: '无法写入输出文件',
          details: error.message,
          stack: error.stack,
        };
      }

      // 通用渲染错误
      return {
        type: RenderErrorType.RENDER_FAILED,
        message: '视频渲染失败',
        details: error.message,
        stack: error.stack,
      };
    }

    // 未知错误类型
    return {
      type: RenderErrorType.RENDER_FAILED,
      message: '渲染过程中发生未知错误',
      details: String(error),
    };
  }

  // ==========================================================================
  // 静态工厂方法
  // ==========================================================================

  /**
   * 创建 RenderEngine 实例
   * @returns RenderEngine 实例
   */
  static create(): RenderEngine {
    return new RenderEngine();
  }
}

// ============================================================================
// 导出单例实例
// ============================================================================

/**
 * 默认的渲染引擎实例
 */
export const renderEngine = RenderEngine.create();
