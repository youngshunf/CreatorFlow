/**
 * 渲染引擎服务
 *
 * 基于 Remotion 的视频渲染核心
 * 始终使用 packages/video/src/Root.tsx 作为统一入口，compositionId 固定为 "SceneComposer"
 * 支持多种输出格式（MP4、WebM、GIF）和质量预设（draft、standard、high）
 *
 * @requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6
 */

import { existsSync, mkdirSync, statSync } from "fs";
import { dirname, resolve } from "path";
import {
  type RenderProgress,
  type RenderResult,
  type QualityPreset,
  type OutputFormat,
  type QualityPresetConfig,
  type Scene,
  type Transition,
  QUALITY_PRESETS,
} from "../types";
import {
  createRenderFailedError,
  createCompositionNotFoundError,
  MCPError,
} from "../types/errors";

// ============================================================================
// 类型定义
// ============================================================================

/**
 * SceneComposer 的 inputProps
 */
export interface SceneComposerInputProps {
  scenes: Scene[];
  transitions: Transition[];
}

/**
 * 渲染配置
 */
export interface RenderConfig {
  /** 输出格式 */
  outputFormat: OutputFormat;
  /** 质量预设 */
  quality: QualityPreset;
  /** SceneComposer 的输入参数（scenes + transitions） */
  inputProps: SceneComposerInputProps;
  /** 输出文件路径 */
  outputPath: string;
}

/**
 * 渲染错误类型
 */
export enum RenderErrorType {
  BUNDLE_FAILED = "BUNDLE_FAILED",
  COMPOSITION_NOT_FOUND = "COMPOSITION_NOT_FOUND",
  RENDER_FAILED = "RENDER_FAILED",
  OUTPUT_WRITE_FAILED = "OUTPUT_WRITE_FAILED",
  CANCELLED = "CANCELLED",
  INVALID_PROJECT = "INVALID_PROJECT",
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
   * 始终使用 packages/video/src/Root.tsx 作为入口，compositionId 固定为 "SceneComposer"。
   * 通过 inputProps 传入 scenes + transitions 数据。
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
    onProgress?: ProgressCallback,
  ): Promise<RenderResult> {
    const { outputFormat, quality, inputProps, outputPath } = config;

    // 重置取消状态
    this.isCancelled = false;
    this.cancelHandle = null;

    // 记录开始时间
    const startTime = Date.now();

    const compositionId = "SceneComposer";

    console.error(
      `[RenderEngine] Starting render: composition=${compositionId}, quality=${quality}, format=${outputFormat}, scenes=${inputProps.scenes.length}`,
    );

    try {
      // 报告 bundling 状态
      this.reportProgress(onProgress, {
        status: "bundling",
        progress: 0,
      });

      // 确保输出目录存在
      const outputDir = dirname(outputPath);
      if (!existsSync(outputDir)) {
        mkdirSync(outputDir, { recursive: true });
      }

      // 获取质量预设配置
      const presetConfig = this.getQualityConfig(quality);
      console.error(
        `[RenderEngine] Using quality preset: crf=${presetConfig.crf}, scale=${presetConfig.scale}, fps=${presetConfig.fps}`,
      );

      // 检查是否已取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 构建 Remotion bundle（使用统一入口）
      const bundleLocation = await this.bundleProject(onProgress);

      // 检查是否已取消
      if (this.isCancelled) {
        throw this.createCancelledError();
      }

      // 报告 preparing 状态
      this.reportProgress(onProgress, {
        status: "preparing",
        progress: 35,
      });

      // 动态导入 Remotion renderer
      const { renderMedia, selectComposition, makeCancelSignal } =
        await import("@remotion/renderer");

      // 创建取消信号
      const { cancelSignal, cancel } = makeCancelSignal();
      this.cancelHandle = { cancelSignal, cancel };

      // 选择 SceneComposer 组合，传入 inputProps 以触发 calculateMetadata
      const composition = await selectComposition({
        serveUrl: bundleLocation,
        id: compositionId,
        inputProps: inputProps as unknown as Record<string, unknown>,
      });

      if (!composition) {
        throw createCompositionNotFoundError(compositionId);
      }

      console.error(
        `[RenderEngine] Selected composition: ${composition.id} (${composition.width}x${composition.height}, ${composition.fps}fps, ${composition.durationInFrames} frames)`,
      );

      // 报告 preparing 完成
      this.reportProgress(onProgress, {
        status: "preparing",
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
      await renderMedia({
        composition: {
          ...composition,
          width: scaledWidth,
          height: scaledHeight,
          fps: presetConfig.fps,
        },
        serveUrl: bundleLocation,
        codec,
        outputLocation: outputPath,
        inputProps: inputProps as unknown as Record<string, unknown>,
        crf: presetConfig.crf,
        onProgress: ({ progress: renderProgress }: { progress: number }) => {
          const progress = 40 + Math.round(renderProgress * 60);
          this.reportProgress(onProgress, {
            status: "rendering",
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
        status: "completed",
        progress: 100,
      });

      // 计算渲染耗时
      const duration = Date.now() - startTime;

      // 获取文件大小
      const fileSize = this.getFileSize(outputPath);

      console.error(
        `[RenderEngine] Render completed: ${outputPath} (${fileSize} bytes, ${duration}ms)`,
      );

      return {
        success: true,
        outputPath,
        duration,
        fileSize,
      };
    } catch (error) {
      // 处理取消
      if (
        this.isCancelled ||
        (error instanceof Error && error.message.includes("cancelled"))
      ) {
        const cancelError = this.createCancelledError();
        this.reportProgress(onProgress, {
          status: "failed",
          progress: 0,
          error: cancelError.message,
        });
        throw createRenderFailedError(cancelError.message);
      }

      // @requirements 5.6 - 返回错误信息和部分进度
      const renderError = this.handleRenderError(error);
      this.reportProgress(onProgress, {
        status: "failed",
        progress: 0,
        error: renderError.message,
      });

      console.error(`[RenderEngine] Render failed:`, renderError);
      throw createRenderFailedError(renderError.message);
    } finally {
      this.cancelHandle = null;
    }
  }

  /**
   * 取消当前渲染操作
   */
  cancel(): void {
    console.error("[RenderEngine] Cancelling render...");
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
   * 获取统一入口文件路径
   *
   * 始终返回 packages/video/src/Root.tsx 的绝对路径
   */
  private getEntryPoint(): string {
    // __dirname 是 packages/video/src/mcp-server/services/
    // Root.tsx 在 packages/video/src/Root.tsx
    return resolve(__dirname, "../../Root.tsx");
  }

  /**
   * 构建 Remotion bundle
   *
   * 使用 packages/video/src/Root.tsx 作为统一入口
   *
   * @param onProgress 进度回调
   * @returns bundle 位置
   */
  private async bundleProject(onProgress?: ProgressCallback): Promise<string> {
    const { bundle } = await import("@remotion/bundler");

    this.reportProgress(onProgress, {
      status: "bundling",
      progress: 10,
    });

    const entryPoint = this.getEntryPoint();
    console.error(`[RenderEngine] Using entry point: ${entryPoint}`);

    const bundleLocation = await bundle({
      entryPoint,
      onProgress: (bundleProgress: number) => {
        const progress = 10 + Math.round(bundleProgress * 20);
        this.reportProgress(onProgress, {
          status: "bundling",
          progress,
        });
      },
    });

    console.error(`[RenderEngine] Bundle created at: ${bundleLocation}`);
    return bundleLocation;
  }

  /**
   * 根据输出格式获取编解码器
   *
   * @param format 输出格式
   * @returns 编解码器
   *
   * @requirements 5.2 - 支持 mp4 (H.264), webm (VP8), gif
   */
  private getCodecForFormat(format: OutputFormat): "h264" | "vp8" | "gif" {
    switch (format) {
      case "mp4":
        return "h264";
      case "webm":
        return "vp8"; // VP8 比 VP9 有更广泛的支持
      case "gif":
        return "gif";
      default:
        return "h264";
    }
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
    progress: RenderProgress,
  ): void {
    if (callback) {
      callback(progress);
    }
    console.error(
      `[RenderEngine] Progress: ${progress.status} - ${progress.progress}%`,
    );
  }

  /**
   * 创建取消错误
   */
  private createCancelledError(): RenderError {
    return {
      type: RenderErrorType.CANCELLED,
      message: "渲染已被用户取消",
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
      if (
        message.includes("bundle") ||
        message.includes("webpack") ||
        message.includes("esbuild")
      ) {
        return {
          type: RenderErrorType.BUNDLE_FAILED,
          message: "视频项目打包失败",
          details: error.message,
          stack: error.stack,
        };
      }

      // 组合未找到
      if (message.includes("composition") && message.includes("not found")) {
        return {
          type: RenderErrorType.COMPOSITION_NOT_FOUND,
          message: "找不到指定的视频组件",
          details: error.message,
          stack: error.stack,
        };
      }

      // 输出写入错误
      if (
        message.includes("write") ||
        message.includes("permission") ||
        message.includes("eacces")
      ) {
        return {
          type: RenderErrorType.OUTPUT_WRITE_FAILED,
          message: "无法写入输出文件",
          details: error.message,
          stack: error.stack,
        };
      }

      // 通用渲染错误
      return {
        type: RenderErrorType.RENDER_FAILED,
        message: "视频渲染失败",
        details: error.message,
        stack: error.stack,
      };
    }

    // 未知错误类型
    return {
      type: RenderErrorType.RENDER_FAILED,
      message: "渲染过程中发生未知错误",
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
