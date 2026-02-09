#!/usr/bin/env bun
/**
 * 独立渲染脚本
 *
 * 通过 Bun 子进程调用，执行 Remotion 视频渲染
 * 通过 stdout 输出 JSON 格式的进度信息
 *
 * 命令行参数：
 *   --project <path>      项目路径（必需）
 *   --composition <id>    组合 ID（必需）
 *   --output <path>       输出文件路径（必需）
 *   --quality <preset>    质量预设：draft, preview, production, high（默认：production）
 *   --format <format>     输出格式：mp4, webm, gif（默认：mp4）
 *
 * 进度输出格式（JSON，每行一个）：
 *   { "status": "bundling" | "preparing" | "rendering" | "completed" | "failed", "progress": 0-100, "error"?: string }
 *
 * 退出码：
 *   0 - 成功
 *   1 - 参数错误
 *   2 - 项目路径不存在
 *   3 - 组合未找到
 *   4 - 打包失败
 *   5 - 渲染失败
 *   6 - 输出写入失败
 *
 * @requirements 4.1, 4.2, 4.3, 4.6
 */

import { existsSync, mkdirSync, statSync } from 'fs';
import { join, dirname } from 'path';
import { parseArgs } from 'util';

// ============================================================================
// 类型定义
// ============================================================================

/**
 * 渲染进度信息
 * @requirements 4.3 - 支持进度回调，报告渲染状态
 */
interface RenderProgress {
  status: 'bundling' | 'preparing' | 'rendering' | 'completed' | 'failed';
  progress: number; // 0-100
  error?: string;
}

/**
 * 质量预设类型
 * @requirements 4.6 - 支持多种质量预设
 */
type QualityPreset = 'draft' | 'preview' | 'production' | 'high';

/**
 * 输出格式类型
 * @requirements 4.6 - 支持多种输出格式
 */
type OutputFormat = 'mp4' | 'webm' | 'gif';

/**
 * 质量预设配置
 */
interface QualityPresetConfig {
  /** Constant Rate Factor (lower = better quality, larger file) */
  crf: number;
  /** Scale factor (1 = original size) */
  scale: number;
  /** Frames per second */
  fps: number;
}

/**
 * 命令行参数
 */
interface RenderArgs {
  projectPath: string;
  compositionId: string;
  outputPath: string;
  quality: QualityPreset;
  format: OutputFormat;
}

// ============================================================================
// 常量定义
// ============================================================================

/**
 * 质量预设配置映射
 * @requirements 4.6 - 支持 draft, preview, production, high 质量预设
 */
const QUALITY_PRESETS: Record<QualityPreset, QualityPresetConfig> = {
  draft: { crf: 32, scale: 0.5, fps: 15 },
  preview: { crf: 28, scale: 0.75, fps: 24 },
  production: { crf: 18, scale: 1, fps: 30 },
  high: { crf: 12, scale: 1, fps: 60 },
};

/**
 * 退出码定义
 */
const EXIT_CODES = {
  SUCCESS: 0,
  INVALID_ARGS: 1,
  PROJECT_NOT_FOUND: 2,
  COMPOSITION_NOT_FOUND: 3,
  BUNDLE_FAILED: 4,
  RENDER_FAILED: 5,
  OUTPUT_WRITE_FAILED: 6,
} as const;

// ============================================================================
// 进度输出函数
// ============================================================================

/**
 * 输出进度信息到 stdout（JSON 格式）
 * @requirements 4.3 - 通过 stdout 输出 JSON 格式的进度信息
 */
function reportProgress(progress: RenderProgress): void {
  console.log(JSON.stringify(progress));
}

/**
 * 输出错误信息到 stderr
 */
function reportError(message: string): void {
  console.error(`[render-script] Error: ${message}`);
}

// ============================================================================
// 参数解析
// ============================================================================

/**
 * 解析命令行参数
 */
function parseArguments(): RenderArgs | null {
  try {
    const { values } = parseArgs({
      options: {
        project: { type: 'string', short: 'p' },
        composition: { type: 'string', short: 'c' },
        output: { type: 'string', short: 'o' },
        quality: { type: 'string', short: 'q' },
        format: { type: 'string', short: 'f' },
      },
      strict: true,
    });

    // 验证必需参数
    if (!values.project) {
      reportError('Missing required argument: --project');
      return null;
    }
    if (!values.composition) {
      reportError('Missing required argument: --composition');
      return null;
    }
    if (!values.output) {
      reportError('Missing required argument: --output');
      return null;
    }

    // 验证质量预设
    const quality = (values.quality || 'production') as QualityPreset;
    if (!['draft', 'preview', 'production', 'high'].includes(quality)) {
      reportError(`Invalid quality preset: ${quality}. Must be one of: draft, preview, production, high`);
      return null;
    }

    // 验证输出格式
    const format = (values.format || 'mp4') as OutputFormat;
    if (!['mp4', 'webm', 'gif'].includes(format)) {
      reportError(`Invalid output format: ${format}. Must be one of: mp4, webm, gif`);
      return null;
    }

    return {
      projectPath: values.project,
      compositionId: values.composition,
      outputPath: values.output,
      quality,
      format,
    };
  } catch (error) {
    reportError(`Failed to parse arguments: ${error instanceof Error ? error.message : String(error)}`);
    return null;
  }
}

// ============================================================================
// 渲染逻辑
// ============================================================================

/**
 * 查找 Remotion 入口文件
 */
function findEntryPoint(projectPath: string): string {
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

  // 尝试使用 @creator-flow/video 包的 Root
  try {
    const videoPackageRoot = require.resolve('@creator-flow/video');
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
 * @requirements 4.6 - 支持 mp4, webm, gif 格式
 */
function getCodecForFormat(format: OutputFormat): 'h264' | 'vp8' | 'gif' {
  switch (format) {
    case 'mp4':
      return 'h264';
    case 'webm':
      return 'vp8';
    case 'gif':
      return 'gif';
    default:
      return 'h264';
  }
}

/**
 * 获取文件大小
 */
function getFileSize(filePath: string): number {
  try {
    const stats = statSync(filePath);
    return stats.size;
  } catch {
    return 0;
  }
}

/**
 * 执行渲染
 * @requirements 4.1, 4.2 - 使用 Bun 运行时调用 Remotion 渲染 API
 */
async function render(args: RenderArgs): Promise<number> {
  const { projectPath, compositionId, outputPath, quality, format } = args;
  const startTime = Date.now();

  try {
    // 报告 bundling 状态
    reportProgress({ status: 'bundling', progress: 0 });

    // 验证项目路径
    if (!existsSync(projectPath)) {
      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Project path does not exist: ${projectPath}`,
      });
      return EXIT_CODES.PROJECT_NOT_FOUND;
    }

    // 确保输出目录存在
    const outputDir = dirname(outputPath);
    if (!existsSync(outputDir)) {
      try {
        mkdirSync(outputDir, { recursive: true });
      } catch (error) {
        reportProgress({
          status: 'failed',
          progress: 0,
          error: `Failed to create output directory: ${outputDir}`,
        });
        return EXIT_CODES.OUTPUT_WRITE_FAILED;
      }
    }

    // 获取质量预设配置
    const presetConfig = QUALITY_PRESETS[quality];

    // 查找入口文件
    const entryPoint = findEntryPoint(projectPath);

    // 动态导入 Remotion bundler
    reportProgress({ status: 'bundling', progress: 5 });
    const { bundle } = await import('@remotion/bundler');

    // 构建 bundle
    let bundleLocation: string;
    try {
      bundleLocation = await bundle({
        entryPoint,
        onProgress: (bundleProgress: number) => {
          // 将 bundle 进度 (0-1) 映射到我们的进度范围 (5-30)
          const progress = 5 + Math.round(bundleProgress * 25);
          reportProgress({ status: 'bundling', progress });
        },
      });
    } catch (error) {
      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Bundle failed: ${error instanceof Error ? error.message : String(error)}`,
      });
      return EXIT_CODES.BUNDLE_FAILED;
    }

    // 报告 preparing 状态
    reportProgress({ status: 'preparing', progress: 35 });

    // 动态导入 Remotion renderer
    const { renderMedia, selectComposition } = await import('@remotion/renderer');

    // 选择组合
    let composition;
    try {
      composition = await selectComposition({
        serveUrl: bundleLocation,
        id: compositionId,
      });
    } catch (error) {
      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Composition not found: ${compositionId}`,
      });
      return EXIT_CODES.COMPOSITION_NOT_FOUND;
    }

    if (!composition) {
      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Composition not found: ${compositionId}`,
      });
      return EXIT_CODES.COMPOSITION_NOT_FOUND;
    }

    // 报告 preparing 完成
    reportProgress({ status: 'preparing', progress: 40 });

    // 获取编解码器
    const codec = getCodecForFormat(format);

    // 计算缩放后的尺寸
    const scaledWidth = Math.round(composition.width * presetConfig.scale);
    const scaledHeight = Math.round(composition.height * presetConfig.scale);

    // 渲染视频
    try {
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
        crf: presetConfig.crf,
        onProgress: ({ progress: renderProgress }: { progress: number }) => {
          // 将渲染进度 (0-1) 映射到我们的进度范围 (40-99)
          const progress = 40 + Math.round(renderProgress * 59);
          reportProgress({ status: 'rendering', progress });
        },
      });
    } catch (error) {
      // 检查是否是输出写入错误
      const errorMessage = error instanceof Error ? error.message : String(error);
      if (errorMessage.includes('write') || errorMessage.includes('permission') || errorMessage.includes('EACCES')) {
        reportProgress({
          status: 'failed',
          progress: 0,
          error: `Output write failed: ${errorMessage}`,
        });
        return EXIT_CODES.OUTPUT_WRITE_FAILED;
      }

      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Render failed: ${errorMessage}`,
      });
      return EXIT_CODES.RENDER_FAILED;
    }

    // 验证输出文件
    if (!existsSync(outputPath)) {
      reportProgress({
        status: 'failed',
        progress: 0,
        error: `Output file was not created: ${outputPath}`,
      });
      return EXIT_CODES.OUTPUT_WRITE_FAILED;
    }

    // 获取文件大小和渲染时间
    const fileSize = getFileSize(outputPath);
    const duration = Date.now() - startTime;

    // 报告完成
    reportProgress({ status: 'completed', progress: 100 });

    // 输出渲染统计信息到 stderr（不影响 JSON 进度输出）
    console.error(`[render-script] Render completed: ${outputPath}`);
    console.error(`[render-script] File size: ${fileSize} bytes`);
    console.error(`[render-script] Duration: ${duration}ms`);

    return EXIT_CODES.SUCCESS;
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : String(error);
    reportProgress({
      status: 'failed',
      progress: 0,
      error: `Unexpected error: ${errorMessage}`,
    });
    return EXIT_CODES.RENDER_FAILED;
  }
}

// ============================================================================
// 主入口
// ============================================================================

async function main(): Promise<void> {
  // 解析命令行参数
  const args = parseArguments();
  if (!args) {
    process.exit(EXIT_CODES.INVALID_ARGS);
  }

  // 执行渲染
  const exitCode = await render(args);
  process.exit(exitCode);
}

// 运行主函数
main().catch((error) => {
  reportError(`Unhandled error: ${error instanceof Error ? error.message : String(error)}`);
  process.exit(EXIT_CODES.RENDER_FAILED);
});
