#!/usr/bin/env bun
/**
 * 独立渲染脚本（v2 — 内置架构）
 *
 * 通过 Bun 子进程调用，执行 Remotion 视频渲染。
 * 数据来源：Electron main 进程通过 --input-props 传入 JSON（scenes + transitions）。
 * 不再依赖 project.json 文件。
 *
 * 命令行参数：
 *   --composition <id>      组合 ID（必需，固定为 SceneComposer）
 *   --output <path>         输出文件路径（必需）
 *   --input-props <json>    SceneComposer 的 inputProps JSON（必需）
 *   --quality <preset>      质量预设：draft, preview, production, high（默认：production）
 *   --format <format>       输出格式：mp4, webm, gif（默认：mp4）
 *
 * 进度输出格式（JSON，每行一个 stdout）：
 *   { "status": "bundling" | "preparing" | "rendering" | "completed" | "failed", "progress": 0-100, "error"?: string }
 *
 * 退出码：
 *   0 - 成功
 *   1 - 参数错误
 *   3 - 组合未找到
 *   4 - 打包失败
 *   5 - 渲染失败
 *   6 - 输出写入失败
 */

import { existsSync, mkdirSync, statSync } from 'fs';
import { join, dirname } from 'path';
import { parseArgs } from 'util';

// ============================================================================
// Types
// ============================================================================

interface RenderProgress {
  status: 'bundling' | 'preparing' | 'rendering' | 'completed' | 'failed';
  progress: number;
  error?: string;
}

type QualityPreset = 'draft' | 'preview' | 'production' | 'high';
type OutputFormat = 'mp4' | 'webm' | 'gif';

interface QualityPresetConfig {
  crf: number;
  scale: number;
  fps: number;
}

interface RenderArgs {
  compositionId: string;
  outputPath: string;
  inputProps: Record<string, unknown>;
  quality: QualityPreset;
  format: OutputFormat;
}

// ============================================================================
// Constants
// ============================================================================

const QUALITY_PRESETS: Record<QualityPreset, QualityPresetConfig> = {
  draft: { crf: 32, scale: 0.5, fps: 15 },
  preview: { crf: 28, scale: 0.75, fps: 24 },
  production: { crf: 18, scale: 1, fps: 30 },
  high: { crf: 12, scale: 1, fps: 60 },
};

const EXIT_CODES = {
  SUCCESS: 0,
  INVALID_ARGS: 1,
  COMPOSITION_NOT_FOUND: 3,
  BUNDLE_FAILED: 4,
  RENDER_FAILED: 5,
  OUTPUT_WRITE_FAILED: 6,
} as const;

// ============================================================================
// Progress output
// ============================================================================

function reportProgress(progress: RenderProgress): void {
  console.log(JSON.stringify(progress));
}

function reportError(message: string): void {
  console.error(`[render-script] Error: ${message}`);
}

// ============================================================================
// Argument parsing
// ============================================================================

function parseArguments(): RenderArgs | null {
  try {
    const { values } = parseArgs({
      options: {
        composition: { type: 'string', short: 'c' },
        output: { type: 'string', short: 'o' },
        quality: { type: 'string', short: 'q' },
        format: { type: 'string', short: 'f' },
        'input-props': { type: 'string' },
      },
      strict: true,
    });

    if (!values.composition) {
      reportError('Missing required argument: --composition');
      return null;
    }
    if (!values.output) {
      reportError('Missing required argument: --output');
      return null;
    }
    if (!values['input-props']) {
      reportError('Missing required argument: --input-props');
      return null;
    }

    // Parse input props JSON
    let inputProps: Record<string, unknown>;
    try {
      inputProps = JSON.parse(values['input-props']);
    } catch {
      reportError('Invalid JSON for --input-props');
      return null;
    }

    const quality = (values.quality || 'production') as QualityPreset;
    if (!['draft', 'preview', 'production', 'high'].includes(quality)) {
      reportError(`Invalid quality preset: ${quality}. Must be one of: draft, preview, production, high`);
      return null;
    }

    const format = (values.format || 'mp4') as OutputFormat;
    if (!['mp4', 'webm', 'gif'].includes(format)) {
      reportError(`Invalid output format: ${format}. Must be one of: mp4, webm, gif`);
      return null;
    }

    return {
      compositionId: values.composition,
      outputPath: values.output,
      inputProps,
      quality,
      format,
    };
  } catch (error) {
    reportError(`Failed to parse arguments: ${error instanceof Error ? error.message : String(error)}`);
    return null;
  }
}

// ============================================================================
// Render logic
// ============================================================================

function getCodecForFormat(format: OutputFormat): 'h264' | 'vp8' | 'gif' {
  switch (format) {
    case 'mp4': return 'h264';
    case 'webm': return 'vp8';
    case 'gif': return 'gif';
    default: return 'h264';
  }
}

/**
 * 查找 Remotion 入口文件（packages/video/src/Root.tsx）
 */
function findEntryPoint(): string {
  // render-script 位于 packages/video/src/render/render-script.ts
  // Root.tsx 位于 packages/video/src/Root.tsx
  const rootFromScript = join(__dirname, '..', 'Root.tsx');
  if (existsSync(rootFromScript)) {
    return rootFromScript;
  }

  // 备用：通过 resolve
  try {
    const videoPackageRoot = require.resolve('@sprouty-ai/video');
    const videoPackageDir = dirname(videoPackageRoot);
    const defaultEntryPoint = join(videoPackageDir, 'Root.tsx');
    if (existsSync(defaultEntryPoint)) {
      return defaultEntryPoint;
    }
  } catch {
    // ignore
  }

  throw new Error('Cannot find Remotion entry point (Root.tsx)');
}

async function render(args: RenderArgs): Promise<number> {
  const { compositionId, outputPath, inputProps, quality, format } = args;
  const startTime = Date.now();

  try {
    reportProgress({ status: 'bundling', progress: 0 });

    // Ensure output directory exists
    const outputDir = dirname(outputPath);
    if (!existsSync(outputDir)) {
      try {
        mkdirSync(outputDir, { recursive: true });
      } catch {
        reportProgress({ status: 'failed', progress: 0, error: `Failed to create output directory: ${outputDir}` });
        return EXIT_CODES.OUTPUT_WRITE_FAILED;
      }
    }

    const presetConfig = QUALITY_PRESETS[quality];
    const entryPoint = findEntryPoint();

    // Bundle
    reportProgress({ status: 'bundling', progress: 5 });
    const { bundle } = await import('@remotion/bundler');

    let bundleLocation: string;
    try {
      bundleLocation = await bundle({
        entryPoint,
        onProgress: (bundleProgress: number) => {
          const progress = 5 + Math.round(bundleProgress * 25);
          reportProgress({ status: 'bundling', progress });
        },
      });
    } catch (error) {
      reportProgress({ status: 'failed', progress: 0, error: `Bundle failed: ${error instanceof Error ? error.message : String(error)}` });
      return EXIT_CODES.BUNDLE_FAILED;
    }

    // Select composition
    reportProgress({ status: 'preparing', progress: 35 });
    const { renderMedia, selectComposition } = await import('@remotion/renderer');

    let composition;
    try {
      composition = await selectComposition({
        serveUrl: bundleLocation,
        id: compositionId,
        inputProps,
      });
    } catch (error) {
      reportProgress({ status: 'failed', progress: 0, error: `Composition not found: ${compositionId}` });
      return EXIT_CODES.COMPOSITION_NOT_FOUND;
    }

    if (!composition) {
      reportProgress({ status: 'failed', progress: 0, error: `Composition not found: ${compositionId}` });
      return EXIT_CODES.COMPOSITION_NOT_FOUND;
    }

    reportProgress({ status: 'preparing', progress: 40 });

    const codec = getCodecForFormat(format);
    const scaledWidth = Math.round(composition.width * presetConfig.scale);
    const scaledHeight = Math.round(composition.height * presetConfig.scale);

    // Render
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
        inputProps,
        onProgress: ({ progress: renderProgress }: { progress: number }) => {
          const progress = 40 + Math.round(renderProgress * 59);
          reportProgress({ status: 'rendering', progress });
        },
      });
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : String(error);
      if (errorMessage.includes('write') || errorMessage.includes('permission') || errorMessage.includes('EACCES')) {
        reportProgress({ status: 'failed', progress: 0, error: `Output write failed: ${errorMessage}` });
        return EXIT_CODES.OUTPUT_WRITE_FAILED;
      }
      reportProgress({ status: 'failed', progress: 0, error: `Render failed: ${errorMessage}` });
      return EXIT_CODES.RENDER_FAILED;
    }

    // Verify output
    if (!existsSync(outputPath)) {
      reportProgress({ status: 'failed', progress: 0, error: `Output file was not created: ${outputPath}` });
      return EXIT_CODES.OUTPUT_WRITE_FAILED;
    }

    const fileSize = statSync(outputPath).size;
    const duration = Date.now() - startTime;

    reportProgress({ status: 'completed', progress: 100 });
    console.error(`[render-script] Render completed: ${outputPath} (${fileSize} bytes, ${duration}ms)`);

    return EXIT_CODES.SUCCESS;
  } catch (error) {
    reportProgress({ status: 'failed', progress: 0, error: `Unexpected error: ${error instanceof Error ? error.message : String(error)}` });
    return EXIT_CODES.RENDER_FAILED;
  }
}

// ============================================================================
// Main
// ============================================================================

async function main(): Promise<void> {
  const args = parseArguments();
  if (!args) {
    process.exit(EXIT_CODES.INVALID_ARGS);
  }
  const exitCode = await render(args);
  process.exit(exitCode);
}

main().catch((error) => {
  reportError(`Unhandled error: ${error instanceof Error ? error.message : String(error)}`);
  process.exit(EXIT_CODES.RENDER_FAILED);
});
