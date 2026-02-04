/**
 * Type definitions for @creator-flow/video
 *
 * This file contains all Zod schemas and TypeScript types
 * for video projects, compositions, assets, and rendering.
 */
import { z } from 'zod';

// ============================================================================
// Video Configuration Schema
// ============================================================================

/**
 * Video configuration schema defining dimensions, frame rate, and duration
 */
export const VideoConfigSchema = z.object({
  /** Video width in pixels */
  width: z.number().int().positive().default(1920),
  /** Video height in pixels */
  height: z.number().int().positive().default(1080),
  /** Frames per second */
  fps: z.number().int().positive().default(30),
  /** Total duration in frames */
  durationInFrames: z.number().int().positive(),
});

export type VideoConfig = z.infer<typeof VideoConfigSchema>;

// ============================================================================
// Composition Schema
// ============================================================================

/**
 * Composition schema representing a video composition unit
 */
export const CompositionSchema = z.object({
  /** Unique identifier for the composition */
  id: z.string(),
  /** Display name of the composition */
  name: z.string(),
  /** React component code (file path or inline code) */
  code: z.string(),
  /** Props to pass to the composition component */
  props: z.record(z.string(), z.any()).default({}),
});

export type Composition = z.infer<typeof CompositionSchema>;

// ============================================================================
// Asset Schema
// ============================================================================

/**
 * Asset type enumeration
 */
export const AssetTypeEnum = z.enum(['image', 'video', 'audio', 'font']);

export type AssetType = z.infer<typeof AssetTypeEnum>;

/**
 * Asset schema representing media resources used in video projects
 */
export const AssetSchema = z.object({
  /** Unique identifier for the asset */
  id: z.string(),
  /** Type of the asset */
  type: AssetTypeEnum,
  /** Display name of the asset */
  name: z.string(),
  /** File path relative to project directory */
  path: z.string(),
});

export type Asset = z.infer<typeof AssetSchema>;

// ============================================================================
// Render History Schema
// ============================================================================

/**
 * Render status enumeration
 */
export const RenderStatusEnum = z.enum([
  'pending',
  'rendering',
  'completed',
  'failed',
]);

export type RenderStatus = z.infer<typeof RenderStatusEnum>;

/**
 * Render history schema tracking render operations
 */
export const RenderHistorySchema = z.object({
  /** Unique identifier for the render operation */
  id: z.string(),
  /** ID of the composition being rendered */
  compositionId: z.string(),
  /** Output file path */
  outputPath: z.string(),
  /** Current render status */
  status: RenderStatusEnum,
  /** Render progress (0-100) */
  progress: z.number().min(0).max(100),
  /** Timestamp when render was initiated */
  createdAt: z.string(),
  /** Error message if render failed */
  error: z.string().optional(),
});

export type RenderHistory = z.infer<typeof RenderHistorySchema>;

// ============================================================================
// Video Project Schema
// ============================================================================

/**
 * Complete video project schema
 */
export const VideoProjectSchema = z.object({
  /** Unique identifier for the project */
  id: z.string(),
  /** Project name */
  name: z.string(),
  /** Optional project description */
  description: z.string().optional(),
  /** Timestamp when project was created */
  createdAt: z.string(),
  /** Timestamp when project was last updated */
  updatedAt: z.string(),
  /** Video configuration */
  config: VideoConfigSchema,
  /** List of compositions in the project */
  compositions: z.array(CompositionSchema).default([]),
  /** List of assets used in the project */
  assets: z.array(AssetSchema).default([]),
  /** Render history */
  renders: z.array(RenderHistorySchema).default([]),
});

export type VideoProject = z.infer<typeof VideoProjectSchema>;

// ============================================================================
// Render Options
// ============================================================================

/**
 * Output format enumeration
 */
export const OutputFormatEnum = z.enum(['mp4', 'webm', 'gif']);

export type OutputFormat = z.infer<typeof OutputFormatEnum>;

/**
 * Quality preset enumeration
 */
export const QualityPresetEnum = z.enum(['draft', 'standard', 'high']);

export type QualityPreset = z.infer<typeof QualityPresetEnum>;

/**
 * Render options schema
 */
export const RenderOptionsSchema = z.object({
  /** Project ID to render */
  projectId: z.string(),
  /** Composition ID to render */
  compositionId: z.string(),
  /** Output format */
  outputFormat: OutputFormatEnum.default('mp4'),
  /** Quality preset */
  quality: QualityPresetEnum.default('standard'),
});

export type RenderOptions = z.infer<typeof RenderOptionsSchema>;

// ============================================================================
// Render Progress
// ============================================================================

/**
 * Render progress status enumeration
 */
export const RenderProgressStatusEnum = z.enum([
  'bundling',
  'preparing',
  'rendering',
  'completed',
  'failed',
]);

export type RenderProgressStatus = z.infer<typeof RenderProgressStatusEnum>;

/**
 * Render progress schema for real-time progress updates
 */
export const RenderProgressSchema = z.object({
  /** Current status */
  status: RenderProgressStatusEnum,
  /** Progress percentage (0-100) */
  progress: z.number().min(0).max(100),
  /** Error message if failed */
  error: z.string().optional(),
});

export type RenderProgress = z.infer<typeof RenderProgressSchema>;

// ============================================================================
// Quality Presets Configuration
// ============================================================================

/**
 * Quality preset configuration
 */
export interface QualityPresetConfig {
  /** Constant Rate Factor (lower = better quality, larger file) */
  crf: number;
  /** Scale factor (1 = original size) */
  scale: number;
  /** Frames per second */
  fps: number;
}

/**
 * Quality presets mapping
 */
export const QUALITY_PRESETS: Record<QualityPreset, QualityPresetConfig> = {
  draft: { crf: 28, scale: 0.5, fps: 15 },
  standard: { crf: 18, scale: 1, fps: 30 },
  high: { crf: 12, scale: 1, fps: 60 },
} as const;

// ============================================================================
// Supported Asset Types
// ============================================================================

/**
 * Supported file extensions for each asset type
 */
export const SUPPORTED_ASSET_EXTENSIONS: Record<AssetType, readonly string[]> =
  {
    image: ['.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'],
    video: ['.mp4', '.webm', '.mov'],
    audio: ['.mp3', '.wav', '.ogg', '.m4a'],
    font: ['.ttf', '.otf', '.woff', '.woff2'],
  } as const;
