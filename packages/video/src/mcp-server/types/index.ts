/**
 * 类型定义
 *
 * MCP Video Server 的核心类型定义
 * 复用 @sprouty-ai/video 的类型，并添加 MCP 特定类型
 */

import { z } from "zod";

// ============================================================================
// 从 @sprouty-ai/video 重新导出核心类型
// ============================================================================

export {
  // Schemas
  VideoConfigSchema,
  SceneSchema,
  TransitionSchema,
  TransitionTypeEnum,
  TransitionDirectionEnum,
  AssetSchema,
  AssetTypeEnum,
  RenderHistorySchema,
  RenderStatusEnum,
  VideoProjectSchema,
  OutputFormatEnum,
  QualityPresetEnum,
  RenderOptionsSchema,
  RenderProgressStatusEnum,
  RenderProgressSchema,
  // Types
  type VideoConfig,
  type Scene,
  type Transition,
  type TransitionType,
  type TransitionDirection,
  type Asset,
  type AssetType,
  type RenderHistory,
  type RenderStatus,
  type VideoProject,
  type OutputFormat,
  type QualityPreset,
  type RenderOptions,
  type RenderProgressStatus,
  type RenderProgress,
  type QualityPresetConfig,
  // Constants
  QUALITY_PRESETS,
  SUPPORTED_ASSET_EXTENSIONS,
} from "../../types";

// ============================================================================
// MCP 错误代码枚举
// ============================================================================

/**
 * MCP 错误代码枚举
 * 用于标识不同类型的错误
 */
export const ErrorCodeEnum = z.enum([
  "PROJECT_NOT_FOUND",
  "ASSET_NOT_FOUND",
  "COMPOSITION_NOT_FOUND",
  "INVALID_INPUT",
  "FILE_NOT_FOUND",
  "UNSUPPORTED_FORMAT",
  "RENDER_FAILED",
  "PREVIEW_FAILED",
  "INTERNAL_ERROR",
]);

export type ErrorCode = z.infer<typeof ErrorCodeEnum>;

// ============================================================================
// MCP 响应类型
// ============================================================================

/**
 * 成功响应 Schema
 */
export const SuccessResponseSchema = <T extends z.ZodTypeAny>(dataSchema: T) =>
  z.object({
    success: z.literal(true),
    data: dataSchema,
  });

/**
 * 成功响应类型
 */
export interface SuccessResponse<T> {
  success: true;
  data: T;
}

/**
 * 错误详情 Schema
 */
export const ErrorDetailsSchema = z.object({
  /** 导致错误的字段 */
  field: z.string().optional(),
  /** 期望的值/格式 */
  expected: z.string().optional(),
  /** 实际收到的值 */
  received: z.string().optional(),
  /** 相关的文件路径 */
  path: z.string().optional(),
});

export type ErrorDetails = z.infer<typeof ErrorDetailsSchema>;

/**
 * 错误响应 Schema
 */
export const ErrorResponseSchema = z.object({
  success: z.literal(false),
  error: z.object({
    code: ErrorCodeEnum,
    message: z.string(),
    details: ErrorDetailsSchema.optional(),
  }),
});

/**
 * 错误响应类型
 */
export interface ErrorResponse {
  success: false;
  error: {
    code: ErrorCode;
    message: string;
    details?: ErrorDetails;
  };
}

/**
 * MCP 响应类型（成功或错误）
 */
export type MCPResponse<T> = SuccessResponse<T> | ErrorResponse;

// ============================================================================
// 项目摘要类型（用于列表查询）
// ============================================================================

/**
 * 项目摘要 Schema
 */
export const ProjectSummarySchema = z.object({
  /** 项目 ID */
  id: z.string(),
  /** 项目名称 */
  name: z.string(),
  /** 项目描述 */
  description: z.string().optional(),
  /** 创建时间 */
  createdAt: z.string(),
  /** 更新时间 */
  updatedAt: z.string(),
  /** 场景数量 */
  sceneCount: z.number().int().nonnegative(),
  /** 素材数量 */
  assetCount: z.number().int().nonnegative(),
});

export type ProjectSummary = z.infer<typeof ProjectSummarySchema>;

// ============================================================================
// 创建项目输入类型
// ============================================================================

/**
 * 创建项目输入 Schema
 */
export const CreateProjectInputSchema = z.object({
  /** 项目名称 */
  name: z.string().min(1, "项目名称不能为空"),
  /** 模板 ID（可选） */
  template: z.string().optional(),
  /** 视频宽度（可选，默认 1920） */
  width: z.number().int().positive().optional(),
  /** 视频高度（可选，默认 1080） */
  height: z.number().int().positive().optional(),
  /** 帧率（可选，默认 30） */
  fps: z.number().int().positive().optional(),
  /** 视频时长（秒） */
  durationInSeconds: z.number().positive(),
  /** 项目描述 */
  description: z.string().optional(),
});

export type CreateProjectInput = z.infer<typeof CreateProjectInputSchema>;

// ============================================================================
// 更新项目输入类型
// ============================================================================

/**
 * 更新项目输入 Schema
 */
export const UpdateProjectInputSchema = z.object({
  /** 项目名称 */
  name: z.string().min(1).optional(),
  /** 项目描述 */
  description: z.string().optional(),
  /** 视频配置 */
  config: z
    .object({
      width: z.number().int().positive().optional(),
      height: z.number().int().positive().optional(),
      fps: z.number().int().positive().optional(),
      durationInFrames: z.number().int().positive().optional(),
    })
    .optional(),
});

export type UpdateProjectInput = z.infer<typeof UpdateProjectInputSchema>;

// ============================================================================
// 预览服务器类型
// ============================================================================

/**
 * 预览实例 Schema
 */
export const PreviewInstanceSchema = z.object({
  /** 项目 ID */
  projectId: z.string(),
  /** 端口号 */
  port: z.number().int().positive(),
  /** 预览 URL */
  url: z.string().url(),
  /** 进程 ID */
  pid: z.number().int().positive(),
});

export type PreviewInstance = z.infer<typeof PreviewInstanceSchema>;

// ============================================================================
// 渲染结果类型
// ============================================================================

/**
 * 渲染结果 Schema
 */
export const RenderResultSchema = z.object({
  /** 是否成功 */
  success: z.boolean(),
  /** 输出文件路径 */
  outputPath: z.string(),
  /** 渲染耗时（毫秒） */
  duration: z.number().nonnegative(),
  /** 文件大小（字节） */
  fileSize: z.number().int().nonnegative(),
});

export type RenderResult = z.infer<typeof RenderResultSchema>;

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * 创建成功响应
 */
export function createSuccessResponse<T>(data: T): SuccessResponse<T> {
  return {
    success: true,
    data,
  };
}

/**
 * 创建错误响应
 */
export function createErrorResponse(
  code: ErrorCode,
  message: string,
  details?: ErrorDetails,
): ErrorResponse {
  return {
    success: false,
    error: {
      code,
      message,
      ...(details && { details }),
    },
  };
}
