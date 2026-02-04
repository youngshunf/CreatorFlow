/**
 * 错误处理模块
 *
 * MCP 错误类和错误工厂函数
 * 提供结构化的错误处理，确保所有错误都有一致的格式
 */

import type { ErrorCode, ErrorDetails, ErrorResponse } from './index';

// ============================================================================
// MCPError 类
// ============================================================================

/**
 * MCP 错误类
 * 用于在 MCP Server 中抛出和处理错误
 */
export class MCPError extends Error {
  /** 错误代码 */
  public readonly code: ErrorCode;
  /** 错误详情 */
  public readonly details?: ErrorDetails;

  constructor(code: ErrorCode, message: string, details?: ErrorDetails) {
    super(message);
    this.name = 'MCPError';
    this.code = code;
    this.details = details;

    // 确保 instanceof 检查正常工作
    Object.setPrototypeOf(this, MCPError.prototype);
  }

  /**
   * 转换为 ErrorResponse 格式
   */
  toResponse(): ErrorResponse {
    return {
      success: false,
      error: {
        code: this.code,
        message: this.message,
        ...(this.details && { details: this.details }),
      },
    };
  }

  /**
   * 从 ErrorResponse 创建 MCPError
   */
  static fromResponse(response: ErrorResponse): MCPError {
    return new MCPError(
      response.error.code,
      response.error.message,
      response.error.details
    );
  }
}

// ============================================================================
// 错误工厂函数
// ============================================================================

/**
 * 创建项目未找到错误
 * @param projectId 项目 ID
 */
export function createProjectNotFoundError(projectId: string): MCPError {
  return new MCPError(
    'PROJECT_NOT_FOUND',
    `项目不存在: ${projectId}`,
    { field: 'projectId', received: projectId }
  );
}

/**
 * 创建素材未找到错误
 * @param assetId 素材 ID
 */
export function createAssetNotFoundError(assetId: string): MCPError {
  return new MCPError(
    'ASSET_NOT_FOUND',
    `素材不存在: ${assetId}`,
    { field: 'assetId', received: assetId }
  );
}

/**
 * 创建组合未找到错误
 * @param compositionId 组合 ID
 */
export function createCompositionNotFoundError(compositionId: string): MCPError {
  return new MCPError(
    'COMPOSITION_NOT_FOUND',
    `组合不存在: ${compositionId}`,
    { field: 'compositionId', received: compositionId }
  );
}

/**
 * 创建输入验证错误
 * @param message 错误消息
 * @param field 导致错误的字段
 * @param expected 期望的值/格式
 * @param received 实际收到的值
 */
export function createValidationError(
  message: string,
  field?: string,
  expected?: string,
  received?: string
): MCPError {
  return new MCPError('INVALID_INPUT', message, {
    ...(field && { field }),
    ...(expected && { expected }),
    ...(received && { received }),
  });
}

/**
 * 创建文件未找到错误
 * @param filePath 文件路径
 */
export function createFileNotFoundError(filePath: string): MCPError {
  return new MCPError(
    'FILE_NOT_FOUND',
    `文件不存在: ${filePath}`,
    { path: filePath }
  );
}

/**
 * 创建不支持的格式错误
 * @param format 收到的格式
 * @param supportedFormats 支持的格式列表
 */
export function createUnsupportedFormatError(
  format: string,
  supportedFormats: string[]
): MCPError {
  return new MCPError(
    'UNSUPPORTED_FORMAT',
    `不支持的格式: ${format}`,
    {
      received: format,
      expected: supportedFormats.join(', '),
    }
  );
}

/**
 * 创建渲染失败错误
 * @param message 错误消息
 * @param details 额外详情
 */
export function createRenderFailedError(
  message: string,
  details?: ErrorDetails
): MCPError {
  return new MCPError('RENDER_FAILED', `渲染失败: ${message}`, details);
}

/**
 * 创建预览失败错误
 * @param message 错误消息
 * @param details 额外详情
 */
export function createPreviewFailedError(
  message: string,
  details?: ErrorDetails
): MCPError {
  return new MCPError('PREVIEW_FAILED', `预览服务器启动失败: ${message}`, details);
}

/**
 * 创建内部错误
 * 用于捕获未预期的错误，不暴露内部实现细节
 * @param originalError 原始错误（仅用于日志记录）
 */
export function createInternalError(originalError?: unknown): MCPError {
  // 记录原始错误用于调试（在生产环境中应该记录到日志系统）
  if (originalError) {
    console.error('[MCP Video Server] Internal error:', originalError);
  }

  return new MCPError(
    'INTERNAL_ERROR',
    '服务器内部错误，请稍后重试'
  );
}

// ============================================================================
// 错误处理辅助函数
// ============================================================================

/**
 * 判断是否为 MCPError
 */
export function isMCPError(error: unknown): error is MCPError {
  return error instanceof MCPError;
}

/**
 * 将任意错误转换为 MCPError
 * 如果已经是 MCPError，直接返回
 * 否则包装为内部错误
 */
export function toMCPError(error: unknown): MCPError {
  if (isMCPError(error)) {
    return error;
  }

  return createInternalError(error);
}

/**
 * 将任意错误转换为 ErrorResponse
 * 确保所有错误都以一致的格式返回
 */
export function toErrorResponse(error: unknown): ErrorResponse {
  return toMCPError(error).toResponse();
}

/**
 * 安全执行函数，捕获错误并转换为 ErrorResponse
 * @param fn 要执行的函数
 */
export async function safeExecute<T>(
  fn: () => Promise<T>
): Promise<{ success: true; data: T } | ErrorResponse> {
  try {
    const data = await fn();
    return { success: true, data };
  } catch (error) {
    return toErrorResponse(error);
  }
}

/**
 * 同步版本的安全执行函数
 * @param fn 要执行的函数
 */
export function safeExecuteSync<T>(
  fn: () => T
): { success: true; data: T } | ErrorResponse {
  try {
    const data = fn();
    return { success: true, data };
  } catch (error) {
    return toErrorResponse(error);
  }
}
