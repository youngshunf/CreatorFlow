/**
 * 验证器工具函数
 *
 * 项目名称、素材类型、工作区路径等验证
 * 提供统一的输入验证，确保数据完整性和安全性
 */

import { existsSync, statSync } from 'fs';
import { extname, isAbsolute } from 'path';
import { SUPPORTED_ASSET_EXTENSIONS, type AssetType } from '../types';
import { createValidationError, createFileNotFoundError, createUnsupportedFormatError } from '../types/errors';

// ============================================================================
// 常量定义
// ============================================================================

/**
 * 项目名称中不允许的字符
 * 包括文件系统保留字符和可能导致问题的特殊字符
 */
const INVALID_PROJECT_NAME_CHARS = /[<>:"/\\|?*\x00-\x1f]/;

/**
 * 项目名称最大长度
 * 考虑到文件系统限制和用户体验
 */
const MAX_PROJECT_NAME_LENGTH = 100;

/**
 * 项目名称最小长度
 */
const MIN_PROJECT_NAME_LENGTH = 1;

// ============================================================================
// 项目名称验证
// ============================================================================

/**
 * 验证项目名称是否有效
 * @param name 项目名称
 * @returns 验证结果对象
 *
 * 验证规则:
 * - 不能为空或仅包含空白字符
 * - 长度在 1-100 字符之间
 * - 不能包含文件系统保留字符 (<>:"/\|?*)
 * - 不能以点号或空格开头/结尾
 * - 不能是保留名称（如 CON, PRN, AUX 等 Windows 保留名）
 *
 * @example
 * validateProjectName('我的视频项目')
 * // => { valid: true }
 *
 * validateProjectName('')
 * // => { valid: false, error: MCPError }
 */
export function validateProjectName(name: string): { valid: true } | { valid: false; error: ReturnType<typeof createValidationError> } {
  // 检查是否为空
  if (!name || typeof name !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '项目名称不能为空',
        'name',
        '非空字符串',
        String(name)
      ),
    };
  }

  // 去除首尾空白后检查
  const trimmedName = name.trim();
  if (trimmedName.length === 0) {
    return {
      valid: false,
      error: createValidationError(
        '项目名称不能仅包含空白字符',
        'name',
        '非空字符串',
        '空白字符串'
      ),
    };
  }

  // 检查长度
  if (trimmedName.length < MIN_PROJECT_NAME_LENGTH) {
    return {
      valid: false,
      error: createValidationError(
        `项目名称长度不能少于 ${MIN_PROJECT_NAME_LENGTH} 个字符`,
        'name',
        `至少 ${MIN_PROJECT_NAME_LENGTH} 个字符`,
        `${trimmedName.length} 个字符`
      ),
    };
  }

  if (trimmedName.length > MAX_PROJECT_NAME_LENGTH) {
    return {
      valid: false,
      error: createValidationError(
        `项目名称长度不能超过 ${MAX_PROJECT_NAME_LENGTH} 个字符`,
        'name',
        `最多 ${MAX_PROJECT_NAME_LENGTH} 个字符`,
        `${trimmedName.length} 个字符`
      ),
    };
  }

  // 检查无效字符
  if (INVALID_PROJECT_NAME_CHARS.test(trimmedName)) {
    return {
      valid: false,
      error: createValidationError(
        '项目名称包含无效字符',
        'name',
        '不包含 <>:"/\\|?* 等特殊字符',
        trimmedName
      ),
    };
  }

  // 检查是否以点号或空格开头/结尾
  if (trimmedName.startsWith('.') || trimmedName.endsWith('.')) {
    return {
      valid: false,
      error: createValidationError(
        '项目名称不能以点号开头或结尾',
        'name',
        '不以点号开头或结尾',
        trimmedName
      ),
    };
  }

  // 检查 Windows 保留名称
  const reservedNames = [
    'CON', 'PRN', 'AUX', 'NUL',
    'COM1', 'COM2', 'COM3', 'COM4', 'COM5', 'COM6', 'COM7', 'COM8', 'COM9',
    'LPT1', 'LPT2', 'LPT3', 'LPT4', 'LPT5', 'LPT6', 'LPT7', 'LPT8', 'LPT9',
  ];
  if (reservedNames.includes(trimmedName.toUpperCase())) {
    return {
      valid: false,
      error: createValidationError(
        '项目名称不能使用系统保留名称',
        'name',
        '非系统保留名称',
        trimmedName
      ),
    };
  }

  return { valid: true };
}

/**
 * 验证项目名称，如果无效则抛出错误
 * @param name 项目名称
 * @throws MCPError 如果名称无效
 */
export function assertValidProjectName(name: string): void {
  const result = validateProjectName(name);
  if (!result.valid) {
    throw result.error;
  }
}

// ============================================================================
// 工作区路径验证
// ============================================================================

/**
 * 验证工作区路径是否有效
 * @param workspacePath 工作区路径
 * @returns 验证结果对象
 *
 * 验证规则:
 * - 不能为空
 * - 必须是绝对路径
 * - 路径必须存在
 * - 必须是目录
 *
 * @example
 * validateWorkspacePath('/Users/user/workspace')
 * // => { valid: true }
 *
 * validateWorkspacePath('./relative/path')
 * // => { valid: false, error: MCPError }
 */
export function validateWorkspacePath(workspacePath: string): { valid: true } | { valid: false; error: ReturnType<typeof createValidationError> | ReturnType<typeof createFileNotFoundError> } {
  // 检查是否为空
  if (!workspacePath || typeof workspacePath !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '工作区路径不能为空',
        'workspacePath',
        '非空字符串',
        String(workspacePath)
      ),
    };
  }

  const trimmedPath = workspacePath.trim();
  if (trimmedPath.length === 0) {
    return {
      valid: false,
      error: createValidationError(
        '工作区路径不能为空',
        'workspacePath',
        '非空字符串',
        '空字符串'
      ),
    };
  }

  // 检查是否为绝对路径
  if (!isAbsolute(trimmedPath)) {
    return {
      valid: false,
      error: createValidationError(
        '工作区路径必须是绝对路径',
        'workspacePath',
        '绝对路径（如 /Users/user/workspace）',
        trimmedPath
      ),
    };
  }

  // 检查路径是否存在
  if (!existsSync(trimmedPath)) {
    return {
      valid: false,
      error: createFileNotFoundError(trimmedPath),
    };
  }

  // 检查是否为目录
  try {
    const stats = statSync(trimmedPath);
    if (!stats.isDirectory()) {
      return {
        valid: false,
        error: createValidationError(
          '工作区路径必须是目录',
          'workspacePath',
          '目录路径',
          '文件路径'
        ),
      };
    }
  } catch {
    return {
      valid: false,
      error: createValidationError(
        '无法访问工作区路径',
        'workspacePath',
        '可访问的目录',
        trimmedPath
      ),
    };
  }

  return { valid: true };
}

/**
 * 验证工作区路径，如果无效则抛出错误
 * @param workspacePath 工作区路径
 * @throws MCPError 如果路径无效
 */
export function assertValidWorkspacePath(workspacePath: string): void {
  const result = validateWorkspacePath(workspacePath);
  if (!result.valid) {
    throw result.error;
  }
}

// ============================================================================
// 素材类型验证
// ============================================================================

/**
 * 验证素材类型是否有效
 * @param assetType 素材类型
 * @returns 验证结果对象
 *
 * @example
 * validateAssetType('image')
 * // => { valid: true }
 *
 * validateAssetType('unknown')
 * // => { valid: false, error: MCPError }
 */
export function validateAssetType(assetType: string): { valid: true; type: AssetType } | { valid: false; error: ReturnType<typeof createValidationError> } {
  const validTypes: AssetType[] = ['image', 'video', 'audio', 'font'];

  if (!validTypes.includes(assetType as AssetType)) {
    return {
      valid: false,
      error: createValidationError(
        `无效的素材类型: ${assetType}`,
        'assetType',
        validTypes.join(', '),
        assetType
      ),
    };
  }

  return { valid: true, type: assetType as AssetType };
}

/**
 * 验证素材文件扩展名是否支持
 * @param filePath 文件路径
 * @param assetType 素材类型
 * @returns 验证结果对象
 *
 * @example
 * validateAssetExtension('/path/to/image.png', 'image')
 * // => { valid: true, extension: '.png' }
 *
 * validateAssetExtension('/path/to/file.txt', 'image')
 * // => { valid: false, error: MCPError }
 */
export function validateAssetExtension(
  filePath: string,
  assetType: AssetType
): { valid: true; extension: string } | { valid: false; error: ReturnType<typeof createUnsupportedFormatError> } {
  const extension = extname(filePath).toLowerCase();
  const supportedExtensions = SUPPORTED_ASSET_EXTENSIONS[assetType];

  if (!supportedExtensions.includes(extension as any)) {
    return {
      valid: false,
      error: createUnsupportedFormatError(
        extension || '(无扩展名)',
        [...supportedExtensions]
      ),
    };
  }

  return { valid: true, extension };
}

/**
 * 从文件扩展名推断素材类型
 * @param filePath 文件路径
 * @returns 素材类型，如果无法推断则返回 null
 *
 * @example
 * inferAssetType('/path/to/image.png')
 * // => 'image'
 *
 * inferAssetType('/path/to/unknown.xyz')
 * // => null
 */
export function inferAssetType(filePath: string): AssetType | null {
  const extension = extname(filePath).toLowerCase();

  for (const [type, extensions] of Object.entries(SUPPORTED_ASSET_EXTENSIONS)) {
    if ((extensions as readonly string[]).includes(extension)) {
      return type as AssetType;
    }
  }

  return null;
}

/**
 * 验证素材文件是否存在且格式正确
 * @param filePath 文件路径
 * @param assetType 素材类型
 * @returns 验证结果对象
 *
 * @example
 * validateAssetFile('/path/to/image.png', 'image')
 * // => { valid: true, extension: '.png' }
 */
export function validateAssetFile(
  filePath: string,
  assetType: AssetType
): { valid: true; extension: string } | { valid: false; error: ReturnType<typeof createValidationError> | ReturnType<typeof createFileNotFoundError> | ReturnType<typeof createUnsupportedFormatError> } {
  // 检查路径是否为空
  if (!filePath || typeof filePath !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '素材文件路径不能为空',
        'assetPath',
        '非空字符串',
        String(filePath)
      ),
    };
  }

  // 检查文件是否存在
  if (!existsSync(filePath)) {
    return {
      valid: false,
      error: createFileNotFoundError(filePath),
    };
  }

  // 检查是否为文件
  try {
    const stats = statSync(filePath);
    if (!stats.isFile()) {
      return {
        valid: false,
        error: createValidationError(
          '素材路径必须是文件',
          'assetPath',
          '文件路径',
          '目录路径'
        ),
      };
    }
  } catch {
    return {
      valid: false,
      error: createValidationError(
        '无法访问素材文件',
        'assetPath',
        '可访问的文件',
        filePath
      ),
    };
  }

  // 验证扩展名
  return validateAssetExtension(filePath, assetType);
}

// ============================================================================
// 通用验证辅助函数
// ============================================================================

/**
 * 验证字符串是否为有效的 UUID v4
 * @param id 要验证的字符串
 * @returns 是否为有效的 UUID
 */
export function isValidUUID(id: string): boolean {
  const uuidRegex = /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/i;
  return uuidRegex.test(id);
}

/**
 * 验证项目 ID 是否有效
 * @param projectId 项目 ID
 * @returns 验证结果对象
 */
export function validateProjectId(projectId: string): { valid: true } | { valid: false; error: ReturnType<typeof createValidationError> } {
  if (!projectId || typeof projectId !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '项目 ID 不能为空',
        'projectId',
        '非空字符串',
        String(projectId)
      ),
    };
  }

  if (!isValidUUID(projectId)) {
    return {
      valid: false,
      error: createValidationError(
        '项目 ID 格式无效',
        'projectId',
        'UUID v4 格式',
        projectId
      ),
    };
  }

  return { valid: true };
}

/**
 * 验证组合 ID 是否有效
 * @param compositionId 组合 ID
 * @returns 验证结果对象
 */
export function validateCompositionId(compositionId: string): { valid: true } | { valid: false; error: ReturnType<typeof createValidationError> } {
  if (!compositionId || typeof compositionId !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '组合 ID 不能为空',
        'compositionId',
        '非空字符串',
        String(compositionId)
      ),
    };
  }

  // 组合 ID 可以是 UUID 或简单的标识符
  const validIdRegex = /^[a-zA-Z0-9_-]+$/;
  if (!validIdRegex.test(compositionId)) {
    return {
      valid: false,
      error: createValidationError(
        '组合 ID 格式无效',
        'compositionId',
        '字母、数字、下划线或连字符',
        compositionId
      ),
    };
  }

  return { valid: true };
}

/**
 * 验证素材 ID 是否有效
 * @param assetId 素材 ID
 * @returns 验证结果对象
 */
export function validateAssetId(assetId: string): { valid: true } | { valid: false; error: ReturnType<typeof createValidationError> } {
  if (!assetId || typeof assetId !== 'string') {
    return {
      valid: false,
      error: createValidationError(
        '素材 ID 不能为空',
        'assetId',
        '非空字符串',
        String(assetId)
      ),
    };
  }

  if (!isValidUUID(assetId)) {
    return {
      valid: false,
      error: createValidationError(
        '素材 ID 格式无效',
        'assetId',
        'UUID v4 格式',
        assetId
      ),
    };
  }

  return { valid: true };
}
