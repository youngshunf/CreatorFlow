/**
 * 路径工具函数
 *
 * 处理视频项目相关的路径生成
 * 所有路径使用中文命名，方便用户直接查看和编辑
 *
 * 文件系统结构:
 * {Workspace_Root}/
 * └── 视频创作/
 *     └── {项目名称}/
 *         ├── project.json          # 项目配置
 *         ├── 素材/                  # 素材文件
 *         │   ├── images/
 *         │   ├── videos/
 *         │   ├── audio/
 *         │   └── fonts/
 *         ├── 组合/                  # 组合代码
 *         │   └── {compositionId}.tsx
 *         └── 输出/                  # 渲染输出
 *             └── {renderName}.mp4
 */

import { join } from 'path';

// ============================================================================
// 常量定义
// ============================================================================

/** 视频项目根目录名称 */
export const VIDEO_PROJECTS_DIR_NAME = '视频创作';

/** 素材目录名称 */
export const ASSETS_DIR_NAME = '素材';

/** 组合目录名称 */
export const COMPOSITIONS_DIR_NAME = '组合';

/** 输出目录名称 */
export const OUTPUT_DIR_NAME = '输出';

/** 项目配置文件名 */
export const PROJECT_CONFIG_FILE = 'project.json';

/** 素材子目录名称 */
export const ASSET_SUBDIRS = {
  images: 'images',
  videos: 'videos',
  audio: 'audio',
  fonts: 'fonts',
} as const;

// ============================================================================
// 路径生成函数
// ============================================================================

/**
 * 获取视频项目根目录路径
 * @param workspacePath 工作区根路径
 * @returns 视频项目目录路径 `{workspacePath}/视频创作/`
 *
 * @example
 * getVideoProjectsDir('/Users/user/workspace')
 * // => '/Users/user/workspace/视频创作'
 */
export function getVideoProjectsDir(workspacePath: string): string {
  return join(workspacePath, VIDEO_PROJECTS_DIR_NAME);
}

/**
 * 获取项目目录路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @returns 项目目录路径 `{workspacePath}/视频创作/{projectName}/`
 *
 * @example
 * getProjectPath('/Users/user/workspace', '我的视频项目')
 * // => '/Users/user/workspace/视频创作/我的视频项目'
 */
export function getProjectPath(workspacePath: string, projectName: string): string {
  return join(getVideoProjectsDir(workspacePath), projectName);
}

/**
 * 获取项目配置文件路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @returns 项目配置文件路径 `{workspacePath}/视频创作/{projectName}/project.json`
 *
 * @example
 * getProjectConfigPath('/Users/user/workspace', '我的视频项目')
 * // => '/Users/user/workspace/视频创作/我的视频项目/project.json'
 */
export function getProjectConfigPath(workspacePath: string, projectName: string): string {
  return join(getProjectPath(workspacePath, projectName), PROJECT_CONFIG_FILE);
}

/**
 * 获取素材目录路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @returns 素材目录路径 `{workspacePath}/视频创作/{projectName}/素材/`
 *
 * @example
 * getAssetsPath('/Users/user/workspace', '我的视频项目')
 * // => '/Users/user/workspace/视频创作/我的视频项目/素材'
 */
export function getAssetsPath(workspacePath: string, projectName: string): string {
  return join(getProjectPath(workspacePath, projectName), ASSETS_DIR_NAME);
}

/**
 * 获取特定类型素材的子目录路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @param assetType 素材类型
 * @returns 素材子目录路径 `{workspacePath}/视频创作/{projectName}/素材/{type}/`
 *
 * @example
 * getAssetTypePath('/Users/user/workspace', '我的视频项目', 'image')
 * // => '/Users/user/workspace/视频创作/我的视频项目/素材/images'
 */
export function getAssetTypePath(
  workspacePath: string,
  projectName: string,
  assetType: 'image' | 'video' | 'audio' | 'font'
): string {
  const subdir = assetType === 'image' ? ASSET_SUBDIRS.images
    : assetType === 'video' ? ASSET_SUBDIRS.videos
    : assetType === 'audio' ? ASSET_SUBDIRS.audio
    : ASSET_SUBDIRS.fonts;
  return join(getAssetsPath(workspacePath, projectName), subdir);
}

/**
 * 获取组合目录路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @returns 组合目录路径 `{workspacePath}/视频创作/{projectName}/组合/`
 *
 * @example
 * getCompositionsPath('/Users/user/workspace', '我的视频项目')
 * // => '/Users/user/workspace/视频创作/我的视频项目/组合'
 */
export function getCompositionsPath(workspacePath: string, projectName: string): string {
  return join(getProjectPath(workspacePath, projectName), COMPOSITIONS_DIR_NAME);
}

/**
 * 获取组合文件路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @param compositionId 组合 ID
 * @returns 组合文件路径 `{workspacePath}/视频创作/{projectName}/组合/{compositionId}.tsx`
 *
 * @example
 * getCompositionFilePath('/Users/user/workspace', '我的视频项目', 'main')
 * // => '/Users/user/workspace/视频创作/我的视频项目/组合/main.tsx'
 */
export function getCompositionFilePath(
  workspacePath: string,
  projectName: string,
  compositionId: string
): string {
  return join(getCompositionsPath(workspacePath, projectName), `${compositionId}.tsx`);
}

/**
 * 获取输出目录路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @returns 输出目录路径 `{workspacePath}/视频创作/{projectName}/输出/`
 *
 * @example
 * getOutputPath('/Users/user/workspace', '我的视频项目')
 * // => '/Users/user/workspace/视频创作/我的视频项目/输出'
 */
export function getOutputPath(workspacePath: string, projectName: string): string {
  return join(getProjectPath(workspacePath, projectName), OUTPUT_DIR_NAME);
}

/**
 * 获取渲染输出文件路径
 * @param workspacePath 工作区根路径
 * @param projectName 项目名称
 * @param fileName 输出文件名（包含扩展名）
 * @returns 输出文件路径 `{workspacePath}/视频创作/{projectName}/输出/{fileName}`
 *
 * @example
 * getOutputFilePath('/Users/user/workspace', '我的视频项目', 'render-001.mp4')
 * // => '/Users/user/workspace/视频创作/我的视频项目/输出/render-001.mp4'
 */
export function getOutputFilePath(
  workspacePath: string,
  projectName: string,
  fileName: string
): string {
  return join(getOutputPath(workspacePath, projectName), fileName);
}

// ============================================================================
// 路径解析函数
// ============================================================================

/**
 * 从项目路径中提取项目名称
 * @param projectPath 项目完整路径
 * @returns 项目名称，如果路径无效则返回 null
 *
 * @example
 * extractProjectName('/Users/user/workspace/视频创作/我的视频项目')
 * // => '我的视频项目'
 */
export function extractProjectName(projectPath: string): string | null {
  const parts = projectPath.split('/').filter(Boolean);
  const videoIndex = parts.indexOf(VIDEO_PROJECTS_DIR_NAME);
  if (videoIndex === -1 || videoIndex >= parts.length - 1) {
    return null;
  }
  return parts[videoIndex + 1] ?? null;
}

/**
 * 获取素材在项目中的相对路径
 * @param assetType 素材类型
 * @param fileName 文件名
 * @returns 相对于项目目录的路径
 *
 * @example
 * getAssetRelativePath('image', 'background.png')
 * // => '素材/images/background.png'
 */
export function getAssetRelativePath(
  assetType: 'image' | 'video' | 'audio' | 'font',
  fileName: string
): string {
  const subdir = assetType === 'image' ? ASSET_SUBDIRS.images
    : assetType === 'video' ? ASSET_SUBDIRS.videos
    : assetType === 'audio' ? ASSET_SUBDIRS.audio
    : ASSET_SUBDIRS.fonts;
  return join(ASSETS_DIR_NAME, subdir, fileName);
}

/**
 * 获取组合在项目中的相对路径
 * @param compositionId 组合 ID
 * @returns 相对于项目目录的路径
 *
 * @example
 * getCompositionRelativePath('main')
 * // => '组合/main.tsx'
 */
export function getCompositionRelativePath(compositionId: string): string {
  return join(COMPOSITIONS_DIR_NAME, `${compositionId}.tsx`);
}

/**
 * 获取输出文件在项目中的相对路径
 * @param fileName 文件名
 * @returns 相对于项目目录的路径
 *
 * @example
 * getOutputRelativePath('render-001.mp4')
 * // => '输出/render-001.mp4'
 */
export function getOutputRelativePath(fileName: string): string {
  return join(OUTPUT_DIR_NAME, fileName);
}
