/**
 * 内容目录路径工具
 *
 * 生成和管理内容的文件目录结构：
 * `项目名称/001_内容标题/视频/脚本|素材|作品`
 */

import { mkdirSync, existsSync } from 'node:fs';
import { join } from 'node:path';
import type { ContentTrack } from '../db/types.ts';

// ============================================================
// 文件名安全处理
// ============================================================

/** 非法文件名字符（Windows + macOS + Linux 通用） */
const ILLEGAL_CHARS = /[<>:"/\\|?*\x00-\x1f]/g;

/**
 * 将字符串处理为安全的文件/目录名
 * - 替换非法字符为下划线
 * - 空格转下划线
 * - 连续下划线合并
 * - 限长 50 字符
 */
export function sanitizeFileName(name: string): string {
  return name
    .replace(ILLEGAL_CHARS, '_')
    .replace(/\s+/g, '_')
    .replace(/_+/g, '_')
    .replace(/^_|_$/g, '')
    .slice(0, 50) || 'untitled';
}

// ============================================================
// 路径生成
// ============================================================

/**
 * 生成内容目录的相对路径
 *
 * @param projectName - 项目名称
 * @param contentNumber - 内容序号（从 1 开始）
 * @param contentTitle - 内容标题
 * @returns 相对路径，如 `我的项目/001_视频标题`
 */
export function generateContentDirPath(
  projectName: string,
  contentNumber: number,
  contentTitle: string,
): string {
  const safeProject = sanitizeFileName(projectName);
  const safeTitle = sanitizeFileName(contentTitle);
  const paddedNumber = String(contentNumber).padStart(3, '0');
  return join(safeProject, `${paddedNumber}_${safeTitle}`);
}

/** 轨道对应的子目录映射 */
const TRACK_SUBDIRS: Record<ContentTrack, string[]> = {
  video: ['视频/脚本', '视频/素材', '视频/作品'],
  article: ['文章/脚本', '文章/素材', '文章/作品'],
};

/**
 * 根据内容轨道返回需要创建的子目录列表
 *
 * @param contentTracks - 轨道字符串，如 'video' | 'article' | 'article,video'
 * @returns 子目录相对路径数组
 */
export function getContentSubDirs(contentTracks: string | null | undefined): string[] {
  const tracks = (contentTracks || 'article,video')
    .split(',')
    .filter((t): t is ContentTrack => t === 'article' || t === 'video');

  return tracks.flatMap(track => TRACK_SUBDIRS[track]);
}

/**
 * 在文件系统中创建内容目录结构
 *
 * @param workspaceRoot - 工作区根目录绝对路径
 * @param contentDirPath - 内容目录相对路径（由 generateContentDirPath 生成）
 * @param contentTracks - 轨道字符串
 */
export function ensureContentDirs(
  workspaceRoot: string,
  contentDirPath: string,
  contentTracks: string | null | undefined,
): void {
  const subDirs = getContentSubDirs(contentTracks);
  for (const sub of subDirs) {
    const fullPath = join(workspaceRoot, contentDirPath, sub);
    if (!existsSync(fullPath)) {
      mkdirSync(fullPath, { recursive: true });
    }
  }
}
