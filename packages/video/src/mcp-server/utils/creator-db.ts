/**
 * Creator DB 轻量读取器
 *
 * 使用 better-sqlite3 只读访问 creator.db，获取活跃项目和内容序号信息。
 * 仅用于路径解析，不做写入操作。
 */

import Database from "better-sqlite3";
import { existsSync } from "fs";
import { join } from "path";

/** DB 中的项目记录 */
export interface CreatorProject {
  id: string;
  name: string;
  is_active: number;
}

/** DB 路径常量 */
const CREATOR_DB_RELATIVE_PATH = join(".sprouty-ai", "db", "creator.db");

/**
 * 获取 creator.db 的完整路径
 */
export function getCreatorDbPath(workspacePath: string): string {
  return join(workspacePath, CREATOR_DB_RELATIVE_PATH);
}

/**
 * 检查 creator.db 是否存在
 */
export function hasCreatorDb(workspacePath: string): boolean {
  return existsSync(getCreatorDbPath(workspacePath));
}

/**
 * 获取当前活跃项目
 *
 * @param dbPath creator.db 的完整路径
 * @returns 活跃项目，如果没有则返回 undefined
 */
export function getActiveProject(dbPath: string): CreatorProject | undefined {
  const db = new Database(dbPath, { readonly: true });
  try {
    const row = db
      .prepare("SELECT id, name, is_active FROM projects WHERE is_active = 1")
      .get() as CreatorProject | undefined;
    return row;
  } finally {
    db.close();
  }
}

/**
 * 获取项目的下一个内容序号
 *
 * 通过统计项目下已有内容数量 + 1 来确定下一个序号。
 * 序号从 1 开始，格式为两位数字（如 01, 02, ...）。
 *
 * @param dbPath creator.db 的完整路径
 * @param projectId 项目 ID
 * @returns 下一个内容序号（从 1 开始）
 */
export function getNextContentIndex(dbPath: string, projectId: string): number {
  const db = new Database(dbPath, { readonly: true });
  try {
    const row = db
      .prepare("SELECT COUNT(*) as count FROM contents WHERE project_id = ?")
      .get(projectId) as { count: number } | undefined;
    return (row?.count ?? 0) + 1;
  } finally {
    db.close();
  }
}

/**
 * 格式化内容序号为两位字符串
 *
 * @example formatContentIndex(1) => "01"
 * @example formatContentIndex(12) => "12"
 */
export function formatContentIndex(index: number): string {
  return String(index).padStart(2, "0");
}
