/**
 * 内容版本管理 Repository — CRUD + 回滚操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { ContentVersion, CreateContentVersionInput, VersionStage, ChangeSource } from '../types.ts';

// 重新导出类型供外部使用
export type { ContentVersion, CreateContentVersionInput, VersionStage, ChangeSource };

/** 创建版本（自动计算 version_number） */
export function createVersion(db: CreatorMediaDB, data: CreateContentVersionInput): ContentVersion {
  // 自动计算下一个版本号
  const row = db.prepare<{ next: number }>(
    'SELECT COALESCE(MAX(version_number), 0) + 1 AS next FROM content_versions WHERE content_id = ?'
  ).get(data.content_id);
  const versionNumber = row?.next ?? 1;

  const stmt = db.prepare(`
    INSERT INTO content_versions (id, content_id, version_number, stage, title, content_snapshot, files_snapshot, change_source, change_description, created_by)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.content_id, versionNumber, data.stage,
    data.title, data.content_snapshot, data.files_snapshot,
    data.change_source, data.change_description, data.created_by,
  );
  return getVersion(db, data.id)!;
}

/** 获取单个版本 */
export function getVersion(db: CreatorMediaDB, id: string): ContentVersion | undefined {
  return db.prepare<ContentVersion>('SELECT * FROM content_versions WHERE id = ?').get(id);
}

/** 按内容查询所有版本（version_number DESC） */
export function listByContent(db: CreatorMediaDB, contentId: string): ContentVersion[] {
  return db.prepare<ContentVersion>(
    'SELECT * FROM content_versions WHERE content_id = ? ORDER BY version_number DESC'
  ).all(contentId);
}

/** 获取最新版本 */
export function getLatestVersion(db: CreatorMediaDB, contentId: string): ContentVersion | undefined {
  return db.prepare<ContentVersion>(
    'SELECT * FROM content_versions WHERE content_id = ? ORDER BY version_number DESC LIMIT 1'
  ).get(contentId);
}

/** 按版本号查询 */
export function getVersionByNumber(db: CreatorMediaDB, contentId: string, versionNumber: number): ContentVersion | undefined {
  return db.prepare<ContentVersion>(
    'SELECT * FROM content_versions WHERE content_id = ? AND version_number = ?'
  ).get(contentId, versionNumber);
}

/** 回滚到指定版本（创建新版本，复制目标版本的快照） */
export function rollbackToVersion(db: CreatorMediaDB, contentId: string, targetVersionNumber: number): ContentVersion {
  const target = getVersionByNumber(db, contentId, targetVersionNumber);
  if (!target) {
    throw new Error(`版本 ${targetVersionNumber} 不存在（content_id=${contentId}）`);
  }

  const id = crypto.randomUUID();
  return createVersion(db, {
    id,
    content_id: contentId,
    version_number: 0, // 由 createVersion 自动计算
    stage: target.stage,
    title: target.title,
    content_snapshot: target.content_snapshot,
    files_snapshot: target.files_snapshot,
    change_source: 'rollback',
    change_description: `回滚到版本 ${targetVersionNumber}`,
    created_by: target.created_by,
  });
}

/** 删除某内容的所有版本 */
export function deleteByContent(db: CreatorMediaDB, contentId: string): number {
  const result = db.prepare('DELETE FROM content_versions WHERE content_id = ?').run(contentId);
  return result.changes;
}

/** 统计某内容的版本数 */
export function countByContent(db: CreatorMediaDB, contentId: string): number {
  const row = db.prepare<{ count: number }>(
    'SELECT COUNT(*) AS count FROM content_versions WHERE content_id = ?'
  ).get(contentId);
  return row?.count ?? 0;
}
