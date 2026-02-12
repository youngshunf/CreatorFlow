/**
 * 内容阶段产出 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { ContentStageRecord, CreateContentStage, UpdateContentStage, ContentStage } from '../types.ts';

/** 创建阶段产出记录 */
export function createContentStage(db: CreatorMediaDB, data: CreateContentStage): ContentStageRecord {
  const stmt = db.prepare(`
    INSERT INTO content_stages (id, content_id, stage, file_path, status, version, source_type, metadata)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.content_id,
    data.stage,
    data.file_path,
    data.status ?? 'draft',
    data.version ?? 1,
    data.source_type,
    data.metadata,
  );
  return getContentStage(db, data.id)!;
}

/** 获取单个阶段产出记录 */
export function getContentStage(db: CreatorMediaDB, id: string): ContentStageRecord | undefined {
  return db.prepare<ContentStageRecord>('SELECT * FROM content_stages WHERE id = ?').get(id);
}

/** 获取内容的所有阶段产出 */
export function listByContent(db: CreatorMediaDB, contentId: string): ContentStageRecord[] {
  return db.prepare<ContentStageRecord>(
    'SELECT * FROM content_stages WHERE content_id = ? ORDER BY created_at DESC'
  ).all(contentId);
}

/** 获取内容的指定阶段的所有版本 */
export function listByContentAndStage(
  db: CreatorMediaDB,
  contentId: string,
  stage: ContentStage
): ContentStageRecord[] {
  return db.prepare<ContentStageRecord>(
    'SELECT * FROM content_stages WHERE content_id = ? AND stage = ? ORDER BY version DESC'
  ).all(contentId, stage);
}

/** 获取内容的指定阶段的最新版本 */
export function getLatestStage(
  db: CreatorMediaDB,
  contentId: string,
  stage: ContentStage
): ContentStageRecord | undefined {
  return db.prepare<ContentStageRecord>(
    'SELECT * FROM content_stages WHERE content_id = ? AND stage = ? ORDER BY version DESC LIMIT 1'
  ).get(contentId, stage);
}

/** 获取内容的指定阶段的下一个版本号 */
export function getNextVersionNumber(
  db: CreatorMediaDB,
  contentId: string,
  stage: ContentStage
): number {
  const result = db.prepare<{ max_version: number | null }>(
    'SELECT MAX(version) as max_version FROM content_stages WHERE content_id = ? AND stage = ?'
  ).get(contentId, stage);
  return (result?.max_version ?? 0) + 1;
}

/** 更新阶段产出记录 */
export function updateContentStage(
  db: CreatorMediaDB,
  id: string,
  data: UpdateContentStage
): ContentStageRecord | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getContentStage(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE content_stages SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getContentStage(db, id);
}

/** 删除阶段产出记录 */
export function deleteContentStage(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM content_stages WHERE id = ?').run(id);
  return result.changes > 0;
}

/** 删除内容的所有阶段产出 */
export function deleteByContent(db: CreatorMediaDB, contentId: string): number {
  const result = db.prepare('DELETE FROM content_stages WHERE content_id = ?').run(contentId);
  return result.changes;
}

/** 解析阶段产出的 metadata JSON */
export function parseStageMetadata(stage: ContentStageRecord): Record<string, unknown> | null {
  if (!stage.metadata) return null;
  try {
    return JSON.parse(stage.metadata) as Record<string, unknown>;
  } catch {
    return null;
  }
}
