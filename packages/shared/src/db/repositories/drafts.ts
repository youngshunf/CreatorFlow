/**
 * 草稿 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { Draft, CreateDraft, UpdateDraft } from '../types.ts';

/** 创建草稿 */
export function createDraft(db: CreatorMediaDB, data: CreateDraft): Draft {
  const stmt = db.prepare(`
    INSERT INTO drafts (id, project_id, title, content, media, tags, target_platforms, metadata)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.project_id,
    data.title,
    data.content,
    data.media ?? '[]',
    data.tags,
    data.target_platforms,
    data.metadata,
  );
  return getDraft(db, data.id)!;
}

/** 获取单个草稿 */
export function getDraft(db: CreatorMediaDB, id: string): Draft | undefined {
  return db.prepare<Draft>('SELECT * FROM drafts WHERE id = ?').get(id);
}

/** 列出项目的草稿 */
export function listByProject(db: CreatorMediaDB, projectId: string): Draft[] {
  return db.prepare<Draft>(
    'SELECT * FROM drafts WHERE project_id = ? ORDER BY updated_at DESC'
  ).all(projectId);
}

/** 更新草稿 */
export function updateDraft(db: CreatorMediaDB, id: string, data: UpdateDraft): Draft | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getDraft(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE drafts SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getDraft(db, id);
}

/** 删除草稿 */
export function deleteDraft(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM drafts WHERE id = ?').run(id);
  return result.changes > 0;
}
