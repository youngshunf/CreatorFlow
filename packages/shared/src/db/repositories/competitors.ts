/**
 * 竞品账号 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { Competitor, CreateCompetitor, UpdateCompetitor } from '../types.ts';

/** 创建竞品 */
export function createCompetitor(db: CreatorMediaDB, data: CreateCompetitor): Competitor {
  const stmt = db.prepare(`
    INSERT INTO competitors (id, project_id, name, platform, url, follower_count, avg_likes, content_style, strengths, notes, tags, last_analyzed)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.project_id, data.name, data.platform,
    data.url, data.follower_count, data.avg_likes,
    data.content_style, data.strengths, data.notes,
    data.tags, data.last_analyzed,
  );
  return db.prepare<Competitor>('SELECT * FROM competitors WHERE id = ?').get(data.id)!;
}

/** 列出项目的所有竞品 */
export function listByProject(db: CreatorMediaDB, projectId: string): Competitor[] {
  return db.prepare<Competitor>(
    'SELECT * FROM competitors WHERE project_id = ? ORDER BY created_at DESC'
  ).all(projectId);
}

/** 更新竞品 */
export function updateCompetitor(db: CreatorMediaDB, id: string, data: UpdateCompetitor): Competitor | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) {
    return db.prepare<Competitor>('SELECT * FROM competitors WHERE id = ?').get(id);
  }

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE competitors SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return db.prepare<Competitor>('SELECT * FROM competitors WHERE id = ?').get(id);
}

/** 删除竞品 */
export function deleteCompetitor(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM competitors WHERE id = ?').run(id);
  return result.changes > 0;
}
