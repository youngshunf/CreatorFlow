/**
 * 选题推荐 Repository — CRUD 操作（精简版）
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { RecommendedTopic, CreateRecommendedTopic, TopicRecommendStatus } from '../types.ts';

/** 精简字段列表 */
const SELECT_COLUMNS = 'id, project_id, title, reason, potential_score, heat_index, status, content_id, md_file_path, source_uid, batch_date, created_at, updated_at';

/** 选题列表过滤条件 */
export interface RecommendedTopicFilters {
  status?: TopicRecommendStatus;
  batchDate?: string;
  limit?: number;
}

/** 创建单条选题 */
export function create(db: CreatorMediaDB, topic: CreateRecommendedTopic): RecommendedTopic {
  const now = new Date().toISOString();
  db.prepare(`
    INSERT INTO recommended_topics (id, project_id, title, reason, potential_score, heat_index, status, content_id, md_file_path, source_uid, batch_date)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `).run(
    topic.id, topic.project_id, topic.title, topic.reason,
    topic.potential_score, topic.heat_index, topic.status,
    topic.content_id, topic.md_file_path, topic.source_uid, topic.batch_date,
  );
  return { ...topic, created_at: now, updated_at: null };
}

/** 批量创建选题（source_uid 去重） */
export function batchCreate(db: CreatorMediaDB, topics: CreateRecommendedTopic[]): void {
  const stmt = db.prepare(`
    INSERT OR IGNORE INTO recommended_topics (id, project_id, title, reason, potential_score, heat_index, status, content_id, md_file_path, source_uid, batch_date)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);

  db.transaction(() => {
    for (const t of topics) {
      stmt.run(
        t.id, t.project_id, t.title, t.reason,
        t.potential_score, t.heat_index, t.status,
        t.content_id, t.md_file_path, t.source_uid, t.batch_date,
      );
    }
  });
}

/** 获取单条选题 */
export function get(db: CreatorMediaDB, id: string): RecommendedTopic | undefined {
  return db.prepare<RecommendedTopic>(`SELECT ${SELECT_COLUMNS} FROM recommended_topics WHERE id = ?`).get(id);
}

/** 按项目查询选题列表 */
export function listByProject(db: CreatorMediaDB, projectId: string, filters?: RecommendedTopicFilters): RecommendedTopic[] {
  let sql = `SELECT ${SELECT_COLUMNS} FROM recommended_topics WHERE project_id = ?`;
  const params: unknown[] = [projectId];

  if (filters?.status !== undefined) {
    sql += ' AND status = ?';
    params.push(filters.status);
  }

  if (filters?.batchDate) {
    sql += ' AND batch_date = ?';
    params.push(filters.batchDate);
  }

  sql += ' ORDER BY potential_score DESC, heat_index DESC';

  if (filters?.limit) {
    sql += ' LIMIT ?';
    params.push(filters.limit);
  }

  return db.prepare<RecommendedTopic>(sql).all(...params);
}

/** 更新选题状态 */
export function updateStatus(db: CreatorMediaDB, id: string, status: TopicRecommendStatus): boolean {
  const now = new Date().toISOString();
  const result = db.prepare('UPDATE recommended_topics SET status = ?, updated_at = ? WHERE id = ?').run(status, now, id);
  return result.changes > 0;
}

/** 批量更新选题状态 */
export function batchUpdateStatus(db: CreatorMediaDB, ids: string[], status: TopicRecommendStatus): number {
  const now = new Date().toISOString();
  const placeholders = ids.map(() => '?').join(', ');
  const result = db.prepare(`UPDATE recommended_topics SET status = ?, updated_at = ? WHERE id IN (${placeholders})`).run(status, now, ...ids);
  return result.changes;
}

/** 关联内容（采纳后） */
export function linkContent(db: CreatorMediaDB, topicId: string, contentId: string): boolean {
  const now = new Date().toISOString();
  const result = db.prepare('UPDATE recommended_topics SET content_id = ?, status = 1, updated_at = ? WHERE id = ?').run(contentId, now, topicId);
  return result.changes > 0;
}

/** 删除选题 */
export function remove(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM recommended_topics WHERE id = ?').run(id);
  return result.changes > 0;
}
