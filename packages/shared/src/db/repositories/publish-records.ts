/**
 * 发布记录 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { PublishRecord, CreatePublishRecord, UpdatePublishRecord, PublishStatus } from '../types.ts';

// 重新导出类型供外部使用
export type { PublishRecord, CreatePublishRecord, UpdatePublishRecord };

/** 发布记录列表过滤条件 */
export interface PublishRecordFilters {
  status?: PublishStatus | PublishStatus[];
  platform?: string;
}

/** 创建发布记录 */
export function createPublishRecord(db: CreatorMediaDB, data: CreatePublishRecord): PublishRecord {
  const stmt = db.prepare(`
    INSERT INTO publish_records (id, content_id, platform_account_id, platform, publish_url, status, method, error_message, published_at, views, likes, comments, shares, favorites, metrics_json, metrics_updated_at)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.content_id, data.platform_account_id, data.platform,
    data.publish_url, data.status ?? 'pending', data.method,
    data.error_message, data.published_at,
    data.views ?? 0, data.likes ?? 0, data.comments ?? 0,
    data.shares ?? 0, data.favorites ?? 0,
    data.metrics_json, data.metrics_updated_at,
  );
  return getPublishRecord(db, data.id)!;
}

/** 获取单个发布记录 */
export function getPublishRecord(db: CreatorMediaDB, id: string): PublishRecord | undefined {
  return db.prepare<PublishRecord>('SELECT * FROM publish_records WHERE id = ?').get(id);
}

/** 列出某内容的发布记录 */
export function listByContent(db: CreatorMediaDB, contentId: string): PublishRecord[] {
  return db.prepare<PublishRecord>('SELECT * FROM publish_records WHERE content_id = ? ORDER BY created_at DESC').all(contentId);
}

/** 列出项目的发布记录（通过 JOIN contents 表） */
export function listByProject(db: CreatorMediaDB, projectId: string, filters?: PublishRecordFilters): PublishRecord[] {
  let sql = `
    SELECT pr.* FROM publish_records pr
    INNER JOIN contents c ON c.id = pr.content_id
    WHERE c.project_id = ?
  `;
  const params: unknown[] = [projectId];

  if (filters?.status) {
    if (Array.isArray(filters.status)) {
      const placeholders = filters.status.map(() => '?').join(', ');
      sql += ` AND pr.status IN (${placeholders})`;
      params.push(...filters.status);
    } else {
      sql += ' AND pr.status = ?';
      params.push(filters.status);
    }
  }

  if (filters?.platform) {
    sql += ' AND pr.platform = ?';
    params.push(filters.platform);
  }

  sql += ' ORDER BY pr.created_at DESC';
  return db.prepare<PublishRecord>(sql).all(...params);
}

/** 更新发布记录 */
export function updatePublishRecord(db: CreatorMediaDB, id: string, data: UpdatePublishRecord): PublishRecord | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getPublishRecord(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE publish_records SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getPublishRecord(db, id);
}

/** 更新发布指标（快捷方法） */
export function updatePublishMetrics(db: CreatorMediaDB, id: string, metrics: {
  views?: number;
  likes?: number;
  comments?: number;
  shares?: number;
  favorites?: number;
  metrics_json?: string | null;
}): PublishRecord | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  if (metrics.views !== undefined) { fields.push('views = ?'); values.push(metrics.views); }
  if (metrics.likes !== undefined) { fields.push('likes = ?'); values.push(metrics.likes); }
  if (metrics.comments !== undefined) { fields.push('comments = ?'); values.push(metrics.comments); }
  if (metrics.shares !== undefined) { fields.push('shares = ?'); values.push(metrics.shares); }
  if (metrics.favorites !== undefined) { fields.push('favorites = ?'); values.push(metrics.favorites); }
  if (metrics.metrics_json !== undefined) { fields.push('metrics_json = ?'); values.push(metrics.metrics_json); }

  if (fields.length === 0) return getPublishRecord(db, id);

  fields.push('metrics_updated_at = CURRENT_TIMESTAMP');
  fields.push('updated_at = CURRENT_TIMESTAMP');

  values.push(id);
  db.prepare(`UPDATE publish_records SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getPublishRecord(db, id);
}

/** 删除发布记录 */
export function deletePublishRecord(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM publish_records WHERE id = ?').run(id);
  return result.changes > 0;
}
