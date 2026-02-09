/**
 * 热点缓存 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { TopicCache, CreateTopicCache, TopicCacheStatus } from '../types.ts';

/** 热点列表过滤条件 */
export interface TopicCacheFilters {
  source?: string;
  source_id?: string;
  status?: TopicCacheStatus | TopicCacheStatus[];
  minRelevance?: number;
}

/** 批量缓存热点 */
export function cacheTopics(db: CreatorMediaDB, topics: CreateTopicCache[]): void {
  const stmt = db.prepare(`
    INSERT OR REPLACE INTO topic_cache (id, source, source_id, title, url, heat_score, relevance_score, category, keywords, status, fetched_at, expires_at)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);

  db.transaction(() => {
    for (const t of topics) {
      stmt.run(
        t.id, t.source, t.source_id, t.title, t.url,
        t.heat_score, t.relevance_score, t.category,
        t.keywords, t.status ?? 'new', t.fetched_at, t.expires_at,
      );
    }
  });
}

/** 列出缓存的热点（支持过滤） */
export function listCachedTopics(db: CreatorMediaDB, filters?: TopicCacheFilters): TopicCache[] {
  let sql = 'SELECT * FROM topic_cache WHERE 1=1';
  const params: unknown[] = [];

  if (filters?.source) {
    sql += ' AND source = ?';
    params.push(filters.source);
  }

  if (filters?.source_id) {
    sql += ' AND source_id = ?';
    params.push(filters.source_id);
  }

  if (filters?.status) {
    if (Array.isArray(filters.status)) {
      const placeholders = filters.status.map(() => '?').join(', ');
      sql += ` AND status IN (${placeholders})`;
      params.push(...filters.status);
    } else {
      sql += ' AND status = ?';
      params.push(filters.status);
    }
  }

  if (filters?.minRelevance !== undefined) {
    sql += ' AND relevance_score >= ?';
    params.push(filters.minRelevance);
  }

  sql += ' ORDER BY relevance_score DESC, heat_score DESC';
  return db.prepare<TopicCache>(sql).all(...params);
}

/** 更新热点状态 */
export function updateTopicStatus(db: CreatorMediaDB, id: string, status: TopicCacheStatus): boolean {
  const result = db.prepare('UPDATE topic_cache SET status = ? WHERE id = ?').run(status, id);
  return result.changes > 0;
}

/** 清理过期热点 */
export function cleanExpiredTopics(db: CreatorMediaDB): number {
  const result = db.prepare(
    "UPDATE topic_cache SET status = 'expired' WHERE expires_at IS NOT NULL AND expires_at < CURRENT_TIMESTAMP AND status != 'expired'"
  ).run();
  return result.changes;
}
