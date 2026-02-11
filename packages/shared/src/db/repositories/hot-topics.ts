/**
 * 热榜快照 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { HotTopic, CreateHotTopic } from '../types.ts';

/** 热榜列表过滤条件 */
export interface HotTopicFilters {
  platformId?: string;
  batchDate?: string;
  limit?: number;
}

/** 批量插入热榜数据（去重） */
export function batchInsert(db: CreatorMediaDB, topics: CreateHotTopic[]): void {
  const stmt = db.prepare(`
    INSERT OR IGNORE INTO hot_topics (id, platform_id, platform_name, title, url, rank, heat_score, fetch_source, fetched_at, batch_date)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);

  db.transaction(() => {
    for (const t of topics) {
      stmt.run(
        t.id, t.platform_id, t.platform_name, t.title, t.url,
        t.rank, t.heat_score, t.fetch_source, t.fetched_at, t.batch_date,
      );
    }
  });
}

/** 按批次查询热榜 */
export function listByBatch(db: CreatorMediaDB, batchDate: string, platformId?: string): HotTopic[] {
  let sql = 'SELECT * FROM hot_topics WHERE batch_date = ?';
  const params: unknown[] = [batchDate];

  if (platformId) {
    sql += ' AND platform_id = ?';
    params.push(platformId);
  }

  sql += ' ORDER BY rank ASC, heat_score DESC';
  return db.prepare<HotTopic>(sql).all(...params);
}

/** 获取最新批次信息 */
export function getLatestBatch(db: CreatorMediaDB): { batchDate: string; fetchedAt: string; count: number } | undefined {
  const row = db.prepare<{ batch_date: string; fetched_at: string; cnt: number }>(
    'SELECT batch_date, MAX(fetched_at) as fetched_at, COUNT(*) as cnt FROM hot_topics GROUP BY batch_date ORDER BY fetched_at DESC LIMIT 1'
  ).get();
  if (!row) return undefined;
  return { batchDate: row.batch_date, fetchedAt: row.fetched_at, count: row.cnt };
}

/** 按过滤条件查询热榜 */
export function listByFilters(db: CreatorMediaDB, filters?: HotTopicFilters): HotTopic[] {
  let sql = 'SELECT * FROM hot_topics WHERE 1=1';
  const params: unknown[] = [];

  if (filters?.platformId) {
    sql += ' AND platform_id = ?';
    params.push(filters.platformId);
  }

  if (filters?.batchDate) {
    sql += ' AND batch_date = ?';
    params.push(filters.batchDate);
  }

  sql += ' ORDER BY rank ASC, heat_score DESC';

  if (filters?.limit) {
    sql += ' LIMIT ?';
    params.push(filters.limit);
  }

  return db.prepare<HotTopic>(sql).all(...params);
}

/** 清理旧数据 */
export function cleanOldData(db: CreatorMediaDB, keepDays: number): number {
  const result = db.prepare(
    "DELETE FROM hot_topics WHERE batch_date < date('now', '-' || ? || ' days')"
  ).run(keepDays);
  return result.changes;
}
