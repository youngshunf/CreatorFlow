/**
 * 采集调度任务 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { ReviewTask, CreateReviewTask, UpdateReviewTask } from '../types.ts';

// 重新导出类型供外部使用
export type { ReviewTask, CreateReviewTask, UpdateReviewTask };

/** 创建采集任务 */
export function createReviewTask(db: CreatorMediaDB, data: CreateReviewTask): ReviewTask {
  const stmt = db.prepare(`
    INSERT INTO review_tasks (id, publish_record_id, scheduled_at, executed_at, status, review_type, result_snapshot, error_message, retry_count)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.publish_record_id, data.scheduled_at,
    data.executed_at, data.status ?? 'pending', data.review_type,
    data.result_snapshot, data.error_message, data.retry_count ?? 0,
  );
  return getReviewTask(db, data.id)!;
}

/** 获取单个采集任务 */
export function getReviewTask(db: CreatorMediaDB, id: string): ReviewTask | undefined {
  return db.prepare<ReviewTask>('SELECT * FROM review_tasks WHERE id = ?').get(id);
}

/** 按发布记录查询采集任务 */
export function listByPublishRecord(db: CreatorMediaDB, publishRecordId: string): ReviewTask[] {
  return db.prepare<ReviewTask>(
    'SELECT * FROM review_tasks WHERE publish_record_id = ? ORDER BY scheduled_at ASC'
  ).all(publishRecordId);
}

/** 查询待执行任务（status='pending', scheduled_at <= now） */
export function listPendingTasks(db: CreatorMediaDB, limit?: number): ReviewTask[] {
  let sql = 'SELECT * FROM review_tasks WHERE status = \'pending\' AND scheduled_at <= CURRENT_TIMESTAMP ORDER BY scheduled_at ASC';
  if (limit) {
    sql += ` LIMIT ${limit}`;
  }
  return db.prepare<ReviewTask>(sql).all();
}

/** 更新采集任务 */
export function updateReviewTask(db: CreatorMediaDB, id: string, data: UpdateReviewTask): ReviewTask | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getReviewTask(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE review_tasks SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getReviewTask(db, id);
}

/** 标记为执行中 */
export function markExecuting(db: CreatorMediaDB, id: string): ReviewTask | undefined {
  db.prepare(
    'UPDATE review_tasks SET status = \'executing\', updated_at = CURRENT_TIMESTAMP WHERE id = ?'
  ).run(id);
  return getReviewTask(db, id);
}

/** 标记完成并保存结果 */
export function markCompleted(db: CreatorMediaDB, id: string, resultSnapshot: string): ReviewTask | undefined {
  db.prepare(
    'UPDATE review_tasks SET status = \'completed\', executed_at = CURRENT_TIMESTAMP, result_snapshot = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?'
  ).run(resultSnapshot, id);
  return getReviewTask(db, id);
}

/** 标记失败，retry_count++ */
export function markFailed(db: CreatorMediaDB, id: string, errorMessage: string): ReviewTask | undefined {
  db.prepare(
    'UPDATE review_tasks SET status = \'failed\', executed_at = CURRENT_TIMESTAMP, error_message = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?'
  ).run(errorMessage, id);
  return getReviewTask(db, id);
}

/** 取消某发布记录的所有待执行任务 */
export function cancelByPublishRecord(db: CreatorMediaDB, publishRecordId: string): number {
  const result = db.prepare(
    'UPDATE review_tasks SET status = \'cancelled\', updated_at = CURRENT_TIMESTAMP WHERE publish_record_id = ? AND status = \'pending\''
  ).run(publishRecordId);
  return result.changes;
}
