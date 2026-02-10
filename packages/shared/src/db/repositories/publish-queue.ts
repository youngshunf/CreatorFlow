/**
 * 发布队列 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { PublishQueueItem, CreatePublishQueueInput, UpdatePublishQueueInput } from '../types.ts';

// 重新导出类型供外部使用
export type { PublishQueueItem, CreatePublishQueueInput, UpdatePublishQueueInput };

/** 入队 */
export function enqueue(db: CreatorMediaDB, data: CreatePublishQueueInput): PublishQueueItem {
  const stmt = db.prepare(`
    INSERT INTO publish_queue (id, content_id, platform_account_id, priority, status, scheduled_at, started_at, completed_at, error_message, retry_count, max_retries)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.content_id, data.platform_account_id,
    data.priority ?? 0, data.status ?? 'queued', data.scheduled_at,
    data.started_at, data.completed_at, data.error_message,
    data.retry_count ?? 0, data.max_retries ?? 3,
  );
  return getQueueItem(db, data.id)!;
}

/** 获取单个队列项 */
export function getQueueItem(db: CreatorMediaDB, id: string): PublishQueueItem | undefined {
  return db.prepare<PublishQueueItem>('SELECT * FROM publish_queue WHERE id = ?').get(id);
}

/** 按内容查询 */
export function listByContent(db: CreatorMediaDB, contentId: string): PublishQueueItem[] {
  return db.prepare<PublishQueueItem>(
    'SELECT * FROM publish_queue WHERE content_id = ? ORDER BY priority DESC, created_at ASC'
  ).all(contentId);
}

/** 按平台账号查询，可选按状态过滤 */
export function listByPlatformAccount(db: CreatorMediaDB, platformAccountId: string, status?: string): PublishQueueItem[] {
  if (status) {
    return db.prepare<PublishQueueItem>(
      'SELECT * FROM publish_queue WHERE platform_account_id = ? AND status = ? ORDER BY priority DESC, scheduled_at ASC'
    ).all(platformAccountId, status);
  }
  return db.prepare<PublishQueueItem>(
    'SELECT * FROM publish_queue WHERE platform_account_id = ? ORDER BY priority DESC, scheduled_at ASC'
  ).all(platformAccountId);
}

/** 获取下一个待处理项（status='queued', 按 priority DESC, scheduled_at ASC） */
export function getNextInQueue(db: CreatorMediaDB, platformAccountId: string): PublishQueueItem | undefined {
  return db.prepare<PublishQueueItem>(
    'SELECT * FROM publish_queue WHERE platform_account_id = ? AND status = \'queued\' AND (scheduled_at IS NULL OR scheduled_at <= CURRENT_TIMESTAMP) ORDER BY priority DESC, scheduled_at ASC LIMIT 1'
  ).get(platformAccountId);
}

/** 标记为处理中（使用 IMMEDIATE 事务防止并发冲突） */
export function markProcessing(db: CreatorMediaDB, id: string): PublishQueueItem | undefined {
  db.transaction(() => {
    db.prepare(
      'UPDATE publish_queue SET status = \'processing\', started_at = CURRENT_TIMESTAMP, updated_at = CURRENT_TIMESTAMP WHERE id = ? AND status = \'queued\''
    ).run(id);
  });
  return getQueueItem(db, id);
}

/** 标记完成 */
export function markCompleted(db: CreatorMediaDB, id: string): PublishQueueItem | undefined {
  db.prepare(
    'UPDATE publish_queue SET status = \'completed\', completed_at = CURRENT_TIMESTAMP, updated_at = CURRENT_TIMESTAMP WHERE id = ?'
  ).run(id);
  return getQueueItem(db, id);
}

/** 标记失败，retry_count++，如果 retry_count < max_retries 则重新入队 */
export function markFailed(db: CreatorMediaDB, id: string, errorMessage: string): PublishQueueItem | undefined {
  db.transaction(() => {
    db.prepare(
      'UPDATE publish_queue SET error_message = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?'
    ).run(errorMessage, id);

    // 如果未超过最大重试次数，重新入队
    db.prepare(
      'UPDATE publish_queue SET status = \'queued\', started_at = NULL WHERE id = ? AND retry_count < max_retries'
    ).run(id);

    // 如果已超过最大重试次数，标记为失败
    db.prepare(
      'UPDATE publish_queue SET status = \'failed\' WHERE id = ? AND retry_count >= max_retries'
    ).run(id);
  });
  return getQueueItem(db, id);
}

/** 取消某内容的所有待处理项 */
export function cancelByContent(db: CreatorMediaDB, contentId: string): number {
  const result = db.prepare(
    'UPDATE publish_queue SET status = \'cancelled\', updated_at = CURRENT_TIMESTAMP WHERE content_id = ? AND status IN (\'queued\', \'processing\')'
  ).run(contentId);
  return result.changes;
}

/** 检查某平台账号是否有正在处理的任务 */
export function isProcessing(db: CreatorMediaDB, platformAccountId: string): boolean {
  const row = db.prepare<{ cnt: number }>(
    'SELECT COUNT(*) as cnt FROM publish_queue WHERE platform_account_id = ? AND status = \'processing\''
  ).get(platformAccountId);
  return (row?.cnt ?? 0) > 0;
}

/** 获取某平台账号最后完成时间（用于 60 秒间隔检查） */
export function getLastCompletedTime(db: CreatorMediaDB, platformAccountId: string): string | null {
  const row = db.prepare<{ completed_at: string | null }>(
    'SELECT MAX(completed_at) as completed_at FROM publish_queue WHERE platform_account_id = ? AND status = \'completed\''
  ).get(platformAccountId);
  return row?.completed_at ?? null;
}

/** 更新队列项 */
export function updateQueueItem(db: CreatorMediaDB, id: string, data: UpdatePublishQueueInput): PublishQueueItem | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getQueueItem(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE publish_queue SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getQueueItem(db, id);
}
