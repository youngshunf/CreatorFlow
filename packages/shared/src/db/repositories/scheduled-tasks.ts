/**
 * 定时任务 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { ScheduledTask, CreateScheduledTask, UpdateScheduledTask } from '../types.ts';

// 重新导出类型供外部使用
export type { ScheduledTask, CreateScheduledTask, UpdateScheduledTask };

/** 创建定时任务 */
export function createScheduledTask(db: CreatorMediaDB, data: CreateScheduledTask): ScheduledTask {
  const stmt = db.prepare(`
    INSERT INTO scheduled_tasks (id, project_id, name, description, task_type, schedule_mode, cron_expression, interval_seconds, scheduled_at, enabled, status, next_run_at, run_count, last_run_at, last_error, payload)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.project_id,
    data.name,
    data.description,
    data.task_type ?? 'custom',
    data.schedule_mode ?? 'cron',
    data.cron_expression,
    data.interval_seconds,
    data.scheduled_at,
    data.enabled ?? 1,
    data.status ?? 'active',
    data.next_run_at,
    data.run_count ?? 0,
    data.last_run_at,
    data.last_error,
    data.payload,
  );
  return getScheduledTask(db, data.id)!;
}

/** 获取单个定时任务 */
export function getScheduledTask(db: CreatorMediaDB, id: string): ScheduledTask | undefined {
  return db.prepare<ScheduledTask>('SELECT * FROM scheduled_tasks WHERE id = ?').get(id);
}

/** 列出定时任务（支持过滤） */
export function listScheduledTasks(
  db: CreatorMediaDB,
  filter?: { project_id?: string; task_type?: string; status?: string },
): ScheduledTask[] {
  const conditions: string[] = [];
  const params: unknown[] = [];

  if (filter?.project_id) {
    conditions.push('project_id = ?');
    params.push(filter.project_id);
  }
  if (filter?.task_type) {
    conditions.push('task_type = ?');
    params.push(filter.task_type);
  }
  if (filter?.status) {
    conditions.push('status = ?');
    params.push(filter.status);
  }

  const where = conditions.length > 0 ? ` WHERE ${conditions.join(' AND ')}` : '';
  return db.prepare<ScheduledTask>(
    `SELECT * FROM scheduled_tasks${where} ORDER BY created_at DESC`
  ).all(...params);
}

/** 更新定时任务 */
export function updateScheduledTask(db: CreatorMediaDB, id: string, data: UpdateScheduledTask): ScheduledTask | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getScheduledTask(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE scheduled_tasks SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getScheduledTask(db, id);
}

/** 删除定时任务 */
export function deleteScheduledTask(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM scheduled_tasks WHERE id = ?').run(id);
  return result.changes > 0;
}

/** 切换启用/禁用 */
export function toggleEnabled(db: CreatorMediaDB, id: string, enabled: boolean): ScheduledTask | undefined {
  db.prepare(
    'UPDATE scheduled_tasks SET enabled = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?'
  ).run(enabled ? 1 : 0, id);
  return getScheduledTask(db, id);
}

/** 标记运行完成：更新 last_run_at、run_count++、next_run_at */
export function markRunCompleted(db: CreatorMediaDB, id: string, nextRunAt?: string | null): ScheduledTask | undefined {
  const setNextRun = nextRunAt !== undefined ? ', next_run_at = ?' : '';
  const params: unknown[] = [id];
  if (nextRunAt !== undefined) params.unshift(nextRunAt);

  db.prepare(
    `UPDATE scheduled_tasks SET last_run_at = CURRENT_TIMESTAMP, run_count = run_count + 1, last_error = NULL, status = 'active', updated_at = CURRENT_TIMESTAMP${setNextRun} WHERE id = ?`
  ).run(...params);
  return getScheduledTask(db, id);
}

/** 标记运行失败：更新 last_error、status='error' */
export function markRunFailed(db: CreatorMediaDB, id: string, error: string): ScheduledTask | undefined {
  db.prepare(
    "UPDATE scheduled_tasks SET last_error = ?, status = 'error', last_run_at = CURRENT_TIMESTAMP, updated_at = CURRENT_TIMESTAMP WHERE id = ?"
  ).run(error, id);
  return getScheduledTask(db, id);
}

/** 查询所有已启用且到期的任务 */
export function listDueTasks(db: CreatorMediaDB, limit?: number): ScheduledTask[] {
  let sql = "SELECT * FROM scheduled_tasks WHERE enabled = 1 AND status = 'active' AND next_run_at <= CURRENT_TIMESTAMP ORDER BY next_run_at ASC";
  if (limit) sql += ` LIMIT ${limit}`;
  return db.prepare<ScheduledTask>(sql).all();
}
