/**
 * 平台账号 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { PlatformAccount, CreatePlatformAccount, UpdatePlatformAccount } from '../types.ts';

/** 创建平台账号 */
export function createPlatformAccount(db: CreatorMediaDB, data: CreatePlatformAccount): PlatformAccount {
  const stmt = db.prepare(`
    INSERT INTO platform_accounts (id, project_id, platform, platform_uid, nickname, avatar_url, bio, home_url, followers, following, total_likes, total_favorites, total_comments, total_posts, metrics_json, metrics_updated_at, auth_status, auth_method, auth_data, auth_expires_at, last_login_at, last_login_check, is_primary, notes)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.project_id, data.platform, data.platform_uid,
    data.nickname, data.avatar_url, data.bio, data.home_url,
    data.followers ?? 0, data.following ?? 0, data.total_likes ?? 0,
    data.total_favorites ?? 0, data.total_comments ?? 0, data.total_posts ?? 0,
    data.metrics_json, data.metrics_updated_at,
    data.auth_status ?? 'not_logged_in', data.auth_method, data.auth_data,
    data.auth_expires_at, data.last_login_at, data.last_login_check,
    data.is_primary ?? 0, data.notes,
  );
  return db.prepare<PlatformAccount>('SELECT * FROM platform_accounts WHERE id = ?').get(data.id)!;
}

/** 列出项目的所有平台账号 */
export function listByProject(db: CreatorMediaDB, projectId: string): PlatformAccount[] {
  return db.prepare<PlatformAccount>(
    'SELECT * FROM platform_accounts WHERE project_id = ? ORDER BY is_primary DESC, created_at ASC'
  ).all(projectId);
}

/** 更新平台账号 */
export function updatePlatformAccount(db: CreatorMediaDB, id: string, data: UpdatePlatformAccount): PlatformAccount | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) {
    return db.prepare<PlatformAccount>('SELECT * FROM platform_accounts WHERE id = ?').get(id);
  }

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE platform_accounts SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return db.prepare<PlatformAccount>('SELECT * FROM platform_accounts WHERE id = ?').get(id);
}

/** 删除平台账号 */
export function deletePlatformAccount(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM platform_accounts WHERE id = ?').run(id);
  return result.changes > 0;
}
