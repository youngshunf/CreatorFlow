/**
 * 账号画像 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { AccountProfile, CreateAccountProfile, UpdateAccountProfile } from '../types.ts';

/** 创建画像 */
export function createProfile(db: CreatorMediaDB, data: CreateAccountProfile): AccountProfile {
  const stmt = db.prepare(`
    INSERT INTO account_profiles (id, project_id, niche, sub_niche, persona, target_audience, tone, keywords, bio, content_pillars, posting_frequency, best_posting_time, style_references, taboo_topics)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.project_id,
    data.niche,
    data.sub_niche,
    data.persona,
    data.target_audience,
    data.tone,
    data.keywords,
    data.bio,
    data.content_pillars,
    data.posting_frequency,
    data.best_posting_time,
    data.style_references,
    data.taboo_topics,
  );
  return getProfileByProject(db, data.project_id)!;
}

/** 获取项目的画像（1:1 关系） */
export function getProfileByProject(db: CreatorMediaDB, projectId: string): AccountProfile | undefined {
  return db.prepare<AccountProfile>(
    'SELECT * FROM account_profiles WHERE project_id = ?'
  ).get(projectId);
}

/** 更新画像 */
export function updateProfile(db: CreatorMediaDB, projectId: string, data: UpdateAccountProfile): AccountProfile | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getProfileByProject(db, projectId);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(projectId);
  db.prepare(`UPDATE account_profiles SET ${fields.join(', ')} WHERE project_id = ?`).run(...values);
  return getProfileByProject(db, projectId);
}

/** 创建或更新画像（upsert） */
export function upsertProfile(db: CreatorMediaDB, data: CreateAccountProfile): AccountProfile {
  const existing = getProfileByProject(db, data.project_id);
  if (existing) {
    const { id: _id, project_id: _pid, ...updateData } = data;
    return updateProfile(db, data.project_id, updateData)!;
  }
  return createProfile(db, data);
}

/** 删除画像 */
export function deleteProfile(db: CreatorMediaDB, projectId: string): boolean {
  const result = db.prepare('DELETE FROM account_profiles WHERE project_id = ?').run(projectId);
  return result.changes > 0;
}
