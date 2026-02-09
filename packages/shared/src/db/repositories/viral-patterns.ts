/**
 * 爆款模式库 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { ViralPattern, CreateViralPattern, UpdateViralPattern, ViralPatternCategory } from '../types.ts';

/** 爆款模式过滤条件 */
export interface ViralPatternFilters {
  category?: ViralPatternCategory;
  platform?: string;
  projectId?: string | null; // null = 仅全局模式
}

/** 列出爆款模式（支持过滤） */
export function listViralPatterns(db: CreatorMediaDB, filters?: ViralPatternFilters): ViralPattern[] {
  let sql = 'SELECT * FROM viral_patterns WHERE 1=1';
  const params: unknown[] = [];

  if (filters?.category) {
    sql += ' AND category = ?';
    params.push(filters.category);
  }

  if (filters?.platform) {
    sql += ' AND (platform = ? OR platform IS NULL)';
    params.push(filters.platform);
  }

  if (filters?.projectId === null) {
    // 仅全局模式
    sql += ' AND project_id IS NULL';
  } else if (filters?.projectId) {
    // 项目级 + 全局
    sql += ' AND (project_id = ? OR project_id IS NULL)';
    params.push(filters.projectId);
  }

  sql += ' ORDER BY usage_count DESC, success_rate DESC NULLS LAST';
  return db.prepare<ViralPattern>(sql).all(...params);
}

/** 获取单个爆款模式 */
export function getViralPattern(db: CreatorMediaDB, id: string): ViralPattern | undefined {
  return db.prepare<ViralPattern>('SELECT * FROM viral_patterns WHERE id = ?').get(id);
}

/** 创建爆款模式 */
export function createViralPattern(db: CreatorMediaDB, data: CreateViralPattern): ViralPattern {
  const stmt = db.prepare(`
    INSERT INTO viral_patterns (id, project_id, platform, category, name, description, template, examples, source, usage_count, success_rate, tags)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.project_id, data.platform, data.category,
    data.name, data.description, data.template, data.examples,
    data.source, data.usage_count ?? 0, data.success_rate, data.tags,
  );
  return getViralPattern(db, data.id)!;
}

/** 更新爆款模式 */
export function updateViralPattern(db: CreatorMediaDB, id: string, data: UpdateViralPattern): ViralPattern | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getViralPattern(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE viral_patterns SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getViralPattern(db, id);
}

/** 增加使用次数 */
export function incrementUsageCount(db: CreatorMediaDB, id: string): ViralPattern | undefined {
  db.prepare('UPDATE viral_patterns SET usage_count = usage_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?').run(id);
  return getViralPattern(db, id);
}

/** 更新成功率 */
export function updateSuccessRate(db: CreatorMediaDB, id: string, rate: number): ViralPattern | undefined {
  db.prepare('UPDATE viral_patterns SET success_rate = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?').run(rate, id);
  return getViralPattern(db, id);
}

/** 删除爆款模式 */
export function deleteViralPattern(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM viral_patterns WHERE id = ?').run(id);
  return result.changes > 0;
}
