/**
 * 项目 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { Project, CreateProject, UpdateProject } from '../types.ts';

/** 创建项目 */
export function createProject(db: CreatorMediaDB, data: CreateProject): Project {
  const stmt = db.prepare<Project>(`
    INSERT INTO projects (id, name, description, platform, platforms, avatar_path, is_active)
    VALUES (?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.name,
    data.description,
    data.platform,
    data.platforms,
    data.avatar_path,
    data.is_active ?? 0,
  );
  return getProject(db, data.id)!;
}

/** 获取单个项目 */
export function getProject(db: CreatorMediaDB, id: string): Project | undefined {
  return db.prepare<Project>('SELECT * FROM projects WHERE id = ?').get(id);
}

/** 列出所有项目 */
export function listProjects(db: CreatorMediaDB): Project[] {
  return db.prepare<Project>('SELECT * FROM projects ORDER BY created_at DESC').all();
}

/** 更新项目 */
export function updateProject(db: CreatorMediaDB, id: string, data: UpdateProject): Project | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getProject(db, id);

  // 自动更新 updated_at
  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE projects SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getProject(db, id);
}

/** 删除项目 */
export function deleteProject(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM projects WHERE id = ?').run(id);
  return result.changes > 0;
}

/** 获取当前激活项目 */
export function getActiveProject(db: CreatorMediaDB): Project | undefined {
  return db.prepare<Project>('SELECT * FROM projects WHERE is_active = 1').get();
}

/** 设置激活项目（事务：先清除旧的，再设置新的） */
export function setActiveProject(db: CreatorMediaDB, projectId: string): Project | undefined {
  db.transaction(() => {
    db.prepare('UPDATE projects SET is_active = 0, updated_at = CURRENT_TIMESTAMP WHERE is_active = 1').run();
    db.prepare('UPDATE projects SET is_active = 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?').run(projectId);
  });
  return getProject(db, projectId);
}

/**
 * 获取项目的所有目标平台（合并 platform 和 platforms，去重）
 * @returns JSON 字符串数组，如 '["xiaohongshu","douyin"]'，如果项目不存在返回 null
 */
export function getProjectTargetPlatforms(db: CreatorMediaDB, projectId: string): string | null {
  const project = getProject(db, projectId);
  if (!project) return null;

  const platforms: string[] = [project.platform];

  if (project.platforms) {
    try {
      const parsed = JSON.parse(project.platforms);
      if (Array.isArray(parsed)) {
        platforms.push(...parsed);
      }
    } catch {
      // 忽略解析错误
    }
  }

  // 去重
  const uniquePlatforms = [...new Set(platforms)];
  return JSON.stringify(uniquePlatforms);
}
