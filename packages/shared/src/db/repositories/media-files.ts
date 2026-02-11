/**
 * 素材库 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { MediaFile, CreateMediaFile, MediaFileType } from '../types.ts';

/** 素材列表过滤条件 */
export interface MediaFileFilters {
  type?: MediaFileType;
  tags?: string[];
}

/** 创建素材记录 */
export function createMediaFile(db: CreatorMediaDB, data: CreateMediaFile): MediaFile {
  const stmt = db.prepare(`
    INSERT INTO media_files (id, project_id, type, path, filename, size, width, height, duration, thumbnail, tags, description)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id,
    data.project_id,
    data.type,
    data.path,
    data.filename,
    data.size,
    data.width,
    data.height,
    data.duration,
    data.thumbnail,
    data.tags,
    data.description,
  );
  return getMediaFile(db, data.id)!;
}

/** 获取单个素材 */
export function getMediaFile(db: CreatorMediaDB, id: string): MediaFile | undefined {
  return db.prepare<MediaFile>('SELECT * FROM media_files WHERE id = ?').get(id);
}

/** 列出项目的素材（支持过滤） */
export function listByProject(db: CreatorMediaDB, projectId: string, filters?: MediaFileFilters): MediaFile[] {
  let sql = 'SELECT * FROM media_files WHERE project_id = ?';
  const params: unknown[] = [projectId];

  if (filters?.type) {
    sql += ' AND type = ?';
    params.push(filters.type);
  }

  if (filters?.tags && filters.tags.length > 0) {
    for (const tag of filters.tags) {
      sql += ' AND tags LIKE ?';
      params.push(`%${tag}%`);
    }
  }

  sql += ' ORDER BY created_at DESC';
  return db.prepare<MediaFile>(sql).all(...params);
}

/** 删除素材 */
export function deleteMediaFile(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM media_files WHERE id = ?').run(id);
  return result.changes > 0;
}
