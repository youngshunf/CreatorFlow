/**
 * 内容创作 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from '../connection.ts';
import type { Content, CreateContent, UpdateContent, ContentStatus, ContentVideoMetadata } from '../types.ts';

/** 内容列表过滤条件 */
export interface ContentFilters {
  status?: ContentStatus | ContentStatus[];
  content_type?: string;
  pipeline_mode?: string;
}

/** 创建内容 */
export function createContent(db: CreatorMediaDB, data: CreateContent): Content {
  const stmt = db.prepare(`
    INSERT INTO contents (id, project_id, title, topic, topic_source, source_topic_id, script_path, status, content_type, target_platforms, pipeline_mode, pipeline_state, viral_pattern_id, tags, scheduled_at, files, metadata)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `);
  stmt.run(
    data.id, data.project_id, data.title, data.topic,
    data.topic_source, data.source_topic_id, data.script_path,
    data.status ?? 'idea', data.content_type, data.target_platforms,
    data.pipeline_mode ?? 'semi-auto', data.pipeline_state,
    data.viral_pattern_id, data.tags, data.scheduled_at,
    data.files, data.metadata,
  );
  return getContent(db, data.id)!;
}

/** 获取单个内容 */
export function getContent(db: CreatorMediaDB, id: string): Content | undefined {
  return db.prepare<Content>('SELECT * FROM contents WHERE id = ?').get(id);
}

/** 列出项目的内容（支持过滤） */
export function listByProject(db: CreatorMediaDB, projectId: string, filters?: ContentFilters): Content[] {
  let sql = 'SELECT * FROM contents WHERE project_id = ?';
  const params: unknown[] = [projectId];

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

  if (filters?.content_type) {
    sql += ' AND content_type = ?';
    params.push(filters.content_type);
  }

  if (filters?.pipeline_mode) {
    sql += ' AND pipeline_mode = ?';
    params.push(filters.pipeline_mode);
  }

  sql += ' ORDER BY updated_at DESC';
  return db.prepare<Content>(sql).all(...params);
}

/** 更新内容 */
export function updateContent(db: CreatorMediaDB, id: string, data: UpdateContent): Content | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getContent(db, id);

  if (!data.updated_at) {
    fields.push('updated_at = CURRENT_TIMESTAMP');
  }

  values.push(id);
  db.prepare(`UPDATE contents SET ${fields.join(', ')} WHERE id = ?`).run(...values);
  return getContent(db, id);
}

/** 更新内容状态（快捷方法） */
export function updateContentStatus(db: CreatorMediaDB, id: string, status: ContentStatus): Content | undefined {
  db.prepare('UPDATE contents SET status = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?').run(status, id);
  return getContent(db, id);
}

/** 删除内容 */
export function deleteContent(db: CreatorMediaDB, id: string): boolean {
  const result = db.prepare('DELETE FROM contents WHERE id = ?').run(id);
  return result.changes > 0;
}

// ============================================================
// 视频内容辅助函数
// ============================================================

/** 解析内容的 metadata JSON */
export function parseContentMetadata(content: Content): Record<string, unknown> | null {
  if (!content.metadata) return null;
  try {
    return JSON.parse(content.metadata) as Record<string, unknown>;
  } catch {
    return null;
  }
}

/** 获取内容的视频元数据 */
export function getContentVideoMetadata(content: Content): ContentVideoMetadata | null {
  const metadata = parseContentMetadata(content);
  if (!metadata) return null;
  return {
    video_project_id: metadata.video_project_id as string | undefined,
    video_project_name: metadata.video_project_name as string | undefined,
    video_template_id: metadata.video_template_id as string | undefined,
    video_render_status: metadata.video_render_status as ContentVideoMetadata['video_render_status'],
    video_output_path: metadata.video_output_path as string | undefined,
    video_duration: metadata.video_duration as number | undefined,
    video_resolution: metadata.video_resolution as ContentVideoMetadata['video_resolution'],
  };
}

/** 合并更新内容的视频元数据 */
export function updateContentVideoMetadata(
  db: CreatorMediaDB,
  contentId: string,
  videoMeta: Partial<ContentVideoMetadata>
): Content | undefined {
  const content = getContent(db, contentId);
  if (!content) return undefined;

  const existingMeta = parseContentMetadata(content) || {};
  const mergedMeta = { ...existingMeta, ...videoMeta };
  const metadataJson = JSON.stringify(mergedMeta);

  db.prepare('UPDATE contents SET metadata = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?').run(metadataJson, contentId);
  return getContent(db, contentId);
}

/** 判断是否为视频类型内容 */
export function isVideoContent(content: Content): boolean {
  return content.content_type === 'video' || content.content_type === 'short-video';
}
