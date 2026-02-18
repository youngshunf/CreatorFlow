/**
 * 视频项目 Repository — CRUD 操作
 */

import type { CreatorMediaDB } from "../connection.ts";
import type {
  VideoProject,
  VideoProjectFull,
  CreateVideoProject,
  UpdateVideoProject,
  VideoScene,
  CreateVideoScene,
  UpdateVideoScene,
  VideoAsset,
  CreateVideoAsset,
} from "../types.ts";

// ============================================================
// 项目 CRUD
// ============================================================

/** 创建视频项目 */
export function createVideoProject(
  db: CreatorMediaDB,
  data: CreateVideoProject,
): VideoProject {
  db.prepare(`
    INSERT INTO video_projects (id, content_id, name, description, width, height, fps, metadata)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
  `).run(
    data.id,
    data.content_id,
    data.name,
    data.description,
    data.width,
    data.height,
    data.fps,
    data.metadata,
  );
  return getVideoProject(db, data.id)!;
}

/** 获取单个视频项目 */
export function getVideoProject(
  db: CreatorMediaDB,
  id: string,
): VideoProject | undefined {
  return db
    .prepare<VideoProject>("SELECT * FROM video_projects WHERE id = ?")
    .get(id);
}

/** 获取完整视频项目（含场景和素材） */
export function getVideoProjectFull(
  db: CreatorMediaDB,
  id: string,
): VideoProjectFull | undefined {
  const project = getVideoProject(db, id);
  if (!project) return undefined;
  return {
    ...project,
    scenes: listVideoScenes(db, id),
    assets: listVideoAssets(db, id),
  };
}

/** 通过 content_id 获取视频项目 */
export function getVideoProjectByContentId(
  db: CreatorMediaDB,
  contentId: string,
): VideoProject | undefined {
  return db
    .prepare<VideoProject>(
      "SELECT * FROM video_projects WHERE content_id = ?",
    )
    .get(contentId);
}

/** 列出项目下的所有视频项目（通过 contents.project_id 关联） */
export function listVideoProjects(
  db: CreatorMediaDB,
  projectId: string,
): VideoProject[] {
  return db
    .prepare<VideoProject>(
      `SELECT vp.* FROM video_projects vp
       INNER JOIN contents c ON c.id = vp.content_id
       WHERE c.project_id = ?
       ORDER BY vp.updated_at DESC`,
    )
    .all(projectId);
}

/** 更新视频项目 */
export function updateVideoProject(
  db: CreatorMediaDB,
  id: string,
  data: UpdateVideoProject,
): VideoProject | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getVideoProject(db, id);

  if (!data.updated_at) {
    fields.push("updated_at = CURRENT_TIMESTAMP");
  }

  values.push(id);
  db.prepare(
    `UPDATE video_projects SET ${fields.join(", ")} WHERE id = ?`,
  ).run(...values);
  return getVideoProject(db, id);
}

/** 删除视频项目 */
export function deleteVideoProject(
  db: CreatorMediaDB,
  id: string,
): boolean {
  const result = db
    .prepare("DELETE FROM video_projects WHERE id = ?")
    .run(id);
  return result.changes > 0;
}

// ============================================================
// 场景 CRUD
// ============================================================

/** 添加视频场景 */
export function addVideoScene(
  db: CreatorMediaDB,
  data: CreateVideoScene,
): VideoScene {
  db.prepare(`
    INSERT INTO video_scenes (id, project_id, composition_id, name, sort_order, duration_in_frames, props, transition_type, transition_duration, transition_direction)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  `).run(
    data.id,
    data.project_id,
    data.composition_id,
    data.name,
    data.sort_order,
    data.duration_in_frames,
    data.props,
    data.transition_type,
    data.transition_duration,
    data.transition_direction,
  );
  return getVideoScene(db, data.id)!;
}

/** 获取单个场景 */
export function getVideoScene(
  db: CreatorMediaDB,
  id: string,
): VideoScene | undefined {
  return db
    .prepare<VideoScene>("SELECT * FROM video_scenes WHERE id = ?")
    .get(id);
}

/** 列出项目的所有场景（按 sort_order 排序） */
export function listVideoScenes(
  db: CreatorMediaDB,
  projectId: string,
): VideoScene[] {
  return db
    .prepare<VideoScene>(
      "SELECT * FROM video_scenes WHERE project_id = ? ORDER BY sort_order ASC",
    )
    .all(projectId);
}

/** 更新场景 */
export function updateVideoScene(
  db: CreatorMediaDB,
  id: string,
  data: UpdateVideoScene,
): VideoScene | undefined {
  const fields: string[] = [];
  const values: unknown[] = [];

  for (const [key, value] of Object.entries(data)) {
    if (value !== undefined) {
      fields.push(`${key} = ?`);
      values.push(value);
    }
  }

  if (fields.length === 0) return getVideoScene(db, id);

  if (!data.updated_at) {
    fields.push("updated_at = CURRENT_TIMESTAMP");
  }

  values.push(id);
  db.prepare(
    `UPDATE video_scenes SET ${fields.join(", ")} WHERE id = ?`,
  ).run(...values);
  return getVideoScene(db, id);
}

/** 删除场景 */
export function removeVideoScene(
  db: CreatorMediaDB,
  id: string,
): boolean {
  const result = db
    .prepare("DELETE FROM video_scenes WHERE id = ?")
    .run(id);
  return result.changes > 0;
}

/** 重新排序场景 */
export function reorderVideoScenes(
  db: CreatorMediaDB,
  projectId: string,
  sceneIds: string[],
): void {
  db.transaction(() => {
    const stmt = db.prepare(
      "UPDATE video_scenes SET sort_order = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ? AND project_id = ?",
    );
    for (let i = 0; i < sceneIds.length; i++) {
      stmt.run(i, sceneIds[i], projectId);
    }
  });
}

/** 获取下一个排序序号 */
export function getNextSortOrder(
  db: CreatorMediaDB,
  projectId: string,
): number {
  const row = db
    .prepare<{ max_order: number | null }>(
      "SELECT MAX(sort_order) as max_order FROM video_scenes WHERE project_id = ?",
    )
    .get(projectId);
  return (row?.max_order ?? -1) + 1;
}

// ============================================================
// 素材 CRUD
// ============================================================

/** 添加视频素材 */
export function addVideoAsset(
  db: CreatorMediaDB,
  data: CreateVideoAsset,
): VideoAsset {
  db.prepare(`
    INSERT INTO video_assets (id, project_id, type, name, file_path, file_size, metadata)
    VALUES (?, ?, ?, ?, ?, ?, ?)
  `).run(
    data.id,
    data.project_id,
    data.type,
    data.name,
    data.file_path,
    data.file_size,
    data.metadata,
  );
  return getVideoAsset(db, data.id)!;
}

/** 获取单个素材 */
export function getVideoAsset(
  db: CreatorMediaDB,
  id: string,
): VideoAsset | undefined {
  return db
    .prepare<VideoAsset>("SELECT * FROM video_assets WHERE id = ?")
    .get(id);
}

/** 列出项目的所有素材 */
export function listVideoAssets(
  db: CreatorMediaDB,
  projectId: string,
): VideoAsset[] {
  return db
    .prepare<VideoAsset>(
      "SELECT * FROM video_assets WHERE project_id = ? ORDER BY created_at DESC",
    )
    .all(projectId);
}

/** 删除素材 */
export function removeVideoAsset(
  db: CreatorMediaDB,
  id: string,
): boolean {
  const result = db
    .prepare("DELETE FROM video_assets WHERE id = ?")
    .run(id);
  return result.changes > 0;
}
