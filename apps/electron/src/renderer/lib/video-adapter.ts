/**
 * 视频适配层 — VideoProjectFull (DB) ↔ VideoProject (旧类型) 转换
 *
 * 策略：在 VideoStudio / VideoEditor 层做一次转换，
 * 子组件（VideoPreview, VideoTimeline, VideoProperties 等）继续使用旧类型，零改动。
 */

import type {
  VideoProjectFull,
  VideoScene,
  VideoAsset,
  UpdateVideoProject,
  UpdateVideoScene,
} from '@sprouty-ai/shared/db/types'
import type {
  VideoProject,
  Scene,
  Transition,
  Asset,
  VideoConfig,
} from '@sprouty-ai/video'

// ============================================================================
// DB → Legacy 转换
// ============================================================================

/**
 * 安全解析 JSON 字符串，失败返回空对象
 */
function safeParseJson(json: string | null | undefined): Record<string, unknown> {
  if (!json) return {}
  try {
    return JSON.parse(json) as Record<string, unknown>
  } catch {
    return {}
  }
}

/**
 * 将 DB VideoScene 转为旧 Scene 类型
 */
export function toLegacyScene(s: VideoScene): Scene {
  return {
    id: s.id,
    name: s.name || '',
    compositionId: s.composition_id,
    durationInFrames: s.duration_in_frames,
    props: safeParseJson(s.props),
  }
}

/**
 * 将 DB VideoAsset 转为旧 Asset 类型
 */
export function toLegacyAsset(a: VideoAsset): Asset {
  return {
    id: a.id,
    type: a.type,
    name: a.name,
    path: a.file_path,
  }
}

/**
 * 从 scenes[1..n] 的 transition_* 字段提取旧 transitions[] 数组
 *
 * 语义偏移：
 * - 新模型：scene.transition_* 表示"进入该场景的过渡效果"
 * - 旧模型：transitions[i] 表示 scenes[i] 和 scenes[i+1] 之间的过渡
 * - 因此 transitions[i] = scenes[i+1].transition_*
 */
function extractTransitions(scenes: VideoScene[]): Transition[] {
  if (scenes.length <= 1) return []

  const transitions: Transition[] = []
  for (let i = 1; i < scenes.length; i++) {
    const s = scenes[i]
    transitions.push({
      type: s.transition_type === 'none' ? 'none' : s.transition_type,
      durationInFrames: s.transition_duration || 15,
      direction: s.transition_direction ?? undefined,
    })
  }
  return transitions
}

/**
 * 计算总帧数：场景帧数之和 - 过渡重叠帧数
 */
function calculateDurationInFrames(scenes: VideoScene[]): number {
  if (scenes.length === 0) return 0

  let total = 0
  for (const s of scenes) {
    total += s.duration_in_frames
  }
  // 减去过渡重叠（从第二个场景开始，每个过渡的 duration 是重叠帧数）
  for (let i = 1; i < scenes.length; i++) {
    if (scenes[i].transition_type !== 'none') {
      total -= scenes[i].transition_duration || 0
    }
  }
  return Math.max(total, 1)
}

/**
 * 主转换函数：VideoProjectFull (DB) → VideoProject (旧类型)
 */
export function toLegacyProject(full: VideoProjectFull): VideoProject {
  // 按 sort_order 排序场景
  const sortedScenes = [...full.scenes].sort((a, b) => a.sort_order - b.sort_order)

  const scenes = sortedScenes.map(toLegacyScene)
  const transitions = extractTransitions(sortedScenes)
  const assets = full.assets.map(toLegacyAsset)

  const config: VideoConfig = {
    width: full.width,
    height: full.height,
    fps: full.fps,
    durationInFrames: calculateDurationInFrames(sortedScenes),
  }

  return {
    id: full.id,
    name: full.name,
    description: full.description ?? undefined,
    createdAt: full.created_at,
    updatedAt: full.updated_at,
    config,
    scenes,
    transitions,
    assets,
    renders: [],
  }
}

// ============================================================================
// Legacy → DB 反向转换（用于更新操作）
// ============================================================================

/**
 * 旧 VideoProject 部分更新 → DB UpdateVideoProject
 */
export function toDbProjectUpdate(
  updates: Partial<Pick<VideoProject, 'name' | 'description' | 'config'>>,
): UpdateVideoProject {
  const result: UpdateVideoProject = {}

  if (updates.name !== undefined) result.name = updates.name
  if (updates.description !== undefined) result.description = updates.description ?? null
  if (updates.config) {
    if (updates.config.width !== undefined) result.width = updates.config.width
    if (updates.config.height !== undefined) result.height = updates.config.height
    if (updates.config.fps !== undefined) result.fps = updates.config.fps
  }

  return result
}

/**
 * 旧 Scene 部分更新 → DB UpdateVideoScene
 */
export function toDbSceneUpdate(updates: Partial<Scene>): UpdateVideoScene {
  const result: UpdateVideoScene = {}

  if (updates.name !== undefined) result.name = updates.name
  if (updates.compositionId !== undefined) result.composition_id = updates.compositionId
  if (updates.durationInFrames !== undefined) result.duration_in_frames = updates.durationInFrames
  if (updates.props !== undefined) result.props = JSON.stringify(updates.props)

  return result
}

/**
 * 旧 Transition 更新 → DB UpdateVideoScene（transition 字段）
 *
 * 注意：旧 transitions[i] 对应新 scenes[i+1] 的 transition_* 字段
 */
export function toDbTransitionUpdate(updates: Partial<Transition>): UpdateVideoScene {
  const result: UpdateVideoScene = {}

  if (updates.type !== undefined) result.transition_type = updates.type
  if (updates.durationInFrames !== undefined) result.transition_duration = updates.durationInFrames
  if (updates.direction !== undefined) result.transition_direction = updates.direction ?? null

  return result
}
