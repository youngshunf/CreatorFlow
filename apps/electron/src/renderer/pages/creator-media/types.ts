/**
 * 创作面板视图相关类型定义
 */

/** 统计数据 */
export interface ContentStats {
  total: number
  idea: number
  researching: number
  scripting: number
  creating: number
  reviewing: number
  scheduled: number
  published: number
  archived: number
}

/** 流水线阶段 */
export type PipelineStage = 'insight' | 'research' | 'create' | 'adapt' | 'ship'

/** 流水线阶段状态 */
export type PipelineStageStatus = 'pending' | 'running' | 'completed' | 'failed' | 'skipped'

/** 流水线状态 */
export interface PipelineState {
  mode: 'auto' | 'semi-auto' | 'manual'
  currentStage: PipelineStage
  startedAt: string
  stages: Record<PipelineStage, { status: PipelineStageStatus; output?: unknown }>
  waitingForUser: boolean
}

/** 发布数据统计 */
export interface PublishStats {
  totalViews: number
  totalLikes: number
  totalComments: number
  totalFavorites: number
  totalShares: number
  successCount: number
  failedCount: number
}
