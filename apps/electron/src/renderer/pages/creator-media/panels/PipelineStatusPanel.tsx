import { useMemo } from 'react'
import { useT } from '@/context/LocaleContext'
import type { Content } from '@sprouty-ai/shared/db/types'
import type { PipelineState, PipelineStage, PipelineStageStatus } from '../types'

/** 流水线阶段定义 */
const PIPELINE_STAGES: { id: PipelineStage; label: string }[] = [
  { id: 'insight', label: '洞察' },
  { id: 'research', label: '研究' },
  { id: 'create', label: '创作' },
  { id: 'adapt', label: '适配' },
  { id: 'ship', label: '分发' },
]

const STAGE_COLORS: Record<PipelineStageStatus, string> = {
  pending: 'bg-muted text-muted-foreground',
  running: 'bg-blue-500 text-white animate-pulse',
  completed: 'bg-green-500 text-white',
  failed: 'bg-red-500 text-white',
  skipped: 'bg-gray-300 text-gray-600 dark:bg-gray-700 dark:text-gray-400',
}

const MODE_LABELS: Record<string, string> = {
  auto: '全自动',
  'semi-auto': '半自动',
  manual: '手动',
}

interface PipelineStatusPanelProps {
  contents: Content[]
}

/**
 * 流水线状态面板 — 显示当前活跃流水线的进度和状态
 * 当前作为可折叠侧边区域使用
 */
export function PipelineStatusPanel({ contents }: PipelineStatusPanelProps) {
  const t = useT()

  /** 解析活跃流水线（非 published/archived 且有 pipeline_state 的内容） */
  const activePipelines = useMemo(() => {
    return contents
      .filter((c) => !['published', 'archived'].includes(c.status) && c.pipeline_state)
      .map((c) => {
        let state: PipelineState | null = null
        try {
          state = JSON.parse(c.pipeline_state!) as PipelineState
        } catch {
          // 解析失败忽略
        }
        return { content: c, state }
      })
      .filter((item) => item.state !== null)
      .slice(0, 3) // 最多显示 3 个
  }, [contents])

  if (activePipelines.length === 0) {
    return (
      <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-4 text-center">
        <p className="text-xs text-muted-foreground">{t('暂无活跃流水线')}</p>
        <p className="mt-1 text-[10px] text-muted-foreground/70">{t('在对话中说"开始创作"启动流水线')}</p>
      </div>
    )
  }

  return (
    <div className="space-y-3">
      {activePipelines.map(({ content, state }) => (
        <div key={content.id} className="rounded-lg border border-border/60 bg-background/40 px-4 py-3 space-y-2">
          {/* 内容标题 */}
          <div className="flex items-center justify-between">
            <span className="text-sm font-medium text-foreground truncate">
              {content.title || t('无标题')}
            </span>
            <span className="inline-flex rounded-full bg-muted/60 px-2 py-0.5 text-[10px] text-muted-foreground flex-shrink-0">
              {t(MODE_LABELS[state!.mode] || state!.mode)}
            </span>
          </div>

          {/* 阶段进度 */}
          <div className="flex items-center gap-1">
            {PIPELINE_STAGES.map((stage, index) => {
              const stageData = state!.stages[stage.id]
              const status = stageData?.status || 'pending'
              const isCurrent = stage.id === state!.currentStage

              return (
                <div key={stage.id} className="flex items-center">
                  <div className="flex flex-col items-center gap-0.5">
                    <div
                      className={`flex h-6 w-6 items-center justify-center rounded-full text-[10px] font-medium ${STAGE_COLORS[status]} ${isCurrent ? 'ring-2 ring-blue-400 ring-offset-1 ring-offset-background' : ''}`}
                    >
                      {status === 'completed' ? '✓' : index + 1}
                    </div>
                    <span className={`text-[9px] ${isCurrent ? 'font-medium text-foreground' : 'text-muted-foreground'}`}>
                      {t(stage.label)}
                    </span>
                  </div>
                  {index < PIPELINE_STAGES.length - 1 && (
                    <div className={`mx-0.5 h-0.5 w-4 ${status === 'completed' ? 'bg-green-500' : 'bg-border'}`} />
                  )}
                </div>
              )
            })}
          </div>

          {/* 等待用户操作提示 */}
          {state!.waitingForUser && (
            <p className="text-[10px] text-amber-600 dark:text-amber-400">{t('等待确认...')}</p>
          )}
        </div>
      ))}
    </div>
  )
}
