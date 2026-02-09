import { useT } from '@/context/LocaleContext'
import type { PipelineStage, PipelineStageStatus } from '../types'

/** 流水线阶段定义 */
const PIPELINE_STAGES: { id: PipelineStage; label: string }[] = [
  { id: 'insight', label: '洞察' },
  { id: 'research', label: '研究' },
  { id: 'create', label: '创作' },
  { id: 'adapt', label: '适配' },
  { id: 'ship', label: '分发' },
]

/** 阶段状态颜色 */
const STAGE_COLORS: Record<PipelineStageStatus, string> = {
  pending: 'bg-muted text-muted-foreground',
  running: 'bg-blue-500 text-white animate-pulse',
  completed: 'bg-green-500 text-white',
  failed: 'bg-red-500 text-white',
  skipped: 'bg-gray-300 text-gray-600 dark:bg-gray-700 dark:text-gray-400',
}

interface PipelineProgressProps {
  stages: Record<PipelineStage, { status: PipelineStageStatus; output?: unknown }>
  currentStage: PipelineStage
}

export function PipelineProgress({ stages, currentStage }: PipelineProgressProps) {
  const t = useT()

  return (
    <div className="flex items-center gap-1">
      {PIPELINE_STAGES.map((stage, index) => {
        const stageData = stages[stage.id]
        const status = stageData?.status || 'pending'
        const isCurrent = stage.id === currentStage

        return (
          <div key={stage.id} className="flex items-center">
            {/* 阶段节点 */}
            <div className="flex flex-col items-center gap-1">
              <div
                className={`flex h-8 w-8 items-center justify-center rounded-full text-xs font-medium ${STAGE_COLORS[status]} ${isCurrent ? 'ring-2 ring-blue-400 ring-offset-1 ring-offset-background' : ''}`}
              >
                {status === 'completed' ? '✓' : index + 1}
              </div>
              <span className={`text-[10px] ${isCurrent ? 'font-medium text-foreground' : 'text-muted-foreground'}`}>
                {t(stage.label)}
              </span>
            </div>
            {/* 连接线 */}
            {index < PIPELINE_STAGES.length - 1 && (
              <div className={`mx-1 h-0.5 w-6 ${status === 'completed' ? 'bg-green-500' : 'bg-border'}`} />
            )}
          </div>
        )
      })}
    </div>
  )
}

/** 空状态的流水线进度（全部 pending） */
export function EmptyPipelineProgress() {
  const t = useT()
  const emptyStages: Record<PipelineStage, { status: PipelineStageStatus }> = {
    insight: { status: 'pending' },
    research: { status: 'pending' },
    create: { status: 'pending' },
    adapt: { status: 'pending' },
    ship: { status: 'pending' },
  }

  return (
    <div className="flex flex-col items-center gap-3 py-4">
      <PipelineProgress stages={emptyStages} currentStage="insight" />
      <p className="text-xs text-muted-foreground">{t('点击"开始创作"启动流水线')}</p>
    </div>
  )
}
