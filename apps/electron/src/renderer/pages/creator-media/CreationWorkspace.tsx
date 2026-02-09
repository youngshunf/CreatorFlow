import { useState } from 'react'
import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { EmptyPipelineProgress } from './components/PipelineProgress'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { ContentTable } from './components/ContentTable'
import { CreateContentDialog } from './components/CreateContentDialog'

/**
 * 创作工作台 — 流水线驱动的创作视图
 * 展示活跃流水线进度、选题池、进行中的内容
 */
export default function CreationWorkspace() {
  const t = useT()
  const {
    projects, activeProject, contents, loading,
    switchProject, createContent, deleteContent, updateContentStatus,
  } = useCreatorMedia()

  const [showCreateContent, setShowCreateContent] = useState(false)

  // 筛选进行中的内容（非 published/archived）
  const activeContents = contents.filter(
    (c) => !['published', 'archived'].includes(c.status)
  )

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* 头部 — relative z-panel 确保在 titlebar 拖拽区域之上 */}
      <div className="relative z-panel flex items-center justify-between px-6 py-4 border-b border-border/40">
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('创作工作台')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('流水线驱动的创作流程')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <button
            type="button"
            onClick={() => setShowCreateContent(true)}
            className="inline-flex items-center gap-1 rounded-md border border-border/60 bg-background px-2.5 py-1 text-xs text-foreground hover:bg-muted/40 transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
            </svg>
            {t('新建内容')}
          </button>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-6">
        {/* 活跃流水线 */}
        <div>
          <h2 className="text-sm font-medium text-foreground mb-3">{t('活跃流水线')}</h2>
          <div className="rounded-lg border border-border/60 bg-background/40 px-4 py-4">
            <EmptyPipelineProgress />
          </div>
        </div>

        {/* 进行中的内容 */}
        <div>
          <h2 className="text-sm font-medium text-foreground mb-3">{t('进行中的内容')}</h2>
          <ContentTable
            contents={activeContents}
            maxItems={10}
            onStatusChange={updateContentStatus}
            onDelete={deleteContent}
          />
        </div>
      </div>

      {/* 新建内容对话框 */}
      {activeProject && (
        <CreateContentDialog
          open={showCreateContent}
          onOpenChange={setShowCreateContent}
          onCreateContent={createContent}
        />
      )}
    </div>
  )
}
