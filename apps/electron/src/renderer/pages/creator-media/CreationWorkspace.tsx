import { useState, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useNavigation, routes } from '@/contexts/NavigationContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { EmptyPipelineProgress } from './components/PipelineProgress'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { ContentTable } from './components/ContentTable'
import { CreateContentDialog } from './components/CreateContentDialog'
import { VersionHistoryDialog } from './components/VersionHistoryDialog'
import { VideoPreviewDialog } from './components/VideoPreviewDialog'
import { useActiveWorkspace } from '@/context/AppShellContext'
import type { Content } from '@sprouty-ai/shared/db/types'

/** 判断是否为视频类型内容（通过 metadata 判断） */
const isVideoType = (c: Content) => {
  if (!c.metadata) return false
  try {
    const meta = JSON.parse(c.metadata)
    return !!meta.video_project_id
  } catch {
    return false
  }
}

/** 从 metadata JSON 中提取视频输出路径 */
function getVideoOutputPath(content: Content): string {
  if (!content.metadata) return ''
  try {
    const meta = JSON.parse(content.metadata)
    return meta.video_output_path || ''
  } catch {
    return ''
  }
}

/** 下一阶段操作映射 */
const NEXT_STAGE_ACTION: Record<string, { skillId: string; label: string; icon: string } | null> = {
  researching: { skillId: 'idea-researcher', label: '灵感调研', icon: '🔍' },
  scripting: { skillId: 'script-create', label: '脚本创作', icon: '✍️' },
  creating: { skillId: 'content-creator', label: '内容创作', icon: '🎨' },
  adapting: { skillId: 'platform-adapter', label: '平台适配', icon: '📱' },
  scheduled: null,
  published: null,
  archived: null,
}

/**
 * 创作工作台 — 流水线驱动的创作视图
 * 展示活跃流水线进度、选题池、进行中的内容
 */
export default function CreationWorkspace() {
  const t = useT()
  const { navigate } = useNavigation()
  const workspace = useActiveWorkspace()
  const {
    projects, activeProject, contents, loading, wsRoot,
    switchProject, createContent, deleteContent, updateContentStatus,
    listContentVersions, rollbackContentVersion,
  } = useCreatorMedia()

  const [showCreateContent, setShowCreateContent] = useState(false)
  const [versionHistoryContent, setVersionHistoryContent] = useState<Content | null>(null)
  const [previewContent, setPreviewContent] = useState<Content | null>(null)

  // 筛选进行中的内容（非 published/archived）
  const activeContents = contents.filter(
    (c) => !['published', 'archived'].includes(c.status)
  )

  /** 点击内容行：视频制作中跳转工作台，审核中打开预览 */
  const handleContentRowClick = useCallback((content: Content) => {
    if (isVideoType(content) && content.status === 'creating') {
      navigate(routes.view.appView('app.creator-media', 'video-studio'))
    } else if (isVideoType(content) && content.status === 'reviewing') {
      setPreviewContent(content)
    }
  }, [navigate])

  // 处理下一阶段操作
  const handleNextStage = useCallback((content: Content, skillId: string) => {
    if (!workspace || !activeProject) return

    // 获取 skill 的显示名称
    const skillAction = NEXT_STAGE_ACTION[content.status]
    const skillLabel = skillAction?.label || '下一步'

    // 构建提示词
    let prompt = `为内容「${content.title}」执行${skillLabel}。\n\n`
    prompt += `内容 ID: ${content.id}\n`

    // 导航到新 session 并激活 skill
    navigate(routes.view.appView('app.creator-media', 'chat', {
      input: `[skill:${workspace.id}:${skillId}] ${prompt}`,
      send: true,
    }))
  }, [workspace, activeProject, navigate])

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
            onRowClick={handleContentRowClick}
            onOpenVideoStudio={() => navigate(routes.view.appView('app.creator-media', 'video-studio'))}
            onStatusChange={updateContentStatus}
            onDelete={deleteContent}
            onVersionHistory={setVersionHistoryContent}
            onNextStage={handleNextStage}
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

      <VersionHistoryDialog
        open={versionHistoryContent !== null}
        onOpenChange={(open) => { if (!open) setVersionHistoryContent(null) }}
        contentId={versionHistoryContent?.id || ''}
        contentTitle={versionHistoryContent?.title || null}
        listContentVersions={listContentVersions}
        rollbackContentVersion={rollbackContentVersion}
      />

      <VideoPreviewDialog
        open={previewContent !== null}
        onOpenChange={(open) => { if (!open) setPreviewContent(null) }}
        videoOutputPath={previewContent ? getVideoOutputPath(previewContent) : ''}
        contentTitle={previewContent?.title || ''}
        onApprove={async () => {
          if (previewContent) {
            await updateContentStatus(previewContent.id, 'scheduled')
            setPreviewContent(null)
          }
        }}
        onReject={async () => {
          if (previewContent) {
            await updateContentStatus(previewContent.id, 'creating')
            setPreviewContent(null)
          }
        }}
      />
    </div>
  )
}
