/**
 * 视频工作台 — 视频创作与编辑视图
 *
 * 三列布局：
 * - 左面板：视频内容列表 + 模板库（双 Tab）
 * - 中心面板：视频预览 + 时间线
 * - 右面板：属性编辑 + 导出 + 完成制作
 *
 * Phase 5 适配：状态使用 VideoProjectFull (DB)，通过适配层转为旧 VideoProject 传给子组件。
 */

import { useState, useCallback, useRef, useEffect, useMemo } from 'react'
import { Clapperboard, Settings, Download, CheckCircle, Sparkles } from 'lucide-react'
import { toast } from 'sonner'
import type { PlayerRef } from '@remotion/player'
import { Button } from '@/components/ui/button'
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { ScrollArea } from '@/components/ui/scroll-area'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useNavigation, routes } from '@/contexts/NavigationContext'
import type { VideoProjectFull } from '@sprouty-ai/shared/db/types'
import type { VideoTemplate } from '@sprouty-ai/video'
import { VideoPreview } from '@/components/video/VideoPreview'
import { VideoTimeline } from '@/components/video/VideoTimeline'
import { VideoProperties } from '@/components/video/VideoProperties'
import { VideoTemplates } from '@/components/video/VideoTemplates'
import { VideoExport } from '@/components/video/VideoExport'
import { CreateVideoProjectDialog } from '@/components/video/CreateVideoProjectDialog'
import { toLegacyProject, toDbProjectUpdate } from '@/lib/video-adapter'
import { EditPopover, getEditConfig } from '@/components/ui/EditPopover'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { VideoContentList } from './components/VideoContentList'
import type { Content, ContentVideoMetadata } from '@sprouty-ai/shared/db/types'

/** 解析 content.metadata JSON */
function parseVideoMetadata(content: Content): ContentVideoMetadata | null {
  if (!content.metadata) return null
  try {
    return JSON.parse(content.metadata) as ContentVideoMetadata
  } catch {
    return null
  }
}

/** 检测是否是临时预览项目（不持久化到数据库） */
function isTemporaryProject(project: VideoProjectFull | null): boolean {
  return project?.id.startsWith('temp-') ?? false
}

export default function VideoStudio() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const { navigate } = useNavigation()
  const playerRef = useRef<PlayerRef>(null)

  const {
    contents, loading, activeProject,
    updateContentStatus, createContent, refreshContents,
  } = useCreatorMedia()

  // 筛选视频类型内容
  const videoContents = useMemo(
    () => contents.filter(c => {
      if (!c.metadata) return false
      try {
        const meta = JSON.parse(c.metadata)
        return !!meta.video_project_id
      } catch {
        return false
      }
    }),
    [contents],
  )

  // 从 URL 参数获取 contentId
  const searchParams = new URLSearchParams(window.location.search)
  const urlContentId = searchParams.get('contentId')

  // 状态 — 使用 DB 类型
  const [activeContentId, setActiveContentId] = useState<string | null>(null)
  const [currentProjectFull, setCurrentProjectFull] = useState<VideoProjectFull | null>(null)
  const [currentFrame, setCurrentFrame] = useState(0)
  const [isPlaying, setIsPlaying] = useState(false)
  const [leftTab, setLeftTab] = useState<'contents' | 'templates'>('contents')
  const [rightTab, setRightTab] = useState<'properties' | 'export'>('properties')

  // 创建项目对话框
  const [createDialogOpen, setCreateDialogOpen] = useState(false)
  const [pendingTemplate, setPendingTemplate] = useState<VideoTemplate | null>(null)

  // 派生旧类型 — 子组件零改动
  const legacyProject = useMemo(
    () => currentProjectFull ? toLegacyProject(currentProjectFull) : null,
    [currentProjectFull],
  )

  /** 从 IPC 刷新项目数据 */
  const refetchProject = useCallback(async (projectId: string) => {
    if (!window.electronAPI?.video?.getProject) return
    try {
      const full = await window.electronAPI.video.getProject(workspaceId, projectId)
      if (full) setCurrentProjectFull(full)
    } catch (err) {
      console.warn('刷新项目失败:', err)
    }
  }, [workspaceId])

  /** 从视频内容加载/创建 VideoProject */
  const handleContentSelect = useCallback(async (content: Content) => {
    setActiveContentId(content.id)
    const meta = parseVideoMetadata(content)

    if (meta?.video_project_id && window.electronAPI?.video?.getProject) {
      try {
        const proj = await window.electronAPI.video.getProject(workspaceId, meta.video_project_id)
        if (proj) {
          setCurrentProjectFull(proj)
          return
        }
      } catch {
        // 加载失败，走创建流程
      }
    }

    // 没有关联项目或加载失败 — 通过 IPC 创建
    if (window.electronAPI?.video?.createProject) {
      try {
        const resolution = meta?.video_resolution || { width: 1080, height: 1920 }
        const created = await window.electronAPI.video.createProject({
          workspaceId,
          projectId: activeProject?.id || '',
          contentId: content.id,
          name: content.title || t('视频项目'),
          width: resolution.width,
          height: resolution.height,
          fps: 30,
        })
        setCurrentProjectFull(created)
      } catch (err) {
        console.error('创建视频项目失败:', err)
        toast.error(t('创建视频项目失败'))
      }
    }
  }, [workspaceId, t])

  /** 模板预览（纯前端内存构建，不持久化） */
  const handleTemplateSelect = useCallback((template: VideoTemplate) => {
    // 构建临时的 VideoProjectFull 对象（不调用 IPC，不持久化到数据库）
    const tempProject: VideoProjectFull = {
      id: `temp-${Date.now()}`, // 临时 ID
      content_id: '', // 临时项目无 content_id
      name: `${template.name} - ${t('预览')}`,
      description: template.description,
      width: template.defaultConfig.width,
      height: template.defaultConfig.height,
      fps: template.defaultConfig.fps,
      metadata: JSON.stringify({ templateId: template.id, isPreview: true }),
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString(),
      scenes: [
        {
          id: `temp-scene-${Date.now()}`,
          project_id: `temp-${Date.now()}`,
          composition_id: template.compositionId,
          name: template.name,
          sort_order: 0,
          duration_in_frames: template.defaultConfig.durationInFrames,
          props: JSON.stringify(template.defaultProps),
          transition_type: 'none',
          transition_duration: 0,
          transition_direction: null,
          created_at: new Date().toISOString(),
          updated_at: new Date().toISOString(),
        },
      ],
      assets: [],
    }
    setCurrentProjectFull(tempProject)
    setActiveContentId(null)
  }, [t])

  /** 从模板创建独立项目 */
  const handleTemplateCreate = useCallback((template: VideoTemplate) => {
    setPendingTemplate(template)
    setCreateDialogOpen(true)
  }, [])

  /** 确认从模板创建 */
  const handleCreateFromTemplate = useCallback(async (data: { name: string; description: string }) => {
    if (!pendingTemplate) return
    if (!window.electronAPI?.video?.createFromTemplate) return

    try {
      const created = await window.electronAPI.video.createFromTemplate({
        workspaceId,
        projectId: activeProject?.id || '',  // 空字符串时后端会自动获取活跃项目
        contentId: activeContentId || '',    // 空字符串时后端会自动创建 content
        templateId: pendingTemplate.id,
        name: data.name,
        description: data.description || undefined,
      })
      setCurrentProjectFull(created)
      setLeftTab('contents')
      setPendingTemplate(null)
      // 刷新内容列表以显示新创建的视频项目
      await refreshContents()
    } catch (err) {
      console.error('从模板创建项目失败:', err)
      toast.error(t('创建项目失败'))
    }
  }, [pendingTemplate, workspaceId, activeContentId, activeProject, t, refreshContents])

  /** 更新项目属性（通过 IPC 持久化 + refetch） */
  const handleProjectUpdate = useCallback(async (updates: Partial<import('@sprouty-ai/video').VideoProject>) => {
    if (!currentProjectFull) return

    // 临时预览项目不持久化，只在内存中更新
    if (isTemporaryProject(currentProjectFull)) {
      setCurrentProjectFull(prev => prev ? { ...prev, ...updates } : null)
      return
    }

    const dbUpdates = toDbProjectUpdate(updates)
    if (Object.keys(dbUpdates).length === 0) return

    try {
      if (window.electronAPI?.video?.updateProject) {
        await window.electronAPI.video.updateProject(workspaceId, currentProjectFull.id, dbUpdates)
        await refetchProject(currentProjectFull.id)
      }
    } catch (err) {
      console.warn('更新项目失败:', err)
    }
  }, [currentProjectFull, workspaceId, refetchProject])

  /** 帧变化 */
  const handleFrameChange = useCallback((frame: number) => {
    setCurrentFrame(frame)
    if (playerRef.current) {
      playerRef.current.seekTo(frame)
    }
  }, [])

  /** 播放/暂停 */
  const handlePlayPause = useCallback(() => {
    if (playerRef.current) {
      if (isPlaying) {
        playerRef.current.pause()
      } else {
        playerRef.current.play()
      }
      setIsPlaying(!isPlaying)
    }
  }, [isPlaying])

  /** 导出渲染 */
  const handleRender = useCallback(async (options: {
    compositionId: string
    outputFormat: 'mp4' | 'webm' | 'gif'
    quality: 'draft' | 'standard' | 'high'
  }) => {
    if (!currentProjectFull) return
    try {
      if (!window.electronAPI?.video?.render) {
        toast.info(t('视频导出功能需要完整的 Electron 环境'))
        return
      }
      const outputPath = await window.electronAPI.video.render({
        projectId: currentProjectFull.id,
        workspaceId,
        outputFormat: options.outputFormat,
        quality: options.quality,
      })
      if (outputPath) {
        toast.success(t('视频导出成功'), { description: outputPath })
      }
    } catch (error) {
      toast.error(t('视频导出失败'), {
        description: error instanceof Error ? error.message : String(error),
      })
    }
  }, [currentProjectFull, workspaceId, t])

  /** 完成制作 — 将关联的 content 状态推进到 adapting */
  const handleFinishCreation = useCallback(async () => {
    if (!activeContentId) {
      toast.info(t('当前项目未关联内容'))
      return
    }
    try {
      await updateContentStatus(activeContentId, 'adapting')
      toast.success(t('视频制作完成，进入平台适配阶段'))
      navigate(routes.view.appView('creator-media', 'dashboard'))
    } catch {
      toast.error(t('状态更新失败'))
    }
  }, [activeContentId, updateContentStatus, navigate, t])

  // 自动加载 URL 参数指定的内容
  useEffect(() => {
    if (urlContentId && videoContents.length > 0) {
      const content = videoContents.find(c => c.id === urlContentId)
      if (content) {
        handleContentSelect(content)
      }
    }
  }, [urlContentId, videoContents, handleContentSelect])

  // 同步播放器状态
  useEffect(() => {
    const player = playerRef.current
    if (!player) return

    const handlePlay = () => setIsPlaying(true)
    const handlePause = () => setIsPlaying(false)
    const handleSeek = (e: { detail: { frame: number } }) => {
      setCurrentFrame(e.detail.frame)
    }

    player.addEventListener('play', handlePlay)
    player.addEventListener('pause', handlePause)
    player.addEventListener('seeked', handleSeek)

    return () => {
      player.removeEventListener('play', handlePlay)
      player.removeEventListener('pause', handlePause)
      player.removeEventListener('seeked', handleSeek)
    }
  }, [])

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  return (
    <div className="flex h-full">
      {/* 左面板 — 视频内容 + 模板库 */}
      <div className="w-64 border-r flex flex-col shrink-0">
        <div className="p-2 border-b flex items-center justify-between titlebar-no-drag relative z-panel h-[40px]">
          <Tabs value={leftTab} onValueChange={(v) => setLeftTab(v as 'contents' | 'templates')}>
            <TabsList className="h-8">
              <TabsTrigger value="contents" className="text-xs px-2">
                <Clapperboard className="h-3.5 w-3.5 mr-1" />
                {t('视频内容')}
              </TabsTrigger>
              <TabsTrigger value="templates" className="text-xs px-2">
                {t('模板')}
              </TabsTrigger>
            </TabsList>
          </Tabs>
          <EditPopover
            trigger={
              <button className="inline-flex items-center gap-1 rounded-md px-1.5 py-1 text-xs text-muted-foreground hover:text-foreground hover:bg-accent transition-colors">
                <Sparkles className="h-3.5 w-3.5" />
              </button>
            }
            {...getEditConfig('creator-media-create-video-project', activeWorkspace?.rootPath || '')}
          />
        </div>

        <ScrollArea className="flex-1">
          {leftTab === 'contents' ? (
            <VideoContentList
              contents={videoContents}
              selected={activeContentId ?? undefined}
              onSelect={handleContentSelect}
            />
          ) : (
            <VideoTemplates onSelect={handleTemplateSelect} onCreate={handleTemplateCreate} />
          )}
        </ScrollArea>
      </div>

      {/* 中心面板 — 预览 + 时间线 */}
      <div className="flex-1 flex flex-col min-w-0">
        <div className="flex-1 min-h-0">
          <VideoPreview
            project={legacyProject}
            playerRef={playerRef}
            onFrameChange={handleFrameChange}
          />
        </div>
        <div className="h-32 border-t shrink-0">
          <VideoTimeline
            project={legacyProject}
            currentFrame={currentFrame}
            onFrameChange={handleFrameChange}
            isPlaying={isPlaying}
            onPlayPause={handlePlayPause}
          />
        </div>
      </div>

      {/* 右面板 — 属性 + 导出 */}
      <div className="w-72 border-l flex flex-col shrink-0">
        <div className="p-2 border-b titlebar-no-drag relative z-panel h-[40px]">
          <Tabs value={rightTab} onValueChange={(v) => setRightTab(v as 'properties' | 'export')}>
            <TabsList className="h-8 w-full">
              <TabsTrigger value="properties" className="text-xs flex-1">
                <Settings className="h-3.5 w-3.5 mr-1" />
                {t('属性')}
              </TabsTrigger>
              <TabsTrigger value="export" className="text-xs flex-1">
                <Download className="h-3.5 w-3.5 mr-1" />
                {t('导出')}
              </TabsTrigger>
            </TabsList>
          </Tabs>
        </div>

        <ScrollArea className="flex-1">
          {rightTab === 'properties' ? (
            <VideoProperties
              project={legacyProject}
              onUpdate={handleProjectUpdate}
              onRender={() => setRightTab('export')}
            />
          ) : legacyProject ? (
            <div className="flex flex-col">
              <VideoExport
                project={legacyProject}
                onExport={handleRender}
                onCancel={() => setRightTab('properties')}
              />
              {/* 完成制作按钮 — 仅当关联了 content 时显示 */}
              {activeContentId && (
                <div className="px-3 pb-3">
                  <Button
                    className="w-full"
                    variant="default"
                    onClick={handleFinishCreation}
                  >
                    <CheckCircle className="h-4 w-4 mr-1" />
                    {t('完成制作')}
                  </Button>
                  <p className="text-[10px] text-muted-foreground text-center mt-1">
                    {t('将内容状态推进到平台适配阶段')}
                  </p>
                </div>
              )}
            </div>
          ) : (
            <div className="p-4 text-center text-muted-foreground text-sm">
              {t('请先选择视频内容或模板')}
            </div>
          )}
        </ScrollArea>
      </div>

      {/* 从模板创建项目对话框 */}
      <CreateVideoProjectDialog
        open={createDialogOpen}
        onOpenChange={setCreateDialogOpen}
        template={pendingTemplate}
        onConfirm={handleCreateFromTemplate}
      />
    </div>
  )
}
