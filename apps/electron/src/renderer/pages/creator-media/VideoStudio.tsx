/**
 * 视频工作台 — 视频创作与编辑视图
 *
 * 三列布局：
 * - 左面板：视频内容列表 + 模板库（双 Tab）
 * - 中心面板：视频预览 + 时间线
 * - 右面板：属性编辑 + 导出 + 完成制作
 *
 * 复用 VideoEditor 子组件（VideoPreview, VideoTimeline, VideoProperties, VideoExport, VideoTemplates）
 */

import { useState, useCallback, useRef, useEffect, useMemo } from 'react'
import { Clapperboard, Film, Settings, Download, Plus, CheckCircle } from 'lucide-react'
import { toast } from 'sonner'
import type { PlayerRef } from '@remotion/player'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { ScrollArea } from '@/components/ui/scroll-area'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useNavigation, routes } from '@/contexts/NavigationContext'
import type { VideoProject, VideoTemplate } from '@sprouty-ai/video'
import { VideoPreview } from '@/components/video/VideoPreview'
import { VideoTimeline } from '@/components/video/VideoTimeline'
import { VideoProperties } from '@/components/video/VideoProperties'
import { VideoTemplates } from '@/components/video/VideoTemplates'
import { VideoExport } from '@/components/video/VideoExport'
import { CreateVideoProjectDialog } from '@/components/video/CreateVideoProjectDialog'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { VideoContentList } from './components/VideoContentList'
import type { Content, ContentVideoMetadata } from '@sprouty-ai/shared/db/types'

/** 模板 ID → Composition 组件 ID 映射 */
const TEMPLATE_COMPOSITION_MAP: Record<string, string> = {
  'social-media-vertical': 'TitleAnimation',
  'social-media-square': 'TitleAnimation',
  'marketing-product': 'ProductShowcase',
  'marketing-promo': 'TitleAnimation',
  'tutorial-steps': 'Slideshow',
  'tutorial-explainer': 'ProductShowcase',
  'tutorial-tips': 'Slideshow',
}

/** 解析 content.metadata JSON */
function parseVideoMetadata(content: Content): ContentVideoMetadata | null {
  if (!content.metadata) return null
  try {
    return JSON.parse(content.metadata) as ContentVideoMetadata
  } catch {
    return null
  }
}

export default function VideoStudio() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const { navigate } = useNavigation()
  const playerRef = useRef<PlayerRef>(null)

  const {
    contents, loading, activeProject,
    updateContentStatus, createContent,
  } = useCreatorMedia()

  // 筛选视频类型内容（通过 metadata 判断）
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

  // 状态
  const [activeContentId, setActiveContentId] = useState<string | null>(null)
  const [currentProject, setCurrentProject] = useState<VideoProject | null>(null)
  const [localProjects, setLocalProjects] = useState<VideoProject[]>([])
  const [currentFrame, setCurrentFrame] = useState(0)
  const [isPlaying, setIsPlaying] = useState(false)
  const [leftTab, setLeftTab] = useState<'contents' | 'templates'>('contents')
  const [rightTab, setRightTab] = useState<'properties' | 'export'>('properties')

  // 创建项目对话框
  const [createDialogOpen, setCreateDialogOpen] = useState(false)
  const [pendingTemplate, setPendingTemplate] = useState<VideoTemplate | null>(null)

  /** 从视频内容加载/创建 VideoProject */
  const handleContentSelect = useCallback((content: Content) => {
    setActiveContentId(content.id)
    const meta = parseVideoMetadata(content)

    if (meta?.video_project_id) {
      // 已有关联项目 — 尝试从 IPC 加载
      if (window.electronAPI?.video?.getProject) {
        window.electronAPI.video.getProject(meta.video_project_id).then((proj: VideoProject | null) => {
          if (proj) {
            setCurrentProject(proj)
          } else {
            // 项目不存在，创建新的
            createProjectFromContent(content, meta)
          }
        }).catch(() => {
          createProjectFromContent(content, meta)
        })
      } else {
        createProjectFromContent(content, meta)
      }
    } else {
      // 没有关联项目 — 自动创建
      createProjectFromContent(content, meta)
    }
  }, [])

  /** 基于内容创建 VideoProject */
  const createProjectFromContent = useCallback((content: Content, meta: ContentVideoMetadata | null) => {
    const now = new Date().toISOString()
    const templateId = meta?.video_template_id || 'social-media-vertical'
    const compositionId = TEMPLATE_COMPOSITION_MAP[templateId] || 'TitleAnimation'
    const resolution = meta?.video_resolution || { width: 1080, height: 1920 }

    const project: VideoProject = {
      id: `proj-${Date.now()}`,
      name: content.title || t('视频项目'),
      createdAt: now,
      updatedAt: now,
      config: {
        width: resolution.width,
        height: resolution.height,
        fps: 30,
        durationInFrames: 300,
      },
      compositions: [
        {
          id: compositionId,
          name: content.title || t('视频组合'),
          code: '',
          props: {
            title: content.title || '',
            subtitle: content.topic || '',
            colors: {
              primary: '#6366f1',
              secondary: '#8b5cf6',
              background: '#1a1a2e',
              text: '#ffffff',
            },
            animationStyle: 'spring',
          },
        },
      ],
      assets: [],
      renders: [],
    }

    setCurrentProject(project)
    setLocalProjects(prev => [...prev, project])
  }, [t])

  /** 模板预览 */
  const handleTemplateSelect = useCallback((template: VideoTemplate) => {
    const now = new Date().toISOString()
    const compositionId = TEMPLATE_COMPOSITION_MAP[template.id] || 'TitleAnimation'
    const project: VideoProject = {
      id: `preview-${Date.now()}`,
      name: `${template.name} - ${t('预览')}`,
      createdAt: now,
      updatedAt: now,
      config: {
        width: template.defaultConfig.width,
        height: template.defaultConfig.height,
        fps: template.defaultConfig.fps,
        durationInFrames: template.defaultConfig.durationInFrames,
      },
      compositions: [
        {
          id: compositionId,
          name: template.name,
          code: template.compositionCode || '',
          props: template.defaultProps || {},
        },
      ],
      assets: [],
      renders: [],
    }
    setCurrentProject(project)
    setActiveContentId(null)
  }, [t])

  /** 从模板创建独立项目 */
  const handleTemplateCreate = useCallback((template: VideoTemplate) => {
    setPendingTemplate(template)
    setCreateDialogOpen(true)
  }, [])

  /** 确认从模板创建 */
  const handleCreateFromTemplate = useCallback((data: { name: string; description: string }) => {
    if (!pendingTemplate) return
    const now = new Date().toISOString()
    const compositionId = TEMPLATE_COMPOSITION_MAP[pendingTemplate.id] || 'TitleAnimation'
    const project: VideoProject = {
      id: `proj-${Date.now()}`,
      name: data.name,
      description: data.description || undefined,
      createdAt: now,
      updatedAt: now,
      config: {
        width: pendingTemplate.defaultConfig.width,
        height: pendingTemplate.defaultConfig.height,
        fps: pendingTemplate.defaultConfig.fps,
        durationInFrames: pendingTemplate.defaultConfig.durationInFrames,
      },
      compositions: [
        {
          id: compositionId,
          name: pendingTemplate.name,
          code: pendingTemplate.compositionCode || '',
          props: pendingTemplate.defaultProps || {},
        },
      ],
      assets: [],
      renders: [],
    }
    setLocalProjects(prev => [...prev, project])
    setCurrentProject(project)
    setLeftTab('contents')
    setPendingTemplate(null)
    setActiveContentId(null)
  }, [pendingTemplate])

  /** 更新项目属性 */
  const handleProjectUpdate = useCallback((updates: Partial<VideoProject>) => {
    if (!currentProject) return
    setCurrentProject({
      ...currentProject,
      ...updates,
      updatedAt: new Date().toISOString(),
    })
  }, [currentProject])

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
    if (!currentProject) return
    try {
      if (!window.electronAPI?.video?.render) {
        toast.info(t('视频导出功能需要完整的 Electron 环境'))
        return
      }
      const outputPath = await window.electronAPI.video.render({
        projectId: currentProject.id,
        compositionId: options.compositionId,
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
  }, [currentProject, t])

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
    player.addEventListener('seeked', handleSeek as EventListener)

    return () => {
      player.removeEventListener('play', handlePlay)
      player.removeEventListener('pause', handlePause)
      player.removeEventListener('seeked', handleSeek as EventListener)
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
            project={currentProject}
            playerRef={playerRef}
            onFrameChange={handleFrameChange}
          />
        </div>
        <div className="h-32 border-t shrink-0">
          <VideoTimeline
            project={currentProject}
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
              project={currentProject}
              onUpdate={handleProjectUpdate}
              onRender={() => setRightTab('export')}
            />
          ) : currentProject ? (
            <div className="flex flex-col">
              <VideoExport
                project={currentProject}
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
