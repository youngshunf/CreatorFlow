/**
 * VideoEditor - Main video editor component with three-column layout
 *
 * Features:
 * - Left panel: Project list and templates
 * - Center panel: Video preview and timeline
 * - Right panel: Properties and export options
 *
 * @requirements 9.1
 */

import * as React from 'react';
import { useState, useCallback, useRef, useEffect } from 'react';
import { Plus, Film, Settings, Download } from 'lucide-react';
import { toast } from 'sonner';
import type { PlayerRef } from '@remotion/player';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Tabs, TabsList, TabsTrigger, TabsContent } from '@/components/ui/tabs';
import { ScrollArea } from '@/components/ui/scroll-area';
import { useT } from '@/context/LocaleContext';
import type { VideoProject, VideoTemplate } from '@creator-flow/video';
import { VideoPreview } from './VideoPreview';
import { VideoTimeline } from './VideoTimeline';
import { VideoProperties } from './VideoProperties';
import { VideoProjectList } from './VideoProjectList';
import { VideoTemplates } from './VideoTemplates';
import { VideoExport } from './VideoExport';
import { CreateVideoProjectDialog } from './CreateVideoProjectDialog';

export interface VideoEditorProps {
  /** Workspace ID for project management */
  workspaceId: string;
  /** Optional class name */
  className?: string;
}

/**
 * Map template IDs to composition component IDs
 */
const TEMPLATE_COMPOSITION_MAP: Record<string, string> = {
  'social-media-vertical': 'TitleAnimation',
  'social-media-square': 'TitleAnimation',
  'marketing-product': 'ProductShowcase',
  'marketing-promo': 'TitleAnimation',
  'tutorial-steps': 'Slideshow',
  'tutorial-explainer': 'ProductShowcase',
  'tutorial-tips': 'Slideshow',
};

/**
 * VideoEditor component
 */
export function VideoEditor({ workspaceId, className }: VideoEditorProps) {
  const t = useT();
  const playerRef = useRef<PlayerRef>(null);

  // State
  const [currentProject, setCurrentProject] = useState<VideoProject | null>(null);
  const [localProjects, setLocalProjects] = useState<VideoProject[]>([]);
  const [currentFrame, setCurrentFrame] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [showExport, setShowExport] = useState(false);
  const [leftTab, setLeftTab] = useState<'projects' | 'templates'>('projects');
  const [rightTab, setRightTab] = useState<'properties' | 'export'>('properties');

  // Dialog state for creating project from template
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [pendingTemplate, setPendingTemplate] = useState<VideoTemplate | null>(null);

  // Handle project selection
  const handleProjectSelect = useCallback((project: VideoProject) => {
    setCurrentProject(project);
    setCurrentFrame(0);
    setIsPlaying(false);
  }, []);

  // Handle template selection - preview only (not added to project list)
  const handleTemplateSelect = useCallback((template: VideoTemplate) => {
    const now = new Date().toISOString();
    const compositionId = TEMPLATE_COMPOSITION_MAP[template.id] || 'TitleAnimation';
    const project: VideoProject = {
      id: `proj-${Date.now()}`,
      name: `${template.name} - ${new Date().toLocaleDateString()}`,
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
    };
    setCurrentProject(project);
  }, []);

  // Handle template create - open dialog to fill project info
  const handleTemplateCreate = useCallback((template: VideoTemplate) => {
    setPendingTemplate(template);
    setCreateDialogOpen(true);
  }, []);

  // Handle confirm from create dialog
  const handleCreateFromTemplate = useCallback((data: { name: string; description: string }) => {
    if (!pendingTemplate) return;
    const now = new Date().toISOString();
    const compositionId = TEMPLATE_COMPOSITION_MAP[pendingTemplate.id] || 'TitleAnimation';
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
    };
    setLocalProjects(prev => [...prev, project]);
    setCurrentProject(project);
    setLeftTab('projects');
    setPendingTemplate(null);
  }, [pendingTemplate]);

  // Handle project update
  const handleProjectUpdate = useCallback((updates: Partial<VideoProject>) => {
    if (!currentProject) return;
    setCurrentProject({
      ...currentProject,
      ...updates,
      updatedAt: new Date().toISOString(),
    });
  }, [currentProject]);

  // Handle frame change
  const handleFrameChange = useCallback((frame: number) => {
    setCurrentFrame(frame);
    if (playerRef.current) {
      playerRef.current.seekTo(frame);
    }
  }, []);

  // Handle play/pause
  const handlePlayPause = useCallback(() => {
    if (playerRef.current) {
      if (isPlaying) {
        playerRef.current.pause();
      } else {
        playerRef.current.play();
      }
      setIsPlaying(!isPlaying);
    }
  }, [isPlaying]);

  // Handle render - connect to real export API
  const handleRender = useCallback(async (options: {
    compositionId: string;
    outputFormat: 'mp4' | 'webm' | 'gif';
    quality: 'draft' | 'standard' | 'high';
  }) => {
    if (!currentProject) return;
    try {
      if (!window.electronAPI?.video?.render) {
        toast.info(t('视频导出功能需要完整的 Electron 环境'));
        return;
      }
      const outputPath = await window.electronAPI.video.render({
        projectId: currentProject.id,
        compositionId: options.compositionId,
        outputFormat: options.outputFormat,
        quality: options.quality,
      });
      if (outputPath) {
        toast.success(t('视频导出成功'), { description: outputPath });
      }
    } catch (error) {
      toast.error(t('视频导出失败'), {
        description: error instanceof Error ? error.message : String(error),
      });
    }
  }, [currentProject, t]);

  // Create new project
  const handleCreateProject = useCallback(() => {
    const now = new Date().toISOString();
    const project: VideoProject = {
      id: `proj-${Date.now()}`,
      name: t('新视频项目'),
      createdAt: now,
      updatedAt: now,
      config: {
        width: 1920,
        height: 1080,
        fps: 30,
        durationInFrames: 300,
      },
      compositions: [
        {
          id: 'TitleAnimation',
          name: t('标题动画'),
          code: '',
          props: {
            title: t('欢迎'),
            subtitle: t('在此输入副标题'),
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
    };
    setCurrentProject(project);
  }, [t]);

  // Sync player state
  useEffect(() => {
    const player = playerRef.current;
    if (!player) return;

    const handlePlay = () => setIsPlaying(true);
    const handlePause = () => setIsPlaying(false);
    const handleSeek = (e: { detail: { frame: number } }) => {
      setCurrentFrame(e.detail.frame);
    };

    player.addEventListener('play', handlePlay);
    player.addEventListener('pause', handlePause);
    player.addEventListener('seeked', handleSeek as EventListener);

    return () => {
      player.removeEventListener('play', handlePlay);
      player.removeEventListener('pause', handlePause);
      player.removeEventListener('seeked', handleSeek as EventListener);
    };
  }, []);

  return (
    <div className={cn('flex h-full', className)}>
      {/* Left Panel - Projects & Templates */}
      <div className="w-64 border-r flex flex-col shrink-0">
        <div className="p-2 border-b flex items-center justify-between titlebar-no-drag relative z-panel h-[40px]">
          <Tabs value={leftTab} onValueChange={(v) => setLeftTab(v as 'projects' | 'templates')}>
            <TabsList className="h-8">
              <TabsTrigger value="projects" className="text-xs px-2">
                <Film className="h-3.5 w-3.5 mr-1" />
                {t('项目')}
              </TabsTrigger>
              <TabsTrigger value="templates" className="text-xs px-2">
                {t('模板')}
              </TabsTrigger>
            </TabsList>
          </Tabs>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleCreateProject}
            title={t('新建项目')}
          >
            <Plus className="h-4 w-4" />
          </Button>
        </div>

        <ScrollArea className="flex-1">
          {leftTab === 'projects' ? (
            <VideoProjectList
              workspaceId={workspaceId}
              selected={currentProject?.id}
              onSelect={handleProjectSelect}
              extraProjects={localProjects}
            />
          ) : (
            <VideoTemplates onSelect={handleTemplateSelect} onCreate={handleTemplateCreate} />
          )}
        </ScrollArea>
      </div>

      {/* Center Panel - Preview & Timeline */}
      <div className="flex-1 flex flex-col min-w-0">
        {/* Preview */}
        <div className="flex-1 min-h-0">
          <VideoPreview
            project={currentProject}
            playerRef={playerRef}
            onFrameChange={handleFrameChange}
          />
        </div>

        {/* Timeline */}
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

      {/* Right Panel - Properties & Export */}
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
            <VideoExport
              project={currentProject}
              onExport={handleRender}
              onCancel={() => setRightTab('properties')}
            />
          ) : (
            <div className="p-4 text-center text-muted-foreground text-sm">
              {t('请先选择或创建一个项目')}
            </div>
          )}
        </ScrollArea>
      </div>

      {/* Create project from template dialog */}
      <CreateVideoProjectDialog
        open={createDialogOpen}
        onOpenChange={setCreateDialogOpen}
        template={pendingTemplate}
        onConfirm={handleCreateFromTemplate}
      />
    </div>
  );
}
