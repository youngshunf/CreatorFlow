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

export interface VideoEditorProps {
  /** Workspace ID for project management */
  workspaceId: string;
  /** Optional class name */
  className?: string;
}

/**
 * VideoEditor component
 */
export function VideoEditor({ workspaceId, className }: VideoEditorProps) {
  const t = useT();
  const playerRef = useRef<PlayerRef>(null);

  // State
  const [currentProject, setCurrentProject] = useState<VideoProject | null>(null);
  const [currentFrame, setCurrentFrame] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [showExport, setShowExport] = useState(false);
  const [leftTab, setLeftTab] = useState<'projects' | 'templates'>('projects');
  const [rightTab, setRightTab] = useState<'properties' | 'export'>('properties');

  // Handle project selection
  const handleProjectSelect = useCallback((project: VideoProject) => {
    setCurrentProject(project);
    setCurrentFrame(0);
    setIsPlaying(false);
  }, []);

  // Handle template selection - create new project from template
  const handleTemplateSelect = useCallback(async (template: VideoTemplate) => {
    try {
      const project = await window.electronAPI.video.createProject({
        name: `${template.name} - ${new Date().toLocaleDateString()}`,
        workspaceId,
        template: template.id,
        config: template.defaultConfig,
      });
      setCurrentProject(project);
      setLeftTab('projects');
    } catch (error) {
      console.error('Failed to create project from template:', error);
    }
  }, [workspaceId]);

  // Handle project update
  const handleProjectUpdate = useCallback(async (updates: Partial<VideoProject>) => {
    if (!currentProject) return;
    try {
      const updated = await window.electronAPI.video.updateProject(
        currentProject.id,
        updates
      );
      setCurrentProject(updated);
    } catch (error) {
      console.error('Failed to update project:', error);
    }
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

  // Handle render
  const handleRender = useCallback(async (options: {
    compositionId: string;
    outputFormat: 'mp4' | 'webm' | 'gif';
    quality: 'draft' | 'standard' | 'high';
  }) => {
    if (!currentProject) return;
    try {
      await window.electronAPI.video.render({
        projectId: currentProject.id,
        ...options,
      });
      setShowExport(false);
    } catch (error) {
      console.error('Failed to render:', error);
    }
  }, [currentProject]);

  // Create new project
  const handleCreateProject = useCallback(async () => {
    try {
      const project = await window.electronAPI.video.createProject({
        name: t('新视频项目'),
        workspaceId,
        config: {
          width: 1920,
          height: 1080,
          fps: 30,
          durationInFrames: 300,
        },
      });
      setCurrentProject(project);
    } catch (error) {
      console.error('Failed to create project:', error);
    }
  }, [workspaceId, t]);

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
        <div className="p-2 border-b flex items-center justify-between">
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
            />
          ) : (
            <VideoTemplates onSelect={handleTemplateSelect} />
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
        <div className="p-2 border-b">
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
    </div>
  );
}
