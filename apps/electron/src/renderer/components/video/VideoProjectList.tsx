/**
 * VideoProjectList - Project list component
 *
 * Features:
 * - List all video projects in workspace
 * - Project selection
 * - Project deletion
 * - Project metadata display
 *
 * @requirements 9.5
 */

import * as React from 'react';
import { useState, useEffect, useCallback, useMemo } from 'react';
import { Film, Trash2, MoreVertical, Calendar, Clock } from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu';
import { useT } from '@/context/LocaleContext';
import type { VideoProject } from '@sprouty-ai/video';

export interface VideoProjectListProps {
  /** Workspace ID to list projects from */
  workspaceId: string;
  /** Currently selected project ID */
  selected?: string;
  /** Callback when project is selected */
  onSelect: (project: VideoProject) => void;
  /** Extra projects created locally (e.g. from templates) */
  extraProjects?: VideoProject[];
  /** Optional class name */
  className?: string;
}

/**
 * Format date for display
 */
function formatDate(dateString: string): string {
  const date = new Date(dateString);
  return date.toLocaleDateString(undefined, {
    month: 'short',
    day: 'numeric',
  });
}

/**
 * Format duration for display
 */
function formatDuration(frames: number, fps: number): string {
  const seconds = Math.floor(frames / fps);
  const minutes = Math.floor(seconds / 60);
  const remainingSeconds = seconds % 60;
  if (minutes > 0) {
    return `${minutes}:${remainingSeconds.toString().padStart(2, '0')}`;
  }
  return `${seconds}s`;
}

/**
 * Project item component
 */
function ProjectItem({
  project,
  isSelected,
  onSelect,
  onDelete,
}: {
  project: VideoProject;
  isSelected: boolean;
  onSelect: () => void;
  onDelete: () => void;
}) {
  const t = useT();

  return (
    <div
      className={cn(
        'group flex items-start gap-3 p-3 rounded-lg cursor-pointer transition-colors',
        'hover:bg-muted/50',
        isSelected && 'bg-primary/10 hover:bg-primary/15'
      )}
      onClick={onSelect}
    >
      {/* Icon */}
      <div
        className={cn(
          'shrink-0 w-10 h-10 rounded-md flex items-center justify-center',
          isSelected ? 'bg-primary text-primary-foreground' : 'bg-muted'
        )}
      >
        <Film className="h-5 w-5" />
      </div>

      {/* Info */}
      <div className="flex-1 min-w-0">
        <div className="font-medium text-sm truncate">{project.name}</div>
        <div className="flex items-center gap-3 mt-1 text-xs text-muted-foreground">
          <span className="flex items-center gap-1">
            <Calendar className="h-3 w-3" />
            {formatDate(project.updatedAt)}
          </span>
          <span className="flex items-center gap-1">
            <Clock className="h-3 w-3" />
            {formatDuration(project.config.durationInFrames, project.config.fps)}
          </span>
        </div>
        {project.description && (
          <div className="text-xs text-muted-foreground mt-1 truncate">
            {project.description}
          </div>
        )}
      </div>

      {/* Actions */}
      <DropdownMenu>
        <DropdownMenuTrigger asChild>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8 opacity-0 group-hover:opacity-100 transition-opacity"
            onClick={(e) => e.stopPropagation()}
          >
            <MoreVertical className="h-4 w-4" />
          </Button>
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end">
          <DropdownMenuItem
            className="text-destructive focus:text-destructive"
            onClick={(e) => {
              e.stopPropagation();
              onDelete();
            }}
          >
            <Trash2 className="h-4 w-4 mr-2" />
            {t('删除')}
          </DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
    </div>
  );
}

/**
 * VideoProjectList component
 */
export function VideoProjectList({
  workspaceId,
  selected,
  onSelect,
  extraProjects = [],
  className,
}: VideoProjectListProps) {
  const t = useT();
  const [remoteProjects, setRemoteProjects] = useState<VideoProject[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // Merge remote projects with extra (locally created) projects, sorted by updatedAt
  const projects = useMemo(() => {
    const merged = [...remoteProjects];
    for (const ep of extraProjects) {
      if (!merged.some((p) => p.id === ep.id)) {
        merged.push(ep);
      }
    }
    return merged.sort(
      (a, b) => new Date(b.updatedAt).getTime() - new Date(a.updatedAt).getTime()
    );
  }, [remoteProjects, extraProjects]);

  // Load projects from backend
  const loadProjects = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    try {
      if (!window.electronAPI?.video?.listProjects) {
        setRemoteProjects([]);
        return;
      }
      const list = await window.electronAPI.video.listProjects(workspaceId);
      setRemoteProjects(list);
    } catch (err) {
      console.error('Failed to load projects:', err);
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsLoading(false);
    }
  }, [workspaceId]);

  // Initial load
  useEffect(() => {
    loadProjects();
  }, [loadProjects]);

  // Handle delete
  const handleDelete = useCallback(
    async (projectId: string) => {
      try {
        if (!window.electronAPI?.video?.deleteProject) return;
        await window.electronAPI.video.deleteProject(projectId);
        setRemoteProjects((prev) => prev.filter((p) => p.id !== projectId));
      } catch (err) {
        console.error('Failed to delete project:', err);
      }
    },
    []
  );

  if (isLoading) {
    return (
      <div className={cn('p-4 text-center text-muted-foreground text-sm', className)}>
        {t('加载中...')}
      </div>
    );
  }

  if (error) {
    return (
      <div className={cn('p-4 text-center', className)}>
        <p className="text-sm text-destructive mb-2">{error}</p>
        <Button variant="link" size="sm" onClick={loadProjects}>
          {t('重试')}
        </Button>
      </div>
    );
  }

  if (projects.length === 0) {
    return (
      <div className={cn('p-4 text-center text-muted-foreground', className)}>
        <Film className="h-12 w-12 mx-auto mb-2 opacity-30" />
        <p className="text-sm">{t('暂无视频项目')}</p>
        <p className="text-xs mt-1">{t('点击 + 创建新项目')}</p>
      </div>
    );
  }

  return (
    <div className={cn('p-2 space-y-1', className)}>
      {projects.map((project) => (
        <ProjectItem
          key={project.id}
          project={project}
          isSelected={project.id === selected}
          onSelect={() => onSelect(project)}
          onDelete={() => handleDelete(project.id)}
        />
      ))}
    </div>
  );
}
