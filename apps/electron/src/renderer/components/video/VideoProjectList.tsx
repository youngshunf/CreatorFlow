/**
 * VideoProjectList - 项目列表组件
 *
 * Phase 5 适配：从 DB 类型 VideoProject (snake_case) 读取，
 * 通过适配层转为旧 VideoProject 传给 onSelect 回调。
 */

import * as React from "react";
import { useState, useEffect, useCallback, useMemo } from "react";
import { Film, Trash2, MoreVertical, Calendar } from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import { useT } from "@/context/LocaleContext";
import type { VideoProject as DbVideoProject } from "@sprouty-ai/shared/db/types";
import type { VideoProject } from "@sprouty-ai/video";

export interface VideoProjectListProps {
  /** 工作区 ID */
  workspaceId: string;
  /** 当前选中的项目 ID */
  selected?: string;
  /** 选中项目回调（传旧类型，调用方通过 ID 重新加载完整数据） */
  onSelect: (project: VideoProject) => void;
  /** 可选 class name */
  className?: string;
}

/** 格式化日期 */
function formatDate(dateString: string): string {
  const date = new Date(dateString);
  return date.toLocaleDateString(undefined, {
    month: "short",
    day: "numeric",
  });
}

/** DB VideoProject → 旧 VideoProject（轻量转换，仅用于列表展示和 onSelect） */
function dbToLegacy(p: DbVideoProject): VideoProject {
  return {
    id: p.id,
    name: p.name,
    description: p.description ?? undefined,
    createdAt: p.created_at,
    updatedAt: p.updated_at,
    config: {
      width: p.width,
      height: p.height,
      fps: p.fps,
      durationInFrames: 0, // 列表不需要精确值，实际编辑时会从 Full 重新计算
    },
    scenes: [],
    transitions: [],
    assets: [],
    renders: [],
  };
}

/** 项目列表项 */
function ProjectItem({
  project,
  isSelected,
  onSelect,
  onDelete,
}: {
  project: DbVideoProject;
  isSelected: boolean;
  onSelect: () => void;
  onDelete: () => void;
}) {
  const t = useT();

  return (
    <div
      className={cn(
        "group flex items-start gap-3 p-3 rounded-lg cursor-pointer transition-colors",
        "hover:bg-muted/50",
        isSelected && "bg-primary/10 hover:bg-primary/15",
      )}
      onClick={onSelect}
    >
      {/* 图标 */}
      <div
        className={cn(
          "shrink-0 w-10 h-10 rounded-md flex items-center justify-center",
          isSelected ? "bg-primary text-primary-foreground" : "bg-muted",
        )}
      >
        <Film className="h-5 w-5" />
      </div>

      {/* 信息 */}
      <div className="flex-1 min-w-0">
        <div className="font-medium text-sm truncate">{project.name}</div>
        <div className="flex items-center gap-3 mt-1 text-xs text-muted-foreground">
          <span className="flex items-center gap-1">
            <Calendar className="h-3 w-3" />
            {formatDate(project.updated_at)}
          </span>
          <span className="text-[10px]">
            {project.width}x{project.height} · {project.fps}fps
          </span>
        </div>
        {project.description && (
          <div className="text-xs text-muted-foreground mt-1 truncate">
            {project.description}
          </div>
        )}
      </div>

      {/* 操作 */}
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
            {t("删除")}
          </DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
    </div>
  );
}

/**
 * VideoProjectList 组件
 */
export function VideoProjectList({
  workspaceId,
  selected,
  onSelect,
  className,
}: VideoProjectListProps) {
  const t = useT();
  const [projects, setProjects] = useState<DbVideoProject[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // 按 updated_at 排序
  const sortedProjects = useMemo(
    () =>
      [...projects].sort(
        (a, b) =>
          new Date(b.updated_at).getTime() - new Date(a.updated_at).getTime(),
      ),
    [projects],
  );

  // 加载项目列表
  const loadProjects = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    try {
      if (!window.electronAPI?.video?.listProjects) {
        setProjects([]);
        return;
      }
      const list = await window.electronAPI.video.listProjects(workspaceId, "");
      setProjects(list);
    } catch (err) {
      console.error("加载项目列表失败:", err);
      setProjects([]);
    } finally {
      setIsLoading(false);
    }
  }, [workspaceId]);

  // 初始加载
  useEffect(() => {
    loadProjects();
  }, [loadProjects]);

  // 删除项目
  const handleDelete = useCallback(
    async (projectId: string) => {
      try {
        if (!window.electronAPI?.video?.deleteProject) return;
        await window.electronAPI.video.deleteProject(workspaceId, projectId);
        setProjects((prev) => prev.filter((p) => p.id !== projectId));
      } catch (err) {
        console.error("删除项目失败:", err);
      }
    },
    [workspaceId],
  );

  if (isLoading) {
    return (
      <div
        className={cn(
          "p-4 text-center text-muted-foreground text-sm",
          className,
        )}
      >
        {t("加载中...")}
      </div>
    );
  }

  if (error) {
    return (
      <div className={cn("p-4 text-center", className)}>
        <p className="text-sm text-destructive mb-2">{error}</p>
        <Button variant="link" size="sm" onClick={loadProjects}>
          {t("重试")}
        </Button>
      </div>
    );
  }

  if (sortedProjects.length === 0) {
    return (
      <div className={cn("p-4 text-center text-muted-foreground", className)}>
        <Film className="h-12 w-12 mx-auto mb-2 opacity-30" />
        <p className="text-sm">{t("暂无视频项目")}</p>
        <p className="text-xs mt-1">{t("点击 + 创建新项目")}</p>
      </div>
    );
  }

  return (
    <div className={cn("p-2 space-y-1", className)}>
      {sortedProjects.map((project) => (
        <ProjectItem
          key={project.id}
          project={project}
          isSelected={project.id === selected}
          onSelect={() => onSelect(dbToLegacy(project))}
          onDelete={() => handleDelete(project.id)}
        />
      ))}
    </div>
  );
}
