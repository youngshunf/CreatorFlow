/**
 * VideoEditor - 视频编辑器主组件（三栏布局）
 *
 * 基于 SceneComposer 架构：
 * - 左侧面板：场景列表（添加/删除/重排）+ 模板
 * - 中央面板：SceneComposer 视频预览 + 时间轴
 * - 右侧面板：属性编辑 + 导出
 *
 * Phase 5 适配：状态使用 VideoProjectFull (DB)，通过适配层转为旧 VideoProject 传给子组件。
 */

import * as React from "react";
import { useState, useCallback, useRef, useEffect, useMemo } from "react";
import {
  Plus,
  Film,
  Settings,
  Download,
  Trash2,
  ChevronUp,
  ChevronDown,
} from "lucide-react";
import { toast } from "sonner";
import type { PlayerRef } from "@remotion/player";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { ScrollArea } from "@/components/ui/scroll-area";
import { useT } from "@/context/LocaleContext";
import type { VideoProjectFull } from "@sprouty-ai/shared/db/types";
import type {
  VideoProject,
  VideoTemplate,
  Scene,
  Transition,
  Asset,
} from "@sprouty-ai/video";
import {
  toLegacyProject,
  toDbProjectUpdate,
  toDbSceneUpdate,
  toDbTransitionUpdate,
} from "@/lib/video-adapter";
import { VideoPreview } from "./VideoPreview";
import { VideoTimeline } from "./VideoTimeline";
import { VideoProperties } from "./VideoProperties";
import { VideoProjectList } from "./VideoProjectList";
import { VideoTemplates } from "./VideoTemplates";
import { VideoExport } from "./VideoExport";
import { CreateVideoProjectDialog } from "./CreateVideoProjectDialog";
import { TransitionIndicator } from "./TransitionEditor";
import { useCreatorMedia } from "@/pages/creator-media/hooks/useCreatorMedia";

/** 检测是否是临时预览项目（不持久化到数据库） */
function isTemporaryProject(project: VideoProjectFull | null): boolean {
  return project?.id.startsWith('temp-') ?? false
}

export interface VideoEditorProps {
  /** 工作区 ID */
  workspaceId: string;
  /** 可选 class name */
  className?: string;
}

/**
 * VideoEditor 组件
 */
export function VideoEditor({ workspaceId, className }: VideoEditorProps) {
  const t = useT();
  const playerRef = useRef<PlayerRef>(null);
  const { activeProject } = useCreatorMedia();

  // 状态 — 使用 DB 类型
  const [currentProjectFull, setCurrentProjectFull] =
    useState<VideoProjectFull | null>(null);
  const [currentFrame, setCurrentFrame] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [selectedSceneId, setSelectedSceneId] = useState<string | null>(null);
  const [selectedTransitionIdx, setSelectedTransitionIdx] = useState<
    number | null
  >(null);
  const [leftTab, setLeftTab] = useState<"projects" | "scenes" | "templates">(
    "projects",
  );
  const [rightTab, setRightTab] = useState<"properties" | "export">(
    "properties",
  );

  // 创建项目对话框
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [pendingTemplate, setPendingTemplate] = useState<VideoTemplate | null>(
    null,
  );

  // 派生旧类型 — 子组件零改动
  const legacyProject = useMemo(
    () => (currentProjectFull ? toLegacyProject(currentProjectFull) : null),
    [currentProjectFull],
  );

  // 当前选中的场景（从旧类型派生）
  const selectedScene =
    legacyProject?.scenes.find((s) => s.id === selectedSceneId) ?? null;

  // 当前选中的过渡
  const selectedTransition =
    selectedTransitionIdx !== null
      ? (legacyProject?.transitions[selectedTransitionIdx] ?? null)
      : null;

  // 排序后的 DB 场景（用于 transition 索引映射）
  const sortedDbScenes = useMemo(
    () =>
      currentProjectFull
        ? [...currentProjectFull.scenes].sort(
            (a, b) => a.sort_order - b.sort_order,
          )
        : [],
    [currentProjectFull],
  );

  /** 从 IPC 刷新项目数据 */
  const refetchProject = useCallback(
    async (projectId: string) => {
      if (!window.electronAPI?.video?.getProject) return;
      try {
        const full = await window.electronAPI.video.getProject(
          workspaceId,
          projectId,
        );
        if (full) setCurrentProjectFull(full);
      } catch (err) {
        console.warn("刷新项目失败:", err);
      }
    },
    [workspaceId],
  );

  // 选择项目（从 VideoProjectList，接收旧类型但通过 IPC 加载完整数据）
  const handleProjectSelect = useCallback(
    async (project: VideoProject) => {
      if (window.electronAPI?.video?.getProject) {
        try {
          const full = await window.electronAPI.video.getProject(
            workspaceId,
            project.id,
          );
          if (full) {
            setCurrentProjectFull(full);
            setCurrentFrame(0);
            setIsPlaying(false);
            const firstScene = [...full.scenes].sort(
              (a, b) => a.sort_order - b.sort_order,
            )[0];
            setSelectedSceneId(firstScene?.id ?? null);
            setLeftTab("scenes");
            return;
          }
        } catch {
          // 降级处理
        }
      }
    },
    [workspaceId],
  );

  // 模板预览（纯前端内存构建，不持久化）
  const handleTemplateSelect = useCallback(
    (template: VideoTemplate) => {
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
      };
      setCurrentProjectFull(tempProject);
      setSelectedSceneId(tempProject.scenes[0]?.id ?? null);
      setLeftTab("scenes");
    },
    [t],
  );

  // 从模板创建项目
  const handleTemplateCreate = useCallback((template: VideoTemplate) => {
    setPendingTemplate(template);
    setCreateDialogOpen(true);
  }, []);

  // 确认创建项目（通过 IPC 持久化）
  const handleCreateFromTemplate = useCallback(
    async (data: { name: string; description: string }) => {
      if (!pendingTemplate) return;
      if (!window.electronAPI?.video?.createFromTemplate) return;

      try {
        const created = await window.electronAPI.video.createFromTemplate({
          workspaceId,
          projectId: activeProject?.id || '',
          contentId: '',
          templateId: pendingTemplate.id,
          name: data.name,
          description: data.description || undefined,
        });
        setCurrentProjectFull(created);
        const firstScene = [...created.scenes].sort(
          (a, b) => a.sort_order - b.sort_order,
        )[0];
        setSelectedSceneId(firstScene?.id ?? null);
        setLeftTab("scenes");
        setPendingTemplate(null);
      } catch (err) {
        console.error("从模板创建项目失败:", err);
        toast.error(t("创建项目失败"));
      }
    },
    [pendingTemplate, workspaceId, t],
  );

  // 更新项目属性（通过 IPC + refetch）
  const handleProjectUpdate = useCallback(
    async (updates: Partial<VideoProject>) => {
      if (!currentProjectFull) return;

      // 临时预览项目不持久化，只在内存中更新
      if (isTemporaryProject(currentProjectFull)) {
        setCurrentProjectFull(prev => prev ? { ...prev, ...updates } : null);
        return;
      }

      const dbUpdates = toDbProjectUpdate(updates);
      if (Object.keys(dbUpdates).length === 0) return;

      try {
        if (window.electronAPI?.video?.updateProject) {
          await window.electronAPI.video.updateProject(
            workspaceId,
            currentProjectFull.id,
            dbUpdates,
          );
          await refetchProject(currentProjectFull.id);
        }
      } catch (err) {
        console.warn("更新项目失败:", err);
      }
    },
    [currentProjectFull, workspaceId, refetchProject],
  );

  // 帧变化
  const handleFrameChange = useCallback((frame: number) => {
    setCurrentFrame(frame);
    if (playerRef.current) {
      playerRef.current.seekTo(frame);
    }
  }, []);

  // 播放/暂停
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

  // 导出渲染
  const handleRender = useCallback(
    async (options: {
      compositionId: string;
      outputFormat: "mp4" | "webm" | "gif";
      quality: "draft" | "standard" | "high";
    }) => {
      if (!currentProjectFull) return;
      try {
        if (!window.electronAPI?.video?.render) {
          toast.info(t("视频导出功能需要完整的 Electron 环境"));
          return;
        }
        const outputPath = await window.electronAPI.video.render({
          projectId: currentProjectFull.id,
          workspaceId,
          outputFormat: options.outputFormat,
          quality: options.quality,
        });
        if (outputPath) {
          toast.success(t("视频导出成功"), { description: outputPath });
        }
      } catch (error) {
        toast.error(t("视频导出失败"), {
          description: error instanceof Error ? error.message : String(error),
        });
      }
    },
    [currentProjectFull, workspaceId, t],
  );

  // 新建空项目（通过 IPC）
  const handleCreateProject = useCallback(async () => {
    if (!window.electronAPI?.video?.createProject) return;
    try {
      const created = await window.electronAPI.video.createProject({
        name: t("新视频项目"),
        workspaceId,
        projectId: activeProject?.id || '',
        contentId: '',
      });
      setCurrentProjectFull(created);
      const firstScene = [...created.scenes].sort(
        (a, b) => a.sort_order - b.sort_order,
      )[0];
      setSelectedSceneId(firstScene?.id ?? null);
      setLeftTab("scenes");
    } catch (err) {
      console.error("创建项目失败:", err);
      toast.error(t("创建项目失败"));
    }
  }, [t, workspaceId]);

  // ========== 场景属性编辑 ==========

  // 更新单个场景属性（通过 IPC）
  const handleSceneUpdate = useCallback(
    async (updates: Partial<Scene>) => {
      if (!currentProjectFull || !selectedSceneId) return;
      if (!window.electronAPI?.video?.updateScene) return;

      // 临时预览项目不持久化，只在内存中更新
      if (isTemporaryProject(currentProjectFull)) {
        setCurrentProjectFull(prev => {
          if (!prev) return null;
          return {
            ...prev,
            scenes: prev.scenes.map(s =>
              s.id === selectedSceneId
                ? { ...s, ...toDbSceneUpdate(updates) }
                : s
            ),
          };
        });
        return;
      }

      const dbUpdates = toDbSceneUpdate(updates);
      try {
        await window.electronAPI.video.updateScene({
          workspaceId,
          sceneId: selectedSceneId,
          updates: dbUpdates,
        });
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("更新场景失败:", err);
      }
    },
    [currentProjectFull, selectedSceneId, workspaceId, refetchProject],
  );

  // 更新过渡效果（通过 IPC 更新对应场景的 transition 字段）
  const handleTransitionUpdate = useCallback(
    async (index: number, updates: Partial<Transition>) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.updateScene) return;

      // transitions[i] 对应 sortedDbScenes[i+1]
      const targetScene = sortedDbScenes[index + 1];
      if (!targetScene) return;

      // 临时预览项目不持久化，只在内存中更新
      if (isTemporaryProject(currentProjectFull)) {
        setCurrentProjectFull(prev => {
          if (!prev) return null;
          return {
            ...prev,
            scenes: prev.scenes.map(s =>
              s.id === targetScene.id
                ? { ...s, ...toDbTransitionUpdate(updates) }
                : s
            ),
          };
        });
        return;
      }

      const dbUpdates = toDbTransitionUpdate(updates);
      try {
        await window.electronAPI.video.updateScene({
          workspaceId,
          sceneId: targetScene.id,
          updates: dbUpdates,
        });
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("更新过渡失败:", err);
      }
    },
    [currentProjectFull, sortedDbScenes, workspaceId, refetchProject],
  );

  // ========== 素材管理 ==========

  // 添加素材（通过 IPC）
  const handleAddAsset = useCallback(
    async (asset: Asset) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.addAsset) return;

      // 临时预览项目不支持编辑，提示用户先保存
      if (isTemporaryProject(currentProjectFull)) {
        toast.info(t("请先保存项目后再编辑"));
        return;
      }

      try {
        await window.electronAPI.video.addAsset({
          workspaceId,
          projectId: currentProjectFull.id,
          filePath: asset.path,
          type: asset.type,
          name: asset.name,
        });
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("添加素材失败:", err);
      }
    },
    [currentProjectFull, workspaceId, refetchProject, t],
  );

  // 删除素材（通过 IPC）
  const handleRemoveAsset = useCallback(
    async (assetId: string) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.removeAsset) return;

      // 临时预览项目不支持编辑，提示用户先保存
      if (isTemporaryProject(currentProjectFull)) {
        toast.info(t("请先保存项目后再编辑"));
        return;
      }

      try {
        await window.electronAPI.video.removeAsset(
          workspaceId,
          currentProjectFull.id,
          assetId,
        );
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("删除素材失败:", err);
      }
    },
    [currentProjectFull, workspaceId, refetchProject, t],
  );

  // ========== 场景操作 ==========

  // 添加场景（通过 IPC）
  const handleAddScene = useCallback(async () => {
    if (!currentProjectFull) return;
    if (!window.electronAPI?.video?.addScene) return;

    // 临时预览项目不支持编辑，提示用户先保存
    if (isTemporaryProject(currentProjectFull)) {
      toast.info(t("请先保存项目后再编辑"));
      return;
    }

    try {
      const { sceneId } = await window.electronAPI.video.addScene({
        workspaceId,
        projectId: currentProjectFull.id,
        compositionId: "TitleAnimation",
        name: t("新场景"),
        durationInFrames: 90,
        props: JSON.stringify({
          title: t("新场景"),
          subtitle: "",
          colors: {
            primary: "#6366f1",
            secondary: "#8b5cf6",
            background: "#1a1a2e",
            text: "#ffffff",
          },
          animationStyle: "fade",
        }),
        transitionType: currentProjectFull.scenes.length > 0 ? "fade" : "none",
        transitionDuration: 15,
      });
      await refetchProject(currentProjectFull.id);
      setSelectedSceneId(sceneId);
    } catch (err) {
      console.warn("添加场景失败:", err);
      toast.error(t("添加场景失败"));
    }
  }, [currentProjectFull, workspaceId, refetchProject, t]);

  // 删除场景（通过 IPC）
  const handleDeleteScene = useCallback(
    async (sceneId: string) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.removeScene) return;

      // 临时预览项目不支持编辑，提示用户先保存
      if (isTemporaryProject(currentProjectFull)) {
        toast.info(t("请先保存项目后再编辑"));
        return;
      }

      try {
        await window.electronAPI.video.removeScene(
          workspaceId,
          currentProjectFull.id,
          sceneId,
        );
        await refetchProject(currentProjectFull.id);

        // 选中相邻场景
        if (selectedSceneId === sceneId) {
          const remaining = sortedDbScenes.filter((s) => s.id !== sceneId);
          setSelectedSceneId(remaining[0]?.id ?? null);
        }
      } catch (err) {
        console.warn("删除场景失败:", err);
      }
    },
    [
      currentProjectFull,
      workspaceId,
      refetchProject,
      selectedSceneId,
      sortedDbScenes,
      t,
    ],
  );

  // 场景上移
  const handleMoveSceneUp = useCallback(
    async (sceneId: string) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.reorderScenes) return;

      // 临时预览项目不支持编辑，提示用户先保存
      if (isTemporaryProject(currentProjectFull)) {
        toast.info(t("请先保存项目后再编辑"));
        return;
      }

      const idx = sortedDbScenes.findIndex((s) => s.id === sceneId);
      if (idx <= 0) return;

      const newOrder = sortedDbScenes.map((s) => s.id);
      [newOrder[idx - 1], newOrder[idx]] = [newOrder[idx], newOrder[idx - 1]];

      try {
        await window.electronAPI.video.reorderScenes(
          workspaceId,
          currentProjectFull.id,
          newOrder,
        );
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("重排场景失败:", err);
      }
    },
    [currentProjectFull, sortedDbScenes, workspaceId, refetchProject, t],
  );

  // 场景下移
  const handleMoveSceneDown = useCallback(
    async (sceneId: string) => {
      if (!currentProjectFull) return;
      if (!window.electronAPI?.video?.reorderScenes) return;

      // 临时预览项目不支持编辑，提示用户先保存
      if (isTemporaryProject(currentProjectFull)) {
        toast.info(t("请先保存项目后再编辑"));
        return;
      }

      const idx = sortedDbScenes.findIndex((s) => s.id === sceneId);
      if (idx === -1 || idx >= sortedDbScenes.length - 1) return;

      const newOrder = sortedDbScenes.map((s) => s.id);
      [newOrder[idx], newOrder[idx + 1]] = [newOrder[idx + 1], newOrder[idx]];

      try {
        await window.electronAPI.video.reorderScenes(
          workspaceId,
          currentProjectFull.id,
          newOrder,
        );
        await refetchProject(currentProjectFull.id);
      } catch (err) {
        console.warn("重排场景失败:", err);
      }
    },
    [currentProjectFull, sortedDbScenes, workspaceId, refetchProject, t],
  );

  // 同步播放器状态
  useEffect(() => {
    const player = playerRef.current;
    if (!player) return;

    const handlePlay = () => setIsPlaying(true);
    const handlePause = () => setIsPlaying(false);
    const handleSeek = (e: { detail: { frame: number } }) => {
      setCurrentFrame(e.detail.frame);
    };

    player.addEventListener("play", handlePlay);
    player.addEventListener("pause", handlePause);
    player.addEventListener("seeked", handleSeek);

    return () => {
      player.removeEventListener("play", handlePlay);
      player.removeEventListener("pause", handlePause);
      player.removeEventListener("seeked", handleSeek);
    };
  }, []);

  return (
    <div className={cn("flex h-full", className)}>
      {/* 左侧面板 - 项目/场景/模板 */}
      <div className="w-64 border-r flex flex-col shrink-0">
        <div className="p-2 border-b flex items-center justify-between titlebar-no-drag relative z-panel h-[40px]">
          <Tabs
            value={leftTab}
            onValueChange={(v) => setLeftTab(v as typeof leftTab)}
          >
            <TabsList className="h-8">
              <TabsTrigger value="projects" className="text-xs px-2">
                <Film className="h-3.5 w-3.5 mr-1" />
                {t("项目")}
              </TabsTrigger>
              <TabsTrigger value="scenes" className="text-xs px-2">
                {t("场景")}
              </TabsTrigger>
              <TabsTrigger value="templates" className="text-xs px-2">
                {t("模板")}
              </TabsTrigger>
            </TabsList>
          </Tabs>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={
              leftTab === "scenes" ? handleAddScene : handleCreateProject
            }
            title={leftTab === "scenes" ? t("添加场景") : t("新建项目")}
          >
            <Plus className="h-4 w-4" />
          </Button>
        </div>

        <ScrollArea className="flex-1">
          {leftTab === "projects" ? (
            <VideoProjectList
              workspaceId={workspaceId}
              selected={currentProjectFull?.id}
              onSelect={handleProjectSelect}
            />
          ) : leftTab === "scenes" ? (
            <SceneList
              scenes={legacyProject?.scenes ?? []}
              transitions={legacyProject?.transitions ?? []}
              selectedId={selectedSceneId}
              selectedTransitionIdx={selectedTransitionIdx}
              onSelect={(id) => {
                setSelectedSceneId(id);
                setSelectedTransitionIdx(null);
              }}
              onTransitionSelect={(idx) => {
                setSelectedTransitionIdx(idx);
                setSelectedSceneId(null);
              }}
              onDelete={handleDeleteScene}
              onMoveUp={handleMoveSceneUp}
              onMoveDown={handleMoveSceneDown}
              fps={legacyProject?.config.fps ?? 30}
            />
          ) : (
            <VideoTemplates
              onSelect={handleTemplateSelect}
              onCreate={handleTemplateCreate}
            />
          )}
        </ScrollArea>
      </div>

      {/* 中央面板 - 预览 + 时间轴 */}
      <div className="flex-1 flex flex-col min-w-0">
        {/* 预览 */}
        <div className="flex-1 min-h-0">
          <VideoPreview
            project={legacyProject}
            playerRef={playerRef}
            onFrameChange={handleFrameChange}
          />
        </div>

        {/* 时间轴 */}
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

      {/* 右侧面板 - 属性 + 导出 */}
      <div className="w-72 border-l flex flex-col shrink-0">
        <div className="p-2 border-b titlebar-no-drag relative z-panel h-[40px]">
          <Tabs
            value={rightTab}
            onValueChange={(v) => setRightTab(v as "properties" | "export")}
          >
            <TabsList className="h-8 w-full">
              <TabsTrigger value="properties" className="text-xs flex-1">
                <Settings className="h-3.5 w-3.5 mr-1" />
                {t("属性")}
              </TabsTrigger>
              <TabsTrigger value="export" className="text-xs flex-1">
                <Download className="h-3.5 w-3.5 mr-1" />
                {t("导出")}
              </TabsTrigger>
            </TabsList>
          </Tabs>
        </div>

        <ScrollArea className="flex-1">
          {rightTab === "properties" ? (
            <VideoProperties
              project={legacyProject}
              onUpdate={handleProjectUpdate}
              onRender={() => setRightTab("export")}
              selectedScene={selectedScene}
              onSceneUpdate={handleSceneUpdate}
              selectedTransition={selectedTransition}
              selectedTransitionIdx={selectedTransitionIdx}
              onTransitionUpdate={handleTransitionUpdate}
              onAddAsset={handleAddAsset}
              onRemoveAsset={handleRemoveAsset}
            />
          ) : legacyProject ? (
            <VideoExport
              project={legacyProject}
              onExport={handleRender}
              onCancel={() => setRightTab("properties")}
            />
          ) : (
            <div className="p-4 text-center text-muted-foreground text-sm">
              {t("请先选择或创建一个项目")}
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
  );
}

// ============================================================================
// SceneList - 场景列表组件
// ============================================================================

interface SceneListProps {
  scenes: Scene[];
  transitions: Transition[];
  selectedId: string | null;
  selectedTransitionIdx: number | null;
  onSelect: (id: string) => void;
  onTransitionSelect: (idx: number) => void;
  onDelete: (id: string) => void;
  onMoveUp: (id: string) => void;
  onMoveDown: (id: string) => void;
  fps: number;
}

function SceneList({
  scenes,
  transitions,
  selectedId,
  selectedTransitionIdx,
  onSelect,
  onTransitionSelect,
  onDelete,
  onMoveUp,
  onMoveDown,
  fps,
}: SceneListProps) {
  const t = useT();

  if (scenes.length === 0) {
    return (
      <div className="p-4 text-center text-muted-foreground">
        <Film className="h-12 w-12 mx-auto mb-2 opacity-30" />
        <p className="text-sm">{t("暂无场景")}</p>
        <p className="text-xs mt-1">{t("点击 + 添加场景")}</p>
      </div>
    );
  }

  return (
    <div className="p-2 space-y-1">
      {scenes.map((scene, index) => (
        <React.Fragment key={scene.id}>
          <div
            className={cn(
              "group flex items-center gap-2 p-2 rounded-lg cursor-pointer transition-colors",
              selectedId === scene.id
                ? "bg-primary/10 border border-primary/30"
                : "hover:bg-muted/50 border border-transparent",
            )}
            onClick={() => onSelect(scene.id)}
          >
            {/* 序号 */}
            <div className="w-5 h-5 rounded bg-muted flex items-center justify-center text-xs font-mono text-muted-foreground shrink-0">
              {index + 1}
            </div>

            {/* 场景信息 */}
            <div className="flex-1 min-w-0">
              <div className="text-sm font-medium truncate">{scene.name}</div>
              <div className="text-[10px] text-muted-foreground">
                {scene.compositionId} ·{" "}
                {(scene.durationInFrames / fps).toFixed(1)}s
              </div>
            </div>

            {/* 操作按钮 */}
            <div className="flex items-center gap-0.5 opacity-0 group-hover:opacity-100 transition-opacity shrink-0">
              <Button
                variant="ghost"
                size="icon"
                className="h-6 w-6"
                onClick={(e) => {
                  e.stopPropagation();
                  onMoveUp(scene.id);
                }}
                disabled={index === 0}
              >
                <ChevronUp className="h-3 w-3" />
              </Button>
              <Button
                variant="ghost"
                size="icon"
                className="h-6 w-6"
                onClick={(e) => {
                  e.stopPropagation();
                  onMoveDown(scene.id);
                }}
                disabled={index === scenes.length - 1}
              >
                <ChevronDown className="h-3 w-3" />
              </Button>
              <Button
                variant="ghost"
                size="icon"
                className="h-6 w-6 text-destructive hover:text-destructive"
                onClick={(e) => {
                  e.stopPropagation();
                  onDelete(scene.id);
                }}
              >
                <Trash2 className="h-3 w-3" />
              </Button>
            </div>
          </div>

          {/* 过渡效果指示条（场景之间） */}
          {index < scenes.length - 1 && transitions[index] && (
            <TransitionIndicator
              transition={transitions[index]}
              onClick={() => onTransitionSelect(index)}
              isSelected={selectedTransitionIdx === index}
            />
          )}
        </React.Fragment>
      ))}
    </div>
  );
}
