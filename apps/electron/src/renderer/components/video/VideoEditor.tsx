/**
 * VideoEditor - 视频编辑器主组件（三栏布局）
 *
 * 基于 SceneComposer 架构：
 * - 左侧面板：场景列表（添加/删除/重排）+ 模板
 * - 中央面板：SceneComposer 视频预览 + 时间轴
 * - 右侧面板：属性编辑 + 导出
 *
 * @requirements 9.1
 */

import * as React from "react";
import { useState, useCallback, useRef, useEffect } from "react";
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
import type {
  VideoProject,
  VideoTemplate,
  Scene,
  Transition,
  Asset,
} from "@sprouty-ai/video";
import { VideoPreview, calculateTotalDuration } from "./VideoPreview";
import { VideoTimeline } from "./VideoTimeline";
import { VideoProperties } from "./VideoProperties";
import { VideoProjectList } from "./VideoProjectList";
import { VideoTemplates } from "./VideoTemplates";
import { VideoExport } from "./VideoExport";
import { CreateVideoProjectDialog } from "./CreateVideoProjectDialog";
import { TransitionIndicator } from "./TransitionEditor";

/** 防抖保存延迟（毫秒） */
const SAVE_DEBOUNCE_MS = 300;

export interface VideoEditorProps {
  /** 工作区 ID */
  workspaceId: string;
  /** 可选 class name */
  className?: string;
}

/**
 * 生成唯一 ID
 */
function generateId(): string {
  return `scene-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;
}

/**
 * 从模板创建默认场景列表
 */
function createScenesFromTemplate(template: VideoTemplate): Scene[] {
  return [
    {
      id: generateId(),
      name: template.name,
      compositionId: template.compositionId,
      durationInFrames: template.defaultConfig.durationInFrames,
      props: template.defaultProps || {},
    },
  ];
}

/**
 * VideoEditor 组件
 */
export function VideoEditor({ workspaceId, className }: VideoEditorProps) {
  const t = useT();
  const playerRef = useRef<PlayerRef>(null);

  // 防抖保存定时器
  const saveTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // 状态
  const [currentProject, setCurrentProject] = useState<VideoProject | null>(
    null,
  );
  const [localProjects, setLocalProjects] = useState<VideoProject[]>([]);
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

  // 清理防抖定时器
  useEffect(() => {
    return () => {
      if (saveTimerRef.current) clearTimeout(saveTimerRef.current);
    };
  }, []);

  /**
   * 防抖保存项目到磁盘（通过 IPC）
   * 如果 electronAPI 不可用则静默跳过
   */
  const debouncedSave = useCallback((project: VideoProject) => {
    if (!window.electronAPI?.video?.updateProject) return;
    if (saveTimerRef.current) clearTimeout(saveTimerRef.current);
    saveTimerRef.current = setTimeout(() => {
      window.electronAPI.video
        .updateProject(project.id, project)
        .catch((err: unknown) => {
          console.warn("自动保存失败:", err);
        });
    }, SAVE_DEBOUNCE_MS);
  }, []);

  // 选择项目时切换到场景面板
  const handleProjectSelect = useCallback((project: VideoProject) => {
    setCurrentProject(project);
    setCurrentFrame(0);
    setIsPlaying(false);
    setSelectedSceneId(project.scenes[0]?.id ?? null);
    setLeftTab("scenes");
  }, []);

  // 模板预览
  const handleTemplateSelect = useCallback((template: VideoTemplate) => {
    const now = new Date().toISOString();
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
      scenes: createScenesFromTemplate(template),
      transitions: [],
      assets: [],
      renders: [],
    };
    setCurrentProject(project);
    setSelectedSceneId(project.scenes[0]?.id ?? null);
    setLeftTab("scenes");
  }, []);

  // 从模板创建项目
  const handleTemplateCreate = useCallback((template: VideoTemplate) => {
    setPendingTemplate(template);
    setCreateDialogOpen(true);
  }, []);

  // 确认创建项目（通过 IPC 持久化 + 本地状态）
  const handleCreateFromTemplate = useCallback(
    async (data: { name: string; description: string }) => {
      if (!pendingTemplate) return;
      const now = new Date().toISOString();
      const scenes = createScenesFromTemplate(pendingTemplate);
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
        scenes,
        transitions: [],
        assets: [],
        renders: [],
      };

      // 尝试通过 IPC 创建持久化项目
      if (window.electronAPI?.video?.createProject) {
        try {
          const created = await window.electronAPI.video.createProject({
            name: data.name,
            workspaceId,
            template: pendingTemplate.id,
            config: pendingTemplate.defaultConfig,
            description: data.description || undefined,
          });
          setLocalProjects((prev) => [...prev, created]);
          setCurrentProject(created);
          setSelectedSceneId(created.scenes[0]?.id ?? null);
        } catch (err) {
          console.warn("IPC 创建项目失败，使用本地模式:", err);
          setLocalProjects((prev) => [...prev, project]);
          setCurrentProject(project);
          setSelectedSceneId(scenes[0]?.id ?? null);
        }
      } else {
        setLocalProjects((prev) => [...prev, project]);
        setCurrentProject(project);
        setSelectedSceneId(scenes[0]?.id ?? null);
      }

      setLeftTab("scenes");
      setPendingTemplate(null);
    },
    [pendingTemplate, workspaceId],
  );

  // 更新项目（乐观更新 UI + 防抖持久化）
  const handleProjectUpdate = useCallback(
    (updates: Partial<VideoProject>) => {
      if (!currentProject) return;
      const updated = {
        ...currentProject,
        ...updates,
        updatedAt: new Date().toISOString(),
      };
      // 同步 config.durationInFrames 为实际总帧数
      if (updates.scenes || updates.transitions) {
        updated.config = {
          ...updated.config,
          durationInFrames: calculateTotalDuration(updated),
        };
      }
      setCurrentProject(updated);
      debouncedSave(updated);
    },
    [currentProject, debouncedSave],
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

  // 导出渲染（不再需要 compositionId，固定使用 SceneComposer）
  const handleRender = useCallback(
    async (options: {
      compositionId: string;
      outputFormat: "mp4" | "webm" | "gif";
      quality: "draft" | "standard" | "high";
    }) => {
      if (!currentProject) return;
      try {
        if (!window.electronAPI?.video?.render) {
          toast.info(t("视频导出功能需要完整的 Electron 环境"));
          return;
        }
        const outputPath = await window.electronAPI.video.render({
          projectId: currentProject.id,
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
    [currentProject, t],
  );

  // 新建空项目（通过 IPC 持久化 + 降级到本地模式）
  const handleCreateProject = useCallback(async () => {
    const now = new Date().toISOString();
    const defaultScene: Scene = {
      id: generateId(),
      name: t("标题动画"),
      compositionId: "TitleAnimation",
      durationInFrames: 150,
      props: {
        title: t("欢迎"),
        subtitle: t("在此输入副标题"),
        colors: {
          primary: "#6366f1",
          secondary: "#8b5cf6",
          background: "#1a1a2e",
          text: "#ffffff",
        },
        animationStyle: "spring",
      },
    };
    const localProject: VideoProject = {
      id: `proj-${Date.now()}`,
      name: t("新视频项目"),
      createdAt: now,
      updatedAt: now,
      config: {
        width: 1920,
        height: 1080,
        fps: 30,
        durationInFrames: 150,
      },
      scenes: [defaultScene],
      transitions: [],
      assets: [],
      renders: [],
    };

    // 尝试通过 IPC 创建持久化项目
    if (window.electronAPI?.video?.createProject) {
      try {
        const created = await window.electronAPI.video.createProject({
          name: t("新视频项目"),
          workspaceId,
        });
        // IPC 创建的项目可能没有默认场景，补充一个
        if (created.scenes.length === 0) {
          const updated = {
            ...created,
            scenes: [defaultScene],
            updatedAt: new Date().toISOString(),
          };
          setCurrentProject(updated);
          debouncedSave(updated);
        } else {
          setCurrentProject(created);
        }
        setSelectedSceneId(created.scenes[0]?.id ?? defaultScene.id);
      } catch (err) {
        console.warn("IPC 创建项目失败，使用本地模式:", err);
        setCurrentProject(localProject);
        setSelectedSceneId(defaultScene.id);
      }
    } else {
      setCurrentProject(localProject);
      setSelectedSceneId(defaultScene.id);
    }

    setLeftTab("scenes");
  }, [t, workspaceId, debouncedSave]);

  // ========== 场景属性编辑 ==========

  // 更新单个场景属性
  const handleSceneUpdate = useCallback(
    (updates: Partial<Scene>) => {
      if (!currentProject || !selectedSceneId) return;
      const newScenes = currentProject.scenes.map((s) =>
        s.id === selectedSceneId ? { ...s, ...updates } : s,
      );
      handleProjectUpdate({ scenes: newScenes });
    },
    [currentProject, selectedSceneId, handleProjectUpdate],
  );

  // 更新过渡效果
  const handleTransitionUpdate = useCallback(
    (index: number, updates: Partial<Transition>) => {
      if (!currentProject) return;
      const newTransitions = currentProject.transitions.map((t, i) =>
        i === index ? { ...t, ...updates } : t,
      );
      handleProjectUpdate({ transitions: newTransitions });
    },
    [currentProject, handleProjectUpdate],
  );

  // ========== 素材管理 ==========

  // 添加素材（乐观更新 UI）
  const handleAddAsset = useCallback(
    (asset: Asset) => {
      if (!currentProject) return;
      handleProjectUpdate({
        assets: [...currentProject.assets, asset],
      });
    },
    [currentProject, handleProjectUpdate],
  );

  // 删除素材
  const handleRemoveAsset = useCallback(
    async (assetId: string) => {
      if (!currentProject) return;
      // 乐观更新 UI
      handleProjectUpdate({
        assets: currentProject.assets.filter((a) => a.id !== assetId),
      });
      // 异步删除文件
      if (window.electronAPI?.video?.removeAsset) {
        try {
          await window.electronAPI.video.removeAsset(
            currentProject.id,
            assetId,
          );
        } catch (err) {
          console.warn("删除素材文件失败:", err);
        }
      }
    },
    [currentProject, handleProjectUpdate],
  );

  // 当前选中的场景
  const selectedScene =
    currentProject?.scenes.find((s) => s.id === selectedSceneId) ?? null;

  // 当前选中的过渡
  const selectedTransition =
    selectedTransitionIdx !== null
      ? (currentProject?.transitions[selectedTransitionIdx] ?? null)
      : null;

  // ========== 场景操作 ==========

  // 添加场景
  const handleAddScene = useCallback(() => {
    if (!currentProject) return;
    const newScene: Scene = {
      id: generateId(),
      name: t("新场景"),
      compositionId: "TitleAnimation",
      durationInFrames: 90,
      props: {
        title: t("新场景"),
        subtitle: "",
        colors: {
          primary: "#6366f1",
          secondary: "#8b5cf6",
          background: "#1a1a2e",
          text: "#ffffff",
        },
        animationStyle: "fade",
      },
    };
    const newScenes = [...currentProject.scenes, newScene];
    // 如果已有场景，添加默认过渡
    const newTransitions =
      currentProject.scenes.length > 0
        ? [
            ...currentProject.transitions,
            { type: "fade" as const, durationInFrames: 15 },
          ]
        : [...currentProject.transitions];
    handleProjectUpdate({ scenes: newScenes, transitions: newTransitions });
    setSelectedSceneId(newScene.id);
  }, [currentProject, handleProjectUpdate, t]);

  // 删除场景
  const handleDeleteScene = useCallback(
    (sceneId: string) => {
      if (!currentProject) return;
      const idx = currentProject.scenes.findIndex((s) => s.id === sceneId);
      if (idx === -1) return;
      const newScenes = currentProject.scenes.filter((s) => s.id !== sceneId);
      // 删除对应的过渡：如果删除的不是最后一个场景，删除该场景后面的过渡；
      // 如果是最后一个，删除前面的过渡
      const newTransitions = [...currentProject.transitions];
      if (newTransitions.length > 0) {
        const transIdx =
          idx < newTransitions.length ? idx : newTransitions.length - 1;
        newTransitions.splice(transIdx, 1);
      }
      handleProjectUpdate({ scenes: newScenes, transitions: newTransitions });
      // 选中相邻场景
      if (selectedSceneId === sceneId) {
        const nextScene = newScenes[Math.min(idx, newScenes.length - 1)];
        setSelectedSceneId(nextScene?.id ?? null);
      }
    },
    [currentProject, handleProjectUpdate, selectedSceneId],
  );

  // 场景上移
  const handleMoveSceneUp = useCallback(
    (sceneId: string) => {
      if (!currentProject) return;
      const idx = currentProject.scenes.findIndex((s) => s.id === sceneId);
      if (idx <= 0) return;
      const newScenes = [...currentProject.scenes];
      [newScenes[idx - 1], newScenes[idx]] = [
        newScenes[idx],
        newScenes[idx - 1],
      ];
      // 同步交换过渡
      const newTransitions = [...currentProject.transitions];
      if (idx - 1 < newTransitions.length && idx < newTransitions.length) {
        [newTransitions[idx - 1], newTransitions[idx]] = [
          newTransitions[idx],
          newTransitions[idx - 1],
        ];
      }
      handleProjectUpdate({ scenes: newScenes, transitions: newTransitions });
    },
    [currentProject, handleProjectUpdate],
  );

  // 场景下移
  const handleMoveSceneDown = useCallback(
    (sceneId: string) => {
      if (!currentProject) return;
      const idx = currentProject.scenes.findIndex((s) => s.id === sceneId);
      if (idx === -1 || idx >= currentProject.scenes.length - 1) return;
      const newScenes = [...currentProject.scenes];
      [newScenes[idx], newScenes[idx + 1]] = [
        newScenes[idx + 1],
        newScenes[idx],
      ];
      // 同步交换过渡
      const newTransitions = [...currentProject.transitions];
      if (idx < newTransitions.length && idx + 1 < newTransitions.length) {
        [newTransitions[idx], newTransitions[idx + 1]] = [
          newTransitions[idx + 1],
          newTransitions[idx],
        ];
      }
      handleProjectUpdate({ scenes: newScenes, transitions: newTransitions });
    },
    [currentProject, handleProjectUpdate],
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
    player.addEventListener("seeked", handleSeek as EventListener);

    return () => {
      player.removeEventListener("play", handlePlay);
      player.removeEventListener("pause", handlePause);
      player.removeEventListener("seeked", handleSeek as EventListener);
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
              selected={currentProject?.id}
              onSelect={handleProjectSelect}
              extraProjects={localProjects}
            />
          ) : leftTab === "scenes" ? (
            <SceneList
              scenes={currentProject?.scenes ?? []}
              transitions={currentProject?.transitions ?? []}
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
              fps={currentProject?.config.fps ?? 30}
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
            project={currentProject}
            playerRef={playerRef}
            onFrameChange={handleFrameChange}
          />
        </div>

        {/* 时间轴 */}
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
              project={currentProject}
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
          ) : currentProject ? (
            <VideoExport
              project={currentProject}
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
