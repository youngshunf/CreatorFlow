/**
 * VideoProperties - Properties panel for video configuration
 *
 * Features:
 * - Project name and description editing
 * - Video configuration (dimensions, fps, duration)
 * - Composition props editing
 * - Asset management
 *
 * @requirements 9.4
 */

import * as React from "react";
import { useState, useCallback } from "react";
import {
  Settings,
  Image,
  Music,
  Type,
  Palette,
  Clock,
  Maximize,
  Plus,
  Trash2,
  Video,
  FileText,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Separator } from "@/components/ui/separator";
import { useT } from "@/context/LocaleContext";
import type {
  VideoProject,
  Scene,
  Transition,
  Asset,
  AssetType,
} from "@sprouty-ai/video";
import { ScenePropertiesEditor } from "./ScenePropertiesEditor";
import { TransitionEditor } from "./TransitionEditor";

export interface VideoPropertiesProps {
  /** Current video project */
  project: VideoProject | null;
  /** Callback when project is updated */
  onUpdate: (updates: Partial<VideoProject>) => void;
  /** Callback to open render dialog */
  onRender: () => void;
  /** 当前选中的场景 */
  selectedScene?: Scene | null;
  /** 场景更新回调 */
  onSceneUpdate?: (updates: Partial<Scene>) => void;
  /** 当前选中的过渡 */
  selectedTransition?: Transition | null;
  /** 当前选中的过渡索引 */
  selectedTransitionIdx?: number | null;
  /** 过渡更新回调 */
  onTransitionUpdate?: (index: number, updates: Partial<Transition>) => void;
  /** 添加素材回调 */
  onAddAsset?: (asset: Asset) => void;
  /** 删除素材回调 */
  onRemoveAsset?: (assetId: string) => void;
  /** Optional class name */
  className?: string;
}

/**
 * Section header component
 */
function SectionHeader({
  icon: Icon,
  title,
}: {
  icon: React.ElementType;
  title: string;
}) {
  return (
    <div className="flex items-center gap-2 text-sm font-medium text-foreground mb-2">
      <Icon className="h-4 w-4 text-muted-foreground" />
      {title}
    </div>
  );
}

/**
 * Property row component
 */
function PropertyRow({
  label,
  children,
}: {
  label: string;
  children: React.ReactNode;
}) {
  return (
    <div className="flex items-center justify-between gap-2 py-1">
      <Label className="text-xs text-muted-foreground shrink-0">{label}</Label>
      <div className="flex-1 max-w-[140px]">{children}</div>
    </div>
  );
}

/**
 * 素材类型图标
 */
function AssetIcon({ type }: { type: string }) {
  switch (type) {
    case "image":
      return <Image className="h-3 w-3 text-blue-500 shrink-0" />;
    case "video":
      return <Video className="h-3 w-3 text-purple-500 shrink-0" />;
    case "audio":
      return <Music className="h-3 w-3 text-green-500 shrink-0" />;
    case "font":
      return <Type className="h-3 w-3 text-orange-500 shrink-0" />;
    default:
      return <FileText className="h-3 w-3 text-muted-foreground shrink-0" />;
  }
}

/**
 * VideoProperties component
 */
export function VideoProperties({
  project,
  onUpdate,
  onRender,
  selectedScene,
  onSceneUpdate,
  selectedTransition,
  selectedTransitionIdx,
  onTransitionUpdate,
  onAddAsset,
  onRemoveAsset,
  className,
}: VideoPropertiesProps) {
  const t = useT();
  const [editingName, setEditingName] = useState(false);
  const [nameValue, setNameValue] = useState("");

  // Handle name edit start
  const handleNameEditStart = useCallback(() => {
    if (project) {
      setNameValue(project.name);
      setEditingName(true);
    }
  }, [project]);

  // Handle name edit submit
  const handleNameEditSubmit = useCallback(() => {
    if (nameValue.trim() && nameValue !== project?.name) {
      onUpdate({ name: nameValue.trim() });
    }
    setEditingName(false);
  }, [nameValue, project?.name, onUpdate]);

  // Handle config change
  const handleConfigChange = useCallback(
    (key: keyof VideoProject["config"], value: number) => {
      if (!project) return;
      onUpdate({
        config: {
          ...project.config,
          [key]: value,
        },
      });
    },
    [project, onUpdate],
  );

  // Handle description change
  const handleDescriptionChange = useCallback(
    (e: React.ChangeEvent<HTMLTextAreaElement>) => {
      onUpdate({ description: e.target.value });
    },
    [onUpdate],
  );

  // Preset dimensions
  const handlePresetChange = useCallback(
    (preset: string) => {
      if (!project) return;
      const presets: Record<string, { width: number; height: number }> = {
        "1080p": { width: 1920, height: 1080 },
        "720p": { width: 1280, height: 720 },
        "4k": { width: 3840, height: 2160 },
        vertical: { width: 1080, height: 1920 },
        square: { width: 1080, height: 1080 },
        portrait: { width: 1080, height: 1350 },
      };
      const dimensions = presets[preset];
      if (dimensions) {
        onUpdate({
          config: {
            ...project.config,
            ...dimensions,
          },
        });
      }
    },
    [project, onUpdate],
  );

  if (!project) {
    return (
      <div
        className={cn(
          "p-4 text-center text-muted-foreground text-sm",
          className,
        )}
      >
        {t("选择一个项目查看属性")}
      </div>
    );
  }

  return (
    <div className={cn("p-3 space-y-4", className)}>
      {/* Project Info */}
      <div>
        <SectionHeader icon={Settings} title={t("项目信息")} />
        <div className="space-y-2">
          {/* Name */}
          <div>
            <Label className="text-xs text-muted-foreground">{t("名称")}</Label>
            {editingName ? (
              <Input
                value={nameValue}
                onChange={(e) => setNameValue(e.target.value)}
                onBlur={handleNameEditSubmit}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleNameEditSubmit();
                  if (e.key === "Escape") setEditingName(false);
                }}
                className="h-8 mt-1"
                autoFocus
              />
            ) : (
              <div
                className="text-sm font-medium cursor-pointer hover:bg-muted/50 rounded px-2 py-1 mt-1"
                onClick={handleNameEditStart}
              >
                {project.name}
              </div>
            )}
          </div>

          {/* Description */}
          <div>
            <Label className="text-xs text-muted-foreground">{t("描述")}</Label>
            <Textarea
              value={project.description || ""}
              onChange={handleDescriptionChange}
              placeholder={t("添加项目描述...")}
              className="h-16 mt-1 text-sm resize-none"
            />
          </div>
        </div>
      </div>

      <Separator />

      {/* Video Config */}
      <div>
        <SectionHeader icon={Maximize} title={t("视频配置")} />
        <div className="space-y-1">
          {/* Preset */}
          <PropertyRow label={t("预设")}>
            <Select onValueChange={handlePresetChange}>
              <SelectTrigger className="h-7 text-xs">
                <SelectValue placeholder={t("选择预设")} />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="1080p">1080p (16:9)</SelectItem>
                <SelectItem value="720p">720p (16:9)</SelectItem>
                <SelectItem value="4k">4K (16:9)</SelectItem>
                <SelectItem value="vertical">{t("竖屏")} (9:16)</SelectItem>
                <SelectItem value="square">{t("方形")} (1:1)</SelectItem>
                <SelectItem value="portrait">{t("肖像")} (4:5)</SelectItem>
              </SelectContent>
            </Select>
          </PropertyRow>

          {/* Width */}
          <PropertyRow label={t("宽度")}>
            <Input
              type="number"
              value={project.config.width}
              onChange={(e) =>
                handleConfigChange("width", parseInt(e.target.value) || 1920)
              }
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* Height */}
          <PropertyRow label={t("高度")}>
            <Input
              type="number"
              value={project.config.height}
              onChange={(e) =>
                handleConfigChange("height", parseInt(e.target.value) || 1080)
              }
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* FPS */}
          <PropertyRow label={t("帧率")}>
            <Select
              value={project.config.fps.toString()}
              onValueChange={(v) => handleConfigChange("fps", parseInt(v))}
            >
              <SelectTrigger className="h-7 text-xs">
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="24">24 fps</SelectItem>
                <SelectItem value="25">25 fps</SelectItem>
                <SelectItem value="30">30 fps</SelectItem>
                <SelectItem value="50">50 fps</SelectItem>
                <SelectItem value="60">60 fps</SelectItem>
              </SelectContent>
            </Select>
          </PropertyRow>

          {/* Duration */}
          <PropertyRow label={t("时长(帧)")}>
            <Input
              type="number"
              value={project.config.durationInFrames}
              onChange={(e) =>
                handleConfigChange(
                  "durationInFrames",
                  parseInt(e.target.value) || 300,
                )
              }
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* Duration in seconds (calculated) */}
          <PropertyRow label={t("时长(秒)")}>
            <div className="text-xs text-muted-foreground text-right">
              {(project.config.durationInFrames / project.config.fps).toFixed(
                2,
              )}
              s
            </div>
          </PropertyRow>
        </div>
      </div>

      <Separator />

      {/* 场景属性编辑器 */}
      {selectedScene && onSceneUpdate && (
        <>
          <ScenePropertiesEditor
            scene={selectedScene}
            fps={project.config.fps}
            onUpdate={onSceneUpdate}
          />
          <Separator />
        </>
      )}

      {/* 过渡效果编辑器 */}
      {selectedTransition &&
        selectedTransitionIdx !== null &&
        selectedTransitionIdx !== undefined &&
        onTransitionUpdate && (
          <>
            <TransitionEditor
              transition={selectedTransition}
              onUpdate={(updates) =>
                onTransitionUpdate(selectedTransitionIdx, updates)
              }
            />
            <Separator />
          </>
        )}

      {/* Assets */}
      <div>
        <div className="flex items-center justify-between mb-2">
          <SectionHeader icon={Image} title={t("素材")} />
          {onAddAsset && (
            <Button
              variant="ghost"
              size="icon"
              className="h-6 w-6"
              onClick={async () => {
                if (!project || !window.electronAPI?.video?.selectAssetFiles)
                  return;
                try {
                  const files =
                    await window.electronAPI.video.selectAssetFiles();
                  for (const { filePath, assetType } of files) {
                    const asset = await window.electronAPI.video.addAsset(
                      project.id,
                      filePath,
                      assetType,
                    );
                    onAddAsset(asset);
                  }
                } catch (err) {
                  console.error("添加素材失败:", err);
                }
              }}
              title={t("添加素材")}
            >
              <Plus className="h-3.5 w-3.5" />
            </Button>
          )}
        </div>
        <div className="space-y-1">
          {project.assets.length === 0 ? (
            <div className="text-xs text-muted-foreground text-center py-3">
              {t("暂无素材")}
              {onAddAsset && (
                <p className="mt-1 text-[10px]">{t("点击 + 添加素材文件")}</p>
              )}
            </div>
          ) : (
            <>
              <div className="text-xs text-muted-foreground mb-1">
                {t("{{count}} 个素材", { count: project.assets.length })}
              </div>
              {project.assets.map((asset) => (
                <div
                  key={asset.id}
                  className="group flex items-center gap-2 text-xs py-1 px-2 rounded hover:bg-muted/50"
                >
                  <AssetIcon type={asset.type} />
                  <span className="truncate flex-1">{asset.name}</span>
                  {onRemoveAsset && (
                    <Button
                      variant="ghost"
                      size="icon"
                      className="h-5 w-5 opacity-0 group-hover:opacity-100 transition-opacity text-destructive hover:text-destructive shrink-0"
                      onClick={() => onRemoveAsset(asset.id)}
                    >
                      <Trash2 className="h-3 w-3" />
                    </Button>
                  )}
                </div>
              ))}
            </>
          )}
        </div>
      </div>

      <Separator />

      {/* Render Button */}
      <Button onClick={onRender} className="w-full">
        {t("导出视频")}
      </Button>
    </div>
  );
}
