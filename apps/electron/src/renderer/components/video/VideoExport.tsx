/**
 * VideoExport - 视频导出组件
 *
 * 基于 SceneComposer 架构，导出整个项目（scenes + transitions）。
 * 不再需要选择单个 composition，固定使用 SceneComposer 渲染。
 *
 * @requirements 9.7
 */

import * as React from "react";
import { useState, useCallback, useEffect } from "react";
import { Download, Loader2, CheckCircle, XCircle } from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Progress } from "@/components/ui/progress";
import { Separator } from "@/components/ui/separator";
import { useT } from "@/context/LocaleContext";
import type {
  VideoProject,
  OutputFormat,
  QualityPreset,
  RenderProgress,
} from "@sprouty-ai/video";
import { calculateTotalDuration } from "./VideoPreview";

export interface VideoExportProps {
  /** 当前视频项目 */
  project: VideoProject;
  /** 导出回调 */
  onExport: (options: {
    compositionId: string;
    outputFormat: OutputFormat;
    quality: QualityPreset;
  }) => void;
  /** 取消/关闭回调 */
  onCancel: () => void;
  /** 可选 class name */
  className?: string;
}

/**
 * 格式描述
 */
const FORMAT_INFO: Record<OutputFormat, { name: string; description: string }> =
  {
    mp4: { name: "MP4", description: "通用格式，兼容性最好" },
    webm: { name: "WebM", description: "网页优化，文件较小" },
    gif: { name: "GIF", description: "动图格式，无声音" },
  };

/**
 * 质量描述
 */
const QUALITY_INFO: Record<
  QualityPreset,
  { name: string; description: string }
> = {
  draft: { name: "草稿", description: "快速预览，低质量" },
  standard: { name: "标准", description: "平衡质量和大小" },
  high: { name: "高质量", description: "最佳质量，文件较大" },
};

/**
 * VideoExport 组件
 */
export function VideoExport({
  project,
  onExport,
  onCancel,
  className,
}: VideoExportProps) {
  const t = useT();

  // 导出选项
  const [outputFormat, setOutputFormat] = useState<OutputFormat>("mp4");
  const [quality, setQuality] = useState<QualityPreset>("standard");

  // 渲染进度
  const [isRendering, setIsRendering] = useState(false);
  const [renderProgress, setRenderProgress] = useState<RenderProgress | null>(
    null,
  );

  // 计算总时长
  const scenes = project.scenes ?? [];
  const totalDuration = calculateTotalDuration(project);
  const durationSeconds = totalDuration / project.config.fps;

  // 订阅渲染进度
  useEffect(() => {
    if (!window.electronAPI?.video?.onRenderProgress) return;
    const unsubscribe = window.electronAPI.video.onRenderProgress(
      (progress: RenderProgress) => {
        setRenderProgress(progress);
        if (progress.status === "completed" || progress.status === "failed") {
          setIsRendering(false);
        }
      },
    );
    return () => unsubscribe();
  }, []);

  // 处理导出 - 固定使用 SceneComposer
  const handleExport = useCallback(() => {
    if (scenes.length === 0) return;
    setIsRendering(true);
    setRenderProgress({ status: "bundling", progress: 0 });
    onExport({ compositionId: "SceneComposer", outputFormat, quality });
  }, [scenes.length, outputFormat, quality, onExport]);

  // 取消渲染
  const handleCancelRender = useCallback(() => {
    if (window.electronAPI?.video?.cancelRender) {
      window.electronAPI.video.cancelRender();
    }
    setIsRendering(false);
    setRenderProgress(null);
  }, []);

  // 状态文本
  const getStatusText = (status: RenderProgress["status"]): string => {
    switch (status) {
      case "bundling":
        return t("打包中...");
      case "preparing":
        return t("准备中...");
      case "rendering":
        return t("渲染中...");
      case "completed":
        return t("完成");
      case "failed":
        return t("失败");
      default:
        return "";
    }
  };

  return (
    <div className={cn("p-3 space-y-4", className)}>
      {/* 标题 */}
      <div className="flex items-center gap-2">
        <Download className="h-5 w-5 text-primary" />
        <h3 className="font-medium">{t("导出视频")}</h3>
      </div>

      {/* 场景概览 */}
      <div className="text-xs text-muted-foreground">
        {scenes.length === 0
          ? t("暂无场景，请先添加场景")
          : t("{{count}} 个场景", { count: scenes.length })}
      </div>

      <Separator />

      {/* 格式选择 */}
      <div>
        <Label className="text-xs text-muted-foreground mb-1.5 block">
          {t("输出格式")}
        </Label>
        <div className="grid grid-cols-3 gap-2">
          {(Object.keys(FORMAT_INFO) as OutputFormat[]).map((format) => (
            <button
              key={format}
              disabled={isRendering}
              className={cn(
                "p-2 rounded-lg border text-center transition-colors",
                "disabled:opacity-50 disabled:cursor-not-allowed",
                outputFormat === format
                  ? "border-primary bg-primary/10"
                  : "border-border hover:border-primary/50",
              )}
              onClick={() => setOutputFormat(format)}
            >
              <div className="font-medium text-sm">
                {FORMAT_INFO[format].name}
              </div>
              <div className="text-[10px] text-muted-foreground mt-0.5">
                {t(FORMAT_INFO[format].description)}
              </div>
            </button>
          ))}
        </div>
      </div>

      {/* 质量选择 */}
      <div>
        <Label className="text-xs text-muted-foreground mb-1.5 block">
          {t("质量预设")}
        </Label>
        <div className="grid grid-cols-3 gap-2">
          {(Object.keys(QUALITY_INFO) as QualityPreset[]).map((q) => (
            <button
              key={q}
              disabled={isRendering}
              className={cn(
                "p-2 rounded-lg border text-center transition-colors",
                "disabled:opacity-50 disabled:cursor-not-allowed",
                quality === q
                  ? "border-primary bg-primary/10"
                  : "border-border hover:border-primary/50",
              )}
              onClick={() => setQuality(q)}
            >
              <div className="font-medium text-sm">
                {t(QUALITY_INFO[q].name)}
              </div>
              <div className="text-[10px] text-muted-foreground mt-0.5">
                {t(QUALITY_INFO[q].description)}
              </div>
            </button>
          ))}
        </div>
      </div>

      <Separator />

      {/* 导出信息 */}
      <div className="text-xs text-muted-foreground space-y-1">
        <div className="flex justify-between">
          <span>{t("分辨率")}</span>
          <span>
            {project.config.width}x{project.config.height}
          </span>
        </div>
        <div className="flex justify-between">
          <span>{t("帧率")}</span>
          <span>{project.config.fps} fps</span>
        </div>
        <div className="flex justify-between">
          <span>{t("时长")}</span>
          <span>
            {durationSeconds.toFixed(1)}s ({totalDuration} {t("帧")})
          </span>
        </div>
        <div className="flex justify-between">
          <span>{t("场景数")}</span>
          <span>{scenes.length}</span>
        </div>
      </div>

      {/* 渲染进度 */}
      {renderProgress && (
        <div className="space-y-2">
          <div className="flex items-center gap-2">
            {renderProgress.status === "completed" ? (
              <CheckCircle className="h-4 w-4 text-green-500" />
            ) : renderProgress.status === "failed" ? (
              <XCircle className="h-4 w-4 text-destructive" />
            ) : (
              <Loader2 className="h-4 w-4 animate-spin" />
            )}
            <span className="text-sm">
              {getStatusText(renderProgress.status)}
            </span>
            <span className="text-sm text-muted-foreground ml-auto">
              {renderProgress.progress.toFixed(0)}%
            </span>
          </div>
          <Progress value={renderProgress.progress} className="h-2" />
          {renderProgress.error && (
            <p className="text-xs text-destructive">{renderProgress.error}</p>
          )}
        </div>
      )}

      {/* 操作按钮 */}
      <div className="flex gap-2">
        <Button variant="outline" className="flex-1" onClick={onCancel}>
          {t("取消")}
        </Button>
        {isRendering ? (
          <Button
            variant="destructive"
            className="flex-1"
            onClick={handleCancelRender}
          >
            {t("停止渲染")}
          </Button>
        ) : (
          <Button
            className="flex-1"
            onClick={handleExport}
            disabled={scenes.length === 0}
          >
            <Download className="h-4 w-4 mr-1" />
            {t("开始导出")}
          </Button>
        )}
      </div>
    </div>
  );
}
