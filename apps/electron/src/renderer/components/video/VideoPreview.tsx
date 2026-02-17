/**
 * VideoPreview - 使用 SceneComposer + @remotion/player 的视频预览组件
 *
 * 通过 SceneComposer 编排多个场景片段，支持过渡效果。
 * Player 的 inputProps = { scenes, transitions }
 *
 * @requirements 9.2
 */

import * as React from "react";
import { useMemo, Component } from "react";
import { Player, type PlayerRef, type RenderLoading } from "@remotion/player";
import { Film, AlertTriangle } from "lucide-react";
import { cn } from "@/lib/utils";
import { useT } from "@/context/LocaleContext";
import type { VideoProject } from "@sprouty-ai/video";
import { SceneComposer } from "@sprouty-ai/video";

/**
 * 错误边界 - 捕获 Player/SceneComposer 渲染错误
 */
interface ErrorBoundaryState {
  hasError: boolean;
  error: Error | null;
}

class PlayerErrorBoundary extends Component<
  { children: React.ReactNode; resetKey: string },
  ErrorBoundaryState
> {
  state: ErrorBoundaryState = { hasError: false, error: null };

  static getDerivedStateFromError(error: Error): ErrorBoundaryState {
    return { hasError: true, error };
  }

  componentDidUpdate(prevProps: { resetKey: string }) {
    if (prevProps.resetKey !== this.props.resetKey && this.state.hasError) {
      this.setState({ hasError: false, error: null });
    }
  }

  render() {
    if (this.state.hasError) {
      return (
        <div className="flex flex-col items-center justify-center h-full text-destructive bg-destructive/5 p-4">
          <AlertTriangle className="h-10 w-10 mb-3 opacity-60" />
          <p className="text-sm font-medium mb-1">渲染出错</p>
          <p className="text-xs text-muted-foreground text-center max-w-[300px]">
            {this.state.error?.message || "未知错误"}
          </p>
        </div>
      );
    }
    return this.props.children;
  }
}

/**
 * Player 加载中占位
 */
const renderLoading: RenderLoading = () => (
  <div className="flex items-center justify-center h-full bg-black/80 text-muted-foreground">
    <div className="text-sm animate-pulse">加载预览...</div>
  </div>
);

/**
 * 计算项目总帧数
 * 总帧数 = sum(scenes.durationInFrames) - sum(有效transitions.durationInFrames)
 */
export function calculateTotalDuration(project: VideoProject): number {
  const scenes = project.scenes ?? [];
  const transitions = project.transitions ?? [];
  const scenesTotal = scenes.reduce((sum, s) => sum + s.durationInFrames, 0);
  const transitionsTotal = transitions
    .filter((t) => t.type !== "none")
    .reduce((sum, t) => sum + t.durationInFrames, 0);
  return Math.max(1, scenesTotal - transitionsTotal);
}

export interface VideoPreviewProps {
  /** 当前视频项目 */
  project: VideoProject | null;
  /** Player ref，用于外部控制 */
  playerRef: React.RefObject<PlayerRef>;
  /** 帧变化回调 */
  onFrameChange?: (frame: number) => void;
  /** 可选 class name */
  className?: string;
}

/**
 * 空状态组件
 */
function EmptyState() {
  const t = useT();
  return (
    <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
      <Film className="h-16 w-16 mb-4 opacity-30" />
      <p className="text-sm">{t("选择一个项目开始预览")}</p>
    </div>
  );
}

/**
 * VideoPreview 组件
 */
export function VideoPreview({
  project,
  playerRef,
  onFrameChange,
  className,
}: VideoPreviewProps) {
  // 计算 SceneComposer 的 inputProps
  const composerProps = useMemo(() => {
    if (!project || !project.scenes || project.scenes.length === 0) return null;
    return {
      scenes: project.scenes,
      transitions: project.transitions ?? [],
    };
  }, [project?.scenes, project?.transitions]);

  // 计算总帧数
  const durationInFrames = useMemo(() => {
    if (!project || !project.scenes || project.scenes.length === 0) return 1;
    return calculateTotalDuration(project);
  }, [project?.scenes, project?.transitions]);

  // 计算播放器宽高比
  const playerStyle = useMemo(() => {
    if (!project) return {};
    const { width, height } = project.config;
    return {
      aspectRatio: `${width} / ${height}`,
    };
  }, [project]);

  if (!project) {
    return (
      <div className={cn("h-full bg-muted/30", className)}>
        <EmptyState />
      </div>
    );
  }

  if (!composerProps) {
    return (
      <div className={cn("h-full bg-muted/30", className)}>
        <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
          <Film className="h-16 w-16 mb-4 opacity-30" />
          <p className="text-sm">暂无场景</p>
        </div>
      </div>
    );
  }

  return (
    <div
      className={cn(
        "h-full flex items-center justify-center bg-black/90 p-4",
        className,
      )}
    >
      <div
        className="relative w-full h-full flex items-center justify-center"
        style={{ maxWidth: "100%", maxHeight: "100%" }}
      >
        <div
          className="relative bg-black rounded-lg overflow-hidden shadow-2xl"
          style={{
            ...playerStyle,
            width: "100%",
            maxHeight: "100%",
          }}
        >
          <PlayerErrorBoundary
            resetKey={`${project.id}-${project.scenes.length}-${project.updatedAt}`}
          >
            <Player
              key={`${project.id}-${project.scenes.length}-${project.updatedAt}`}
              ref={playerRef}
              component={SceneComposer}
              inputProps={composerProps}
              durationInFrames={durationInFrames}
              fps={project.config.fps}
              compositionWidth={project.config.width}
              compositionHeight={project.config.height}
              style={{
                width: "100%",
                height: "100%",
              }}
              acknowledgeRemotionLicense
              renderLoading={renderLoading}
              controls
              autoPlay={false}
              loop
              clickToPlay
              doubleClickToFullscreen
              spaceKeyToPlayOrPause
            />
          </PlayerErrorBoundary>
        </div>
      </div>
    </div>
  );
}
