/**
 * VideoTimeline - Timeline component for video editing
 *
 * Features:
 * - Frame position display and scrubbing
 * - Play/pause controls
 * - Time display (current / total)
 * - Visual timeline with frame markers
 *
 * @requirements 9.3
 */

import * as React from 'react';
import { useCallback, useMemo, useRef } from 'react';
import { Play, Pause, SkipBack, SkipForward, Rewind, FastForward } from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Slider } from '@/components/ui/slider';
import { useT } from '@/context/LocaleContext';
import type { VideoProject } from '@creator-flow/video';

export interface VideoTimelineProps {
  /** Current video project */
  project: VideoProject | null;
  /** Current frame position */
  currentFrame: number;
  /** Callback when frame changes */
  onFrameChange: (frame: number) => void;
  /** Whether video is playing */
  isPlaying: boolean;
  /** Callback to toggle play/pause */
  onPlayPause: () => void;
  /** Optional class name */
  className?: string;
}

/**
 * Format frame number to time string (MM:SS:FF)
 */
function formatTime(frame: number, fps: number): string {
  const totalSeconds = Math.floor(frame / fps);
  const minutes = Math.floor(totalSeconds / 60);
  const seconds = totalSeconds % 60;
  const frames = frame % fps;
  return `${minutes.toString().padStart(2, '0')}:${seconds.toString().padStart(2, '0')}:${frames.toString().padStart(2, '0')}`;
}

/**
 * VideoTimeline component
 */
export function VideoTimeline({
  project,
  currentFrame,
  onFrameChange,
  isPlaying,
  onPlayPause,
  className,
}: VideoTimelineProps) {
  const t = useT();
  const timelineRef = useRef<HTMLDivElement>(null);

  const fps = project?.config.fps || 30;
  const totalFrames = project?.config.durationInFrames || 300;

  // Time display
  const currentTime = useMemo(() => formatTime(currentFrame, fps), [currentFrame, fps]);
  const totalTime = useMemo(() => formatTime(totalFrames, fps), [totalFrames, fps]);

  // Handle slider change
  const handleSliderChange = useCallback(
    (value: number[]) => {
      onFrameChange(value[0]);
    },
    [onFrameChange]
  );

  // Jump to start
  const handleJumpToStart = useCallback(() => {
    onFrameChange(0);
  }, [onFrameChange]);

  // Jump to end
  const handleJumpToEnd = useCallback(() => {
    onFrameChange(totalFrames - 1);
  }, [onFrameChange, totalFrames]);

  // Step backward (1 second)
  const handleStepBackward = useCallback(() => {
    const newFrame = Math.max(0, currentFrame - fps);
    onFrameChange(newFrame);
  }, [currentFrame, fps, onFrameChange]);

  // Step forward (1 second)
  const handleStepForward = useCallback(() => {
    const newFrame = Math.min(totalFrames - 1, currentFrame + fps);
    onFrameChange(newFrame);
  }, [currentFrame, fps, totalFrames, onFrameChange]);

  // Frame step backward
  const handleFrameBackward = useCallback(() => {
    const newFrame = Math.max(0, currentFrame - 1);
    onFrameChange(newFrame);
  }, [currentFrame, onFrameChange]);

  // Frame step forward
  const handleFrameForward = useCallback(() => {
    const newFrame = Math.min(totalFrames - 1, currentFrame + 1);
    onFrameChange(newFrame);
  }, [currentFrame, totalFrames, onFrameChange]);

  if (!project) {
    return (
      <div className={cn('h-full flex items-center justify-center text-muted-foreground', className)}>
        <p className="text-sm">{t('无项目')}</p>
      </div>
    );
  }

  return (
    <div className={cn('h-full flex flex-col p-2', className)}>
      {/* Controls */}
      <div className="flex items-center gap-2 mb-2">
        {/* Transport controls */}
        <div className="flex items-center gap-1">
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleJumpToStart}
            title={t('跳到开头')}
          >
            <SkipBack className="h-4 w-4" />
          </Button>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleStepBackward}
            title={t('后退1秒')}
          >
            <Rewind className="h-4 w-4" />
          </Button>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleFrameBackward}
            title={t('上一帧')}
          >
            <SkipBack className="h-3 w-3" />
          </Button>
          <Button
            variant="default"
            size="icon"
            className="h-9 w-9"
            onClick={onPlayPause}
            title={isPlaying ? t('暂停') : t('播放')}
          >
            {isPlaying ? (
              <Pause className="h-4 w-4" />
            ) : (
              <Play className="h-4 w-4 ml-0.5" />
            )}
          </Button>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleFrameForward}
            title={t('下一帧')}
          >
            <SkipForward className="h-3 w-3" />
          </Button>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleStepForward}
            title={t('前进1秒')}
          >
            <FastForward className="h-4 w-4" />
          </Button>
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={handleJumpToEnd}
            title={t('跳到结尾')}
          >
            <SkipForward className="h-4 w-4" />
          </Button>
        </div>

        {/* Time display */}
        <div className="flex items-center gap-2 ml-auto text-sm font-mono">
          <span className="text-foreground">{currentTime}</span>
          <span className="text-muted-foreground">/</span>
          <span className="text-muted-foreground">{totalTime}</span>
        </div>

        {/* Frame display */}
        <div className="text-xs text-muted-foreground">
          {t('帧')}: {currentFrame} / {totalFrames}
        </div>
      </div>

      {/* Timeline slider */}
      <div ref={timelineRef} className="flex-1 flex items-center">
        <Slider
          value={[currentFrame]}
          min={0}
          max={totalFrames - 1}
          step={1}
          onValueChange={handleSliderChange}
          className="w-full"
        />
      </div>

      {/* Timeline markers */}
      <div className="flex justify-between text-xs text-muted-foreground mt-1">
        <span>0:00</span>
        <span>{formatTime(Math.floor(totalFrames / 4), fps)}</span>
        <span>{formatTime(Math.floor(totalFrames / 2), fps)}</span>
        <span>{formatTime(Math.floor((totalFrames * 3) / 4), fps)}</span>
        <span>{totalTime}</span>
      </div>
    </div>
  );
}
