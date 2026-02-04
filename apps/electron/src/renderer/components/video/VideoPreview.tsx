/**
 * VideoPreview - Video preview component with Remotion Player
 *
 * Features:
 * - Integrates @remotion/player for video preview
 * - Displays composition preview
 * - Handles frame synchronization
 *
 * @requirements 9.2
 */

import * as React from 'react';
import { useMemo } from 'react';
import { Player, type PlayerRef } from '@remotion/player';
import { Film } from 'lucide-react';
import { cn } from '@/lib/utils';
import { useT } from '@/context/LocaleContext';
import type { VideoProject } from '@creator-flow/video';
import { RemotionRoot } from '@creator-flow/video';

export interface VideoPreviewProps {
  /** Current video project */
  project: VideoProject | null;
  /** Player ref for external control */
  playerRef: React.RefObject<PlayerRef>;
  /** Callback when frame changes */
  onFrameChange?: (frame: number) => void;
  /** Optional class name */
  className?: string;
}

/**
 * Empty state component
 */
function EmptyState() {
  const t = useT();
  return (
    <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
      <Film className="h-16 w-16 mb-4 opacity-30" />
      <p className="text-sm">{t('选择一个项目开始预览')}</p>
    </div>
  );
}

/**
 * VideoPreview component
 */
export function VideoPreview({
  project,
  playerRef,
  onFrameChange,
  className,
}: VideoPreviewProps) {
  const t = useT();

  // Get the first composition or default
  const composition = useMemo(() => {
    if (!project || project.compositions.length === 0) {
      return null;
    }
    return project.compositions[0];
  }, [project]);

  // Calculate player dimensions to fit container while maintaining aspect ratio
  const playerStyle = useMemo(() => {
    if (!project) return {};
    const { width, height } = project.config;
    return {
      aspectRatio: `${width} / ${height}`,
    };
  }, [project]);

  if (!project) {
    return (
      <div className={cn('h-full bg-muted/30', className)}>
        <EmptyState />
      </div>
    );
  }

  return (
    <div className={cn('h-full flex items-center justify-center bg-black/90 p-4', className)}>
      <div
        className="relative w-full h-full flex items-center justify-center"
        style={{ maxWidth: '100%', maxHeight: '100%' }}
      >
        <div
          className="relative bg-black rounded-lg overflow-hidden shadow-2xl"
          style={{
            ...playerStyle,
            width: '100%',
            maxHeight: '100%',
          }}
        >
          <Player
            ref={playerRef}
            component={RemotionRoot}
            inputProps={{
              compositionId: composition?.id || 'TitleAnimation',
              ...composition?.props,
            }}
            durationInFrames={project.config.durationInFrames}
            fps={project.config.fps}
            compositionWidth={project.config.width}
            compositionHeight={project.config.height}
            style={{
              width: '100%',
              height: '100%',
            }}
            controls
            autoPlay={false}
            loop
            clickToPlay
            doubleClickToFullscreen
            spaceKeyToPlayOrPause
          />
        </div>
      </div>
    </div>
  );
}
