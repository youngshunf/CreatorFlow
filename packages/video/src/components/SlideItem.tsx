/**
 * SlideItem Component
 *
 * Displays a single slide with image/gradient background and optional title overlay.
 * Designed to be used within a Sequence for slideshow animations.
 */
import React from 'react';
import { useCurrentFrame, useVideoConfig, spring, Img, staticFile } from 'remotion';

export type SlideItemProps = {
  /** Slide data */
  slide: {
    title: string;
    description?: string;
    image?: string;
  };
  /** Color scheme */
  colors: {
    primary: string;
    secondary: string;
    text: string;
  };
  /** Show title overlay */
  showTitle?: boolean;
};

/**
 * SlideItem - Animated slide component
 *
 * This component handles the animation of a single slide.
 * It automatically fades in at the start and fades out at the end.
 */
export const SlideItem: React.FC<SlideItemProps> = ({
  slide,
  colors,
  showTitle = true,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Title animation
  const titleOpacity = spring({
    frame,
    fps,
    config: {
      damping: 200,
      stiffness: 100,
      mass: 0.5,
    },
  });

  const titleTranslateY = spring({
    frame,
    fps,
    config: {
      damping: 200,
      stiffness: 100,
      mass: 0.5,
    },
    from: 30,
    to: 0,
  });

  return (
    <div
      style={{
        width: '100%',
        height: '100%',
        position: 'relative',
        overflow: 'hidden',
      }}
    >
      {/* Background Image */}
      {slide.image ? (
        <Img
          src={staticFile(slide.image)}
          style={{
            width: '100%',
            height: '100%',
            objectFit: 'cover',
          }}
        />
      ) : (
        <div
          style={{
            width: '100%',
            height: '100%',
            background: `linear-gradient(135deg, ${colors.primary} 0%, ${colors.secondary} 100%)`,
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
          }}
        >
          <span
            style={{
              fontSize: 120,
              color: 'rgba(255, 255, 255, 0.2)',
              fontWeight: 'bold',
            }}
          >
            {slide.title}
          </span>
        </div>
      )}

      {/* Title Overlay */}
      {showTitle && (
        <div
          style={{
            position: 'absolute',
            bottom: 0,
            left: 0,
            right: 0,
            padding: '60px',
            background: 'linear-gradient(transparent, rgba(0, 0, 0, 0.7))',
          }}
        >
          <h2
            style={{
              fontSize: 48,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 10,
              opacity: titleOpacity,
              transform: `translateY(${titleTranslateY}px)`,
            }}
          >
            {slide.title}
          </h2>
          {slide.description && (
            <p
              style={{
                fontSize: 24,
                color: 'rgba(255, 255, 255, 0.8)',
                margin: 0,
                opacity: titleOpacity,
                transform: `translateY(${titleTranslateY}px)`,
              }}
            >
              {slide.description}
            </p>
          )}
        </div>
      )}
    </div>
  );
};
