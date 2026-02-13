/**
 * SubtitleText Component
 *
 * Displays animated subtitle text with fade-in and fade-out effects.
 * Designed to be used within a Sequence for time-based animations.
 */
import React from 'react';
import { useCurrentFrame, useVideoConfig, interpolate, spring } from 'remotion';

export type SubtitleTextProps = {
  /** Subtitle text to display */
  text: string;
  /** Font size in pixels */
  fontSize?: number;
  /** Text color */
  color?: string;
  /** Animation style */
  animationStyle?: 'fade' | 'slide' | 'scale' | 'spring';
};

/**
 * SubtitleText - Animated subtitle component
 *
 * This component handles the animation of a single subtitle text element.
 * It automatically fades in at the start and fades out at the end.
 */
export const SubtitleText: React.FC<SubtitleTextProps> = ({
  text,
  fontSize = 36,
  color = '#6366f1',
  animationStyle = 'spring',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Calculate opacity based on animation style
  const getOpacity = () => {
    if (animationStyle === 'spring') {
      const fadeIn = spring({
        frame,
        fps,
        config: {
          damping: 200,
          stiffness: 100,
          mass: 0.5,
        },
      });
      const fadeOut = interpolate(
        frame,
        [durationInFrames - 30, durationInFrames],
        [1, 0],
        {
          extrapolateLeft: 'clamp',
          extrapolateRight: 'clamp',
        }
      );
      return Math.min(fadeIn, fadeOut);
    }

    // Fade animation
    const fadeIn = interpolate(frame, [0, 30], [0, 1], {
      extrapolateRight: 'clamp',
    });
    const fadeOut = interpolate(
      frame,
      [durationInFrames - 30, durationInFrames],
      [1, 0],
      {
        extrapolateLeft: 'clamp',
        extrapolateRight: 'clamp',
      }
    );
    return Math.min(fadeIn, fadeOut);
  };

  // Calculate transform based on animation style
  const getTransform = () => {
    if (animationStyle === 'slide') {
      const translateY = interpolate(frame, [0, 30], [30, 0], {
        extrapolateRight: 'clamp',
      });
      return `translateY(${translateY}px)`;
    }
    if (animationStyle === 'scale') {
      const scale = interpolate(frame, [0, 30], [0.8, 1], {
        extrapolateRight: 'clamp',
      });
      return `scale(${scale})`;
    }
    if (animationStyle === 'spring') {
      const translateY = spring({
        frame,
        fps,
        config: {
          damping: 200,
          stiffness: 100,
          mass: 0.5,
        },
        from: 20,
        to: 0,
      });
      return `translateY(${translateY}px)`;
    }
    return 'none';
  };

  return (
    <p
      style={{
        fontSize,
        color,
        margin: 0,
        opacity: getOpacity(),
        transform: getTransform(),
        textAlign: 'center',
        maxWidth: '70%',
      }}
    >
      {text}
    </p>
  );
};
