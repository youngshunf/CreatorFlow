/**
 * TitleText Component
 *
 * Displays animated title text with fade-in and fade-out effects.
 * Designed to be used within a Sequence for time-based animations.
 */
import React from 'react';
import { useCurrentFrame, useVideoConfig, interpolate, spring } from 'remotion';

export type TitleTextProps = {
  /** Title text to display */
  text: string;
  /** Font size in pixels */
  fontSize?: number;
  /** Text color */
  color?: string;
  /** Animation style */
  animationStyle?: 'fade' | 'slide' | 'scale' | 'spring';
};

/**
 * TitleText - Animated title component
 *
 * This component handles the animation of a single title text element.
 * It automatically fades in at the start and fades out at the end.
 */
export const TitleText: React.FC<TitleTextProps> = ({
  text,
  fontSize = 72,
  color = '#ffffff',
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
      const translateY = interpolate(frame, [0, 30], [50, 0], {
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
        from: 30,
        to: 0,
      });
      return `translateY(${translateY}px)`;
    }
    return 'none';
  };

  return (
    <h1
      style={{
        fontSize,
        fontWeight: 'bold',
        color,
        margin: 0,
        opacity: getOpacity(),
        transform: getTransform(),
        textAlign: 'center',
        maxWidth: '80%',
      }}
    >
      {text}
    </h1>
  );
};
