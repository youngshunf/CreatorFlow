/**
 * TitleAnimation Composition
 *
 * A video composition that displays a title and subtitle with fade-in animation effects.
 * Uses Remotion's interpolate and spring APIs for smooth animations.
 *
 * @requirements 3.1, 3.5, 3.6
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';
import type { TemplateProps } from '../templates/types';

/**
 * Props for TitleAnimation composition
 */
export interface TitleAnimationProps extends TemplateProps {
  /** Main title text */
  title?: string;
  /** Subtitle text */
  subtitle?: string;
  /** Animation style */
  animationStyle?: 'fade' | 'slide' | 'scale' | 'spring';
  /** Title font size in pixels */
  titleFontSize?: number;
  /** Subtitle font size in pixels */
  subtitleFontSize?: number;
}

/**
 * Default colors for the composition
 */
const DEFAULT_COLORS = {
  primary: '#6366f1',
  secondary: '#8b5cf6',
  background: '#1a1a2e',
  text: '#ffffff',
};

/**
 * TitleAnimation - Displays animated title and subtitle
 *
 * This composition creates a professional title card with customizable
 * fade-in animations for both the title and subtitle.
 */
export const TitleAnimation: React.FC<TitleAnimationProps> = ({
  title = 'Welcome',
  subtitle = 'Your subtitle here',
  colors = DEFAULT_COLORS,
  animationStyle = 'spring',
  titleFontSize = 72,
  subtitleFontSize = 36,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Calculate animation progress based on style
  const getTitleOpacity = () => {
    if (animationStyle === 'spring') {
      return spring({
        frame,
        fps,
        config: {
          damping: 200,
          stiffness: 100,
          mass: 0.5,
        },
      });
    }
    return interpolate(frame, [0, 30], [0, 1], {
      extrapolateRight: 'clamp',
    });
  };

  const getSubtitleOpacity = () => {
    if (animationStyle === 'spring') {
      return spring({
        frame: frame - 15, // Delay subtitle animation
        fps,
        config: {
          damping: 200,
          stiffness: 100,
          mass: 0.5,
        },
      });
    }
    return interpolate(frame, [15, 45], [0, 1], {
      extrapolateRight: 'clamp',
    });
  };

  // Calculate slide/scale transforms
  const getTitleTransform = () => {
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

  const getSubtitleTransform = () => {
    if (animationStyle === 'slide') {
      const translateY = interpolate(frame, [15, 45], [30, 0], {
        extrapolateRight: 'clamp',
      });
      return `translateY(${translateY}px)`;
    }
    if (animationStyle === 'scale') {
      const scale = interpolate(frame, [15, 45], [0.8, 1], {
        extrapolateRight: 'clamp',
      });
      return `scale(${scale})`;
    }
    if (animationStyle === 'spring') {
      const translateY = spring({
        frame: frame - 15,
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

  // Fade out animation near the end
  const fadeOutOpacity = interpolate(
    frame,
    [durationInFrames - 30, durationInFrames],
    [1, 0],
    {
      extrapolateLeft: 'clamp',
      extrapolateRight: 'clamp',
    }
  );

  const titleOpacity = Math.min(getTitleOpacity(), fadeOutOpacity);
  const subtitleOpacity = Math.min(
    Math.max(0, getSubtitleOpacity()),
    fadeOutOpacity
  );

  return (
    <AbsoluteFill
      style={{
        backgroundColor: colors.background,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
      }}
    >
      {/* Optional Logo */}
      {logo && (
        <img
          src={logo}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 60,
            left: 60,
            height: 60,
            opacity: titleOpacity,
          }}
        />
      )}

      {/* Title Container */}
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          gap: 20,
        }}
      >
        {/* Main Title */}
        <h1
          style={{
            fontSize: titleFontSize,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleOpacity,
            transform: getTitleTransform(),
            textAlign: 'center',
            maxWidth: '80%',
          }}
        >
          {title}
        </h1>

        {/* Subtitle */}
        {subtitle && (
          <p
            style={{
              fontSize: subtitleFontSize,
              color: colors.primary,
              margin: 0,
              opacity: subtitleOpacity,
              transform: getSubtitleTransform(),
              textAlign: 'center',
              maxWidth: '70%',
            }}
          >
            {subtitle}
          </p>
        )}
      </div>
    </AbsoluteFill>
  );
};

export default TitleAnimation;
