/**
 * AnimatedText Component
 *
 * A reusable component for displaying text with various animation effects.
 * Supports fade-in, slide-in, scale, and spring animations.
 *
 * @requirements 4.1, 4.5
 */
import React from 'react';
import {
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

/**
 * Animation type options
 */
export type AnimationType = 'fade' | 'slideUp' | 'slideDown' | 'slideLeft' | 'slideRight' | 'scale' | 'spring' | 'typewriter';

/**
 * Props for AnimatedText component
 */
export interface AnimatedTextProps {
  /** Text content to display */
  text: string;
  /** Animation type */
  animation?: AnimationType;
  /** Delay before animation starts (in frames) */
  delay?: number;
  /** Duration of the animation (in frames) */
  duration?: number;
  /** Font size in pixels */
  fontSize?: number;
  /** Font weight */
  fontWeight?: 'normal' | 'bold' | number;
  /** Text color */
  color?: string;
  /** Text alignment */
  textAlign?: 'left' | 'center' | 'right';
  /** Font family */
  fontFamily?: string;
  /** Line height */
  lineHeight?: number;
  /** Letter spacing */
  letterSpacing?: number;
  /** Text shadow */
  textShadow?: string;
  /** Additional CSS styles */
  style?: React.CSSProperties;
  /** Spring configuration for spring animation */
  springConfig?: {
    damping?: number;
    stiffness?: number;
    mass?: number;
  };
  /** Distance for slide animations (in pixels) */
  slideDistance?: number;
  /** Scale range for scale animation [from, to] */
  scaleRange?: [number, number];
}

/**
 * Default spring configuration
 */
const DEFAULT_SPRING_CONFIG = {
  damping: 200,
  stiffness: 100,
  mass: 0.5,
};

/**
 * AnimatedText - Displays text with customizable animation effects
 *
 * This component provides various animation options for text elements,
 * making it easy to create engaging text animations in video compositions.
 *
 * @example
 * ```tsx
 * <AnimatedText
 *   text="Hello World"
 *   animation="spring"
 *   fontSize={48}
 *   color="#ffffff"
 *   delay={15}
 * />
 * ```
 */
export const AnimatedText: React.FC<AnimatedTextProps> = ({
  text,
  animation = 'fade',
  delay = 0,
  duration = 30,
  fontSize = 24,
  fontWeight = 'normal',
  color = '#ffffff',
  textAlign = 'center',
  fontFamily = 'system-ui, -apple-system, sans-serif',
  lineHeight = 1.4,
  letterSpacing = 0,
  textShadow,
  style = {},
  springConfig = DEFAULT_SPRING_CONFIG,
  slideDistance = 50,
  scaleRange = [0.8, 1],
}) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();

  // Adjusted frame accounting for delay
  const adjustedFrame = Math.max(0, frame - delay);

  /**
   * Calculate opacity based on animation type
   */
  const getOpacity = (): number => {
    if (frame < delay) return 0;

    if (animation === 'spring') {
      return spring({
        frame: adjustedFrame,
        fps,
        config: springConfig,
      });
    }

    if (animation === 'typewriter') {
      // Typewriter doesn't fade, it reveals characters
      return 1;
    }

    return interpolate(adjustedFrame, [0, duration], [0, 1], {
      extrapolateRight: 'clamp',
    });
  };

  /**
   * Calculate transform based on animation type
   */
  const getTransform = (): string => {
    if (frame < delay) {
      // Initial state before animation starts
      switch (animation) {
        case 'slideUp':
          return `translateY(${slideDistance}px)`;
        case 'slideDown':
          return `translateY(-${slideDistance}px)`;
        case 'slideLeft':
          return `translateX(${slideDistance}px)`;
        case 'slideRight':
          return `translateX(-${slideDistance}px)`;
        case 'scale':
          return `scale(${scaleRange[0]})`;
        case 'spring':
          return `translateY(${slideDistance}px)`;
        default:
          return 'none';
      }
    }

    switch (animation) {
      case 'slideUp': {
        const translateY = interpolate(adjustedFrame, [0, duration], [slideDistance, 0], {
          extrapolateRight: 'clamp',
        });
        return `translateY(${translateY}px)`;
      }
      case 'slideDown': {
        const translateY = interpolate(adjustedFrame, [0, duration], [-slideDistance, 0], {
          extrapolateRight: 'clamp',
        });
        return `translateY(${translateY}px)`;
      }
      case 'slideLeft': {
        const translateX = interpolate(adjustedFrame, [0, duration], [slideDistance, 0], {
          extrapolateRight: 'clamp',
        });
        return `translateX(${translateX}px)`;
      }
      case 'slideRight': {
        const translateX = interpolate(adjustedFrame, [0, duration], [-slideDistance, 0], {
          extrapolateRight: 'clamp',
        });
        return `translateX(${translateX}px)`;
      }
      case 'scale': {
        const scale = interpolate(adjustedFrame, [0, duration], scaleRange, {
          extrapolateRight: 'clamp',
        });
        return `scale(${scale})`;
      }
      case 'spring': {
        const translateY = spring({
          frame: adjustedFrame,
          fps,
          config: springConfig,
          from: slideDistance,
          to: 0,
        });
        return `translateY(${translateY}px)`;
      }
      default:
        return 'none';
    }
  };

  /**
   * Render typewriter effect
   */
  const renderTypewriter = (): React.ReactNode => {
    if (frame < delay) return null;

    // Calculate how many characters to show
    const charsPerFrame = text.length / duration;
    const visibleChars = Math.min(
      Math.floor(adjustedFrame * charsPerFrame),
      text.length
    );

    return (
      <span>
        {text.slice(0, visibleChars)}
        {visibleChars < text.length && (
          <span
            style={{
              opacity: adjustedFrame % 15 < 7.5 ? 1 : 0,
              marginLeft: 2,
            }}
          >
            |
          </span>
        )}
      </span>
    );
  };

  const opacity = getOpacity();
  const transform = getTransform();

  return (
    <div
      style={{
        fontSize,
        fontWeight,
        color,
        textAlign,
        fontFamily,
        lineHeight,
        letterSpacing,
        textShadow,
        opacity,
        transform,
        willChange: 'opacity, transform',
        ...style,
      }}
    >
      {animation === 'typewriter' ? renderTypewriter() : text}
    </div>
  );
};

export default AnimatedText;
