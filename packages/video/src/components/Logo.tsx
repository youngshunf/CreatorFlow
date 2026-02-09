/**
 * Logo Component
 *
 * A reusable component for displaying logos with animation effects.
 * Supports various entrance animations and positioning options.
 *
 * @requirements 4.4, 4.5
 */
import React from 'react';
import {
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Img,
} from 'remotion';

/**
 * Logo animation type options
 */
export type LogoAnimation =
  | 'none'
  | 'fade'
  | 'scale'
  | 'slideUp'
  | 'slideDown'
  | 'slideLeft'
  | 'slideRight'
  | 'spring'
  | 'bounce'
  | 'rotate'
  | 'flip';

/**
 * Logo position options
 */
export type LogoPosition =
  | 'top-left'
  | 'top-center'
  | 'top-right'
  | 'center-left'
  | 'center'
  | 'center-right'
  | 'bottom-left'
  | 'bottom-center'
  | 'bottom-right';

/**
 * Props for Logo component
 */
export interface LogoProps {
  /** Logo image URL */
  src: string;
  /** Logo width in pixels (height auto-calculated to maintain aspect ratio) */
  width?: number;
  /** Logo height in pixels (if specified, overrides auto-calculation) */
  height?: number;
  /** Animation type */
  animation?: LogoAnimation;
  /** Delay before animation starts (in frames) */
  delay?: number;
  /** Duration of the animation (in frames) */
  duration?: number;
  /** Position of the logo */
  position?: LogoPosition;
  /** Padding from edges (in pixels) */
  padding?: number;
  /** Opacity (0-1) */
  opacity?: number;
  /** Enable exit animation at the end */
  exitAnimation?: boolean;
  /** Frames before end to start exit animation */
  exitOffset?: number;
  /** Spring configuration for spring-based animations */
  springConfig?: {
    damping?: number;
    stiffness?: number;
    mass?: number;
  };
  /** Additional CSS styles */
  style?: React.CSSProperties;
  /** Alt text for accessibility */
  alt?: string;
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
 * Get position styles based on position prop
 */
const getPositionStyles = (
  position: LogoPosition,
  padding: number
): React.CSSProperties => {
  const baseStyles: React.CSSProperties = {
    position: 'absolute',
  };

  switch (position) {
    case 'top-left':
      return { ...baseStyles, top: padding, left: padding };
    case 'top-center':
      return { ...baseStyles, top: padding, left: '50%', transform: 'translateX(-50%)' };
    case 'top-right':
      return { ...baseStyles, top: padding, right: padding };
    case 'center-left':
      return { ...baseStyles, top: '50%', left: padding, transform: 'translateY(-50%)' };
    case 'center':
      return { ...baseStyles, top: '50%', left: '50%', transform: 'translate(-50%, -50%)' };
    case 'center-right':
      return { ...baseStyles, top: '50%', right: padding, transform: 'translateY(-50%)' };
    case 'bottom-left':
      return { ...baseStyles, bottom: padding, left: padding };
    case 'bottom-center':
      return { ...baseStyles, bottom: padding, left: '50%', transform: 'translateX(-50%)' };
    case 'bottom-right':
      return { ...baseStyles, bottom: padding, right: padding };
    default:
      return { ...baseStyles, top: padding, left: padding };
  }
};

/**
 * Logo - Displays a logo with customizable animation effects
 *
 * This component provides various animation options for logo elements,
 * making it easy to create professional logo animations in video compositions.
 *
 * @example
 * ```tsx
 * <Logo
 *   src="/logo.png"
 *   width={120}
 *   animation="spring"
 *   position="top-left"
 *   padding={40}
 * />
 * ```
 */
export const Logo: React.FC<LogoProps> = ({
  src,
  width = 100,
  height,
  animation = 'fade',
  delay = 0,
  duration = 30,
  position = 'top-left',
  padding = 40,
  opacity: maxOpacity = 1,
  exitAnimation = false,
  exitOffset = 30,
  springConfig = DEFAULT_SPRING_CONFIG,
  style = {},
  alt = 'Logo',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Adjusted frame accounting for delay
  const adjustedFrame = Math.max(0, frame - delay);

  // Exit animation timing
  const exitStartFrame = durationInFrames - exitOffset;
  const isExiting = exitAnimation && frame >= exitStartFrame;
  const exitProgress = isExiting
    ? interpolate(frame, [exitStartFrame, durationInFrames], [0, 1], {
        extrapolateRight: 'clamp',
      })
    : 0;

  /**
   * Calculate entrance opacity
   */
  const getEntranceOpacity = (): number => {
    if (frame < delay) return 0;

    if (animation === 'spring' || animation === 'bounce') {
      return spring({
        frame: adjustedFrame,
        fps,
        config: springConfig,
      });
    }

    if (animation === 'none') return 1;

    return interpolate(adjustedFrame, [0, duration], [0, 1], {
      extrapolateRight: 'clamp',
    });
  };

  /**
   * Calculate final opacity including exit
   */
  const getOpacity = (): number => {
    const entranceOpacity = getEntranceOpacity();

    if (isExiting) {
      return entranceOpacity * (1 - exitProgress) * maxOpacity;
    }

    return entranceOpacity * maxOpacity;
  };

  /**
   * Calculate transform based on animation type
   */
  const getTransform = (): string => {
    const transforms: string[] = [];

    // Handle entrance animation
    if (frame >= delay) {
      switch (animation) {
        case 'scale': {
          const scale = interpolate(adjustedFrame, [0, duration], [0.5, 1], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`scale(${scale})`);
          break;
        }

        case 'slideUp': {
          const translateY = interpolate(adjustedFrame, [0, duration], [50, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`translateY(${translateY}px)`);
          break;
        }

        case 'slideDown': {
          const translateY = interpolate(adjustedFrame, [0, duration], [-50, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`translateY(${translateY}px)`);
          break;
        }

        case 'slideLeft': {
          const translateX = interpolate(adjustedFrame, [0, duration], [50, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`translateX(${translateX}px)`);
          break;
        }

        case 'slideRight': {
          const translateX = interpolate(adjustedFrame, [0, duration], [-50, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`translateX(${translateX}px)`);
          break;
        }

        case 'spring': {
          const translateY = spring({
            frame: adjustedFrame,
            fps,
            config: springConfig,
            from: 30,
            to: 0,
          });
          transforms.push(`translateY(${translateY}px)`);
          break;
        }

        case 'bounce': {
          const translateY = spring({
            frame: adjustedFrame,
            fps,
            config: {
              damping: 10,
              stiffness: 200,
              mass: 0.5,
            },
            from: -50,
            to: 0,
          });
          transforms.push(`translateY(${translateY}px)`);
          break;
        }

        case 'rotate': {
          const rotate = interpolate(adjustedFrame, [0, duration], [180, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`rotate(${rotate}deg)`);
          break;
        }

        case 'flip': {
          const rotateY = interpolate(adjustedFrame, [0, duration], [90, 0], {
            extrapolateRight: 'clamp',
          });
          transforms.push(`perspective(500px) rotateY(${rotateY}deg)`);
          break;
        }
      }
    } else {
      // Initial state before animation
      switch (animation) {
        case 'scale':
          transforms.push('scale(0.5)');
          break;
        case 'slideUp':
          transforms.push('translateY(50px)');
          break;
        case 'slideDown':
          transforms.push('translateY(-50px)');
          break;
        case 'slideLeft':
          transforms.push('translateX(50px)');
          break;
        case 'slideRight':
          transforms.push('translateX(-50px)');
          break;
        case 'spring':
          transforms.push('translateY(30px)');
          break;
        case 'bounce':
          transforms.push('translateY(-50px)');
          break;
        case 'rotate':
          transforms.push('rotate(180deg)');
          break;
        case 'flip':
          transforms.push('perspective(500px) rotateY(90deg)');
          break;
      }
    }

    // Handle exit animation
    if (isExiting) {
      const exitScale = interpolate(exitProgress, [0, 1], [1, 0.8]);
      const exitTranslateY = interpolate(exitProgress, [0, 1], [0, -20]);
      transforms.push(`scale(${exitScale})`);
      transforms.push(`translateY(${exitTranslateY}px)`);
    }

    return transforms.length > 0 ? transforms.join(' ') : 'none';
  };

  const positionStyles = getPositionStyles(position, padding);
  const opacity = getOpacity();
  const transform = getTransform();

  // Merge position transform with animation transform
  const finalTransform = positionStyles.transform
    ? `${positionStyles.transform} ${transform !== 'none' ? transform : ''}`
    : transform;

  return (
    <div
      style={{
        ...positionStyles,
        transform: finalTransform,
        opacity,
        willChange: 'opacity, transform',
        ...style,
      }}
    >
      <Img
        src={src}
        alt={alt}
        style={{
          width,
          height: height || 'auto',
          display: 'block',
        }}
      />
    </div>
  );
};

/**
 * AnimatedLogo - Convenience component with preset animations
 */
export interface AnimatedLogoProps extends Omit<LogoProps, 'animation'> {
  /** Preset animation style */
  preset?: 'subtle' | 'dynamic' | 'professional' | 'playful';
}

export const AnimatedLogo: React.FC<AnimatedLogoProps> = ({
  preset = 'professional',
  ...props
}) => {
  const presetConfig: Record<string, Partial<LogoProps>> = {
    subtle: {
      animation: 'fade',
      duration: 45,
    },
    dynamic: {
      animation: 'spring',
      duration: 30,
      springConfig: { damping: 15, stiffness: 150, mass: 0.5 },
    },
    professional: {
      animation: 'scale',
      duration: 30,
    },
    playful: {
      animation: 'bounce',
      duration: 45,
    },
  };

  return <Logo {...presetConfig[preset]} {...props} />;
};

export default Logo;
