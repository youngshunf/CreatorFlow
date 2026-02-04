/**
 * Background Component
 *
 * A reusable component for rendering various background types.
 * Supports solid colors, gradients, images, and animated backgrounds.
 *
 * @requirements 4.3, 4.5
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  Img,
} from 'remotion';

/**
 * Background type options
 */
export type BackgroundType = 'solid' | 'gradient' | 'image' | 'video' | 'animated';

/**
 * Gradient direction options
 */
export type GradientDirection =
  | 'to-right'
  | 'to-left'
  | 'to-top'
  | 'to-bottom'
  | 'to-top-right'
  | 'to-top-left'
  | 'to-bottom-right'
  | 'to-bottom-left'
  | 'radial';

/**
 * Animation type for animated backgrounds
 */
export type BackgroundAnimation =
  | 'none'
  | 'pulse'
  | 'gradientShift'
  | 'zoom'
  | 'pan'
  | 'parallax';

/**
 * Props for Background component
 */
export interface BackgroundProps {
  /** Type of background */
  type?: BackgroundType;
  /** Solid color (for 'solid' type) */
  color?: string;
  /** Gradient colors (for 'gradient' type) */
  gradientColors?: string[];
  /** Gradient direction */
  gradientDirection?: GradientDirection;
  /** Image URL (for 'image' type) */
  imageUrl?: string;
  /** Image fit mode */
  imageFit?: 'cover' | 'contain' | 'fill' | 'none';
  /** Image position */
  imagePosition?: string;
  /** Background animation type */
  animation?: BackgroundAnimation;
  /** Animation speed (frames per cycle) */
  animationSpeed?: number;
  /** Overlay color with opacity */
  overlay?: string;
  /** Blur amount in pixels */
  blur?: number;
  /** Brightness adjustment (0-2, 1 is normal) */
  brightness?: number;
  /** Children to render on top of background */
  children?: React.ReactNode;
  /** Additional CSS styles */
  style?: React.CSSProperties;
}

/**
 * Convert gradient direction to CSS value
 */
const getGradientAngle = (direction: GradientDirection): string => {
  switch (direction) {
    case 'to-right':
      return '90deg';
    case 'to-left':
      return '270deg';
    case 'to-top':
      return '0deg';
    case 'to-bottom':
      return '180deg';
    case 'to-top-right':
      return '45deg';
    case 'to-top-left':
      return '315deg';
    case 'to-bottom-right':
      return '135deg';
    case 'to-bottom-left':
      return '225deg';
    case 'radial':
      return 'radial';
    default:
      return '180deg';
  }
};

/**
 * Background - Renders various background types with optional animations
 *
 * This component provides a flexible way to create backgrounds for video
 * compositions, supporting solid colors, gradients, images, and animations.
 *
 * @example
 * ```tsx
 * // Solid color background
 * <Background type="solid" color="#1a1a2e" />
 *
 * // Gradient background
 * <Background
 *   type="gradient"
 *   gradientColors={['#6366f1', '#8b5cf6', '#a855f7']}
 *   gradientDirection="to-bottom-right"
 * />
 *
 * // Image background with overlay
 * <Background
 *   type="image"
 *   imageUrl="https://example.com/bg.jpg"
 *   overlay="rgba(0, 0, 0, 0.5)"
 * />
 * ```
 */
export const Background: React.FC<BackgroundProps> = ({
  type = 'solid',
  color = '#1a1a2e',
  gradientColors = ['#6366f1', '#8b5cf6'],
  gradientDirection = 'to-bottom',
  imageUrl,
  imageFit = 'cover',
  imagePosition = 'center',
  animation = 'none',
  animationSpeed = 120,
  overlay,
  blur = 0,
  brightness = 1,
  children,
  style = {},
}) => {
  const frame = useCurrentFrame();
  const { fps, width, height } = useVideoConfig();

  /**
   * Get background style based on type
   */
  const getBackgroundStyle = (): React.CSSProperties => {
    switch (type) {
      case 'solid':
        return { backgroundColor: color };

      case 'gradient': {
        const angle = getGradientAngle(gradientDirection);
        const colorStops = gradientColors.join(', ');

        if (gradientDirection === 'radial') {
          return {
            background: `radial-gradient(circle at center, ${colorStops})`,
          };
        }

        return {
          background: `linear-gradient(${angle}, ${colorStops})`,
        };
      }

      case 'image':
        // Image is rendered separately
        return {};

      default:
        return { backgroundColor: color };
    }
  };

  /**
   * Get animation transform
   */
  const getAnimationTransform = (): string => {
    const cycleProgress = (frame % animationSpeed) / animationSpeed;

    switch (animation) {
      case 'pulse': {
        const scale = interpolate(
          Math.sin(cycleProgress * Math.PI * 2),
          [-1, 1],
          [1, 1.05]
        );
        return `scale(${scale})`;
      }

      case 'zoom': {
        const scale = interpolate(frame, [0, animationSpeed * 2], [1, 1.2], {
          extrapolateRight: 'clamp',
        });
        return `scale(${scale})`;
      }

      case 'pan': {
        const translateX = interpolate(
          Math.sin(cycleProgress * Math.PI * 2),
          [-1, 1],
          [-5, 5]
        );
        const translateY = interpolate(
          Math.cos(cycleProgress * Math.PI * 2),
          [-1, 1],
          [-3, 3]
        );
        return `translate(${translateX}%, ${translateY}%)`;
      }

      default:
        return 'none';
    }
  };

  /**
   * Get animated gradient for gradientShift animation
   */
  const getAnimatedGradient = (): string | undefined => {
    if (animation !== 'gradientShift' || type !== 'gradient') return undefined;

    const cycleProgress = (frame % animationSpeed) / animationSpeed;
    const hueShift = cycleProgress * 360;

    // Shift the gradient angle
    const baseAngle = parseInt(getGradientAngle(gradientDirection)) || 180;
    const animatedAngle = (baseAngle + cycleProgress * 360) % 360;

    const colorStops = gradientColors.join(', ');
    return `linear-gradient(${animatedAngle}deg, ${colorStops})`;
  };

  /**
   * Get filter string
   */
  const getFilter = (): string => {
    const filters: string[] = [];

    if (blur > 0) {
      filters.push(`blur(${blur}px)`);
    }

    if (brightness !== 1) {
      filters.push(`brightness(${brightness})`);
    }

    return filters.length > 0 ? filters.join(' ') : 'none';
  };

  const backgroundStyle = getBackgroundStyle();
  const animatedGradient = getAnimatedGradient();
  const transform = getAnimationTransform();
  const filter = getFilter();

  return (
    <AbsoluteFill style={style}>
      {/* Background Layer */}
      <AbsoluteFill
        style={{
          ...backgroundStyle,
          ...(animatedGradient ? { background: animatedGradient } : {}),
          transform: animation !== 'none' ? transform : undefined,
          filter: filter !== 'none' ? filter : undefined,
          willChange: animation !== 'none' ? 'transform' : undefined,
        }}
      >
        {/* Image Background */}
        {type === 'image' && imageUrl && (
          <Img
            src={imageUrl}
            style={{
              width: '100%',
              height: '100%',
              objectFit: imageFit,
              objectPosition: imagePosition,
            }}
          />
        )}
      </AbsoluteFill>

      {/* Overlay Layer */}
      {overlay && (
        <AbsoluteFill
          style={{
            backgroundColor: overlay,
          }}
        />
      )}

      {/* Content Layer */}
      {children}
    </AbsoluteFill>
  );
};

/**
 * GradientBackground - Convenience component for gradient backgrounds
 */
export interface GradientBackgroundProps {
  /** Array of colors for the gradient */
  colors: string[];
  /** Direction of the gradient */
  direction?: GradientDirection;
  /** Enable gradient animation */
  animated?: boolean;
  /** Animation speed in frames */
  animationSpeed?: number;
  /** Children to render */
  children?: React.ReactNode;
}

export const GradientBackground: React.FC<GradientBackgroundProps> = ({
  colors,
  direction = 'to-bottom',
  animated = false,
  animationSpeed = 120,
  children,
}) => (
  <Background
    type="gradient"
    gradientColors={colors}
    gradientDirection={direction}
    animation={animated ? 'gradientShift' : 'none'}
    animationSpeed={animationSpeed}
  >
    {children}
  </Background>
);

/**
 * ImageBackground - Convenience component for image backgrounds
 */
export interface ImageBackgroundProps {
  /** Image URL */
  src: string;
  /** Image fit mode */
  fit?: 'cover' | 'contain' | 'fill' | 'none';
  /** Image position */
  position?: string;
  /** Overlay color */
  overlay?: string;
  /** Blur amount */
  blur?: number;
  /** Enable Ken Burns effect (slow zoom) */
  kenBurns?: boolean;
  /** Children to render */
  children?: React.ReactNode;
}

export const ImageBackground: React.FC<ImageBackgroundProps> = ({
  src,
  fit = 'cover',
  position = 'center',
  overlay,
  blur = 0,
  kenBurns = false,
  children,
}) => (
  <Background
    type="image"
    imageUrl={src}
    imageFit={fit}
    imagePosition={position}
    overlay={overlay}
    blur={blur}
    animation={kenBurns ? 'zoom' : 'none'}
    animationSpeed={300}
  >
    {children}
  </Background>
);

export default Background;
