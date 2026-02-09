/**
 * Transition Component
 *
 * A reusable component for scene transition effects.
 * Supports various transition styles like fade, wipe, slide, zoom, and dissolve.
 *
 * @requirements 4.2, 4.5
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

/**
 * Transition type options
 */
export type TransitionType =
  | 'fade'
  | 'wipeLeft'
  | 'wipeRight'
  | 'wipeUp'
  | 'wipeDown'
  | 'slideLeft'
  | 'slideRight'
  | 'slideUp'
  | 'slideDown'
  | 'zoom'
  | 'dissolve'
  | 'iris'
  | 'clock';

/**
 * Props for Transition component
 */
export interface TransitionProps {
  /** Type of transition effect */
  type?: TransitionType;
  /** Duration of the transition in frames */
  duration?: number;
  /** Frame at which the transition starts */
  startFrame?: number;
  /** Direction of the transition: 'in' (appearing) or 'out' (disappearing) */
  direction?: 'in' | 'out';
  /** Color for fade/wipe transitions */
  color?: string;
  /** Children to wrap with transition effect */
  children?: React.ReactNode;
  /** Spring configuration for spring-based transitions */
  springConfig?: {
    damping?: number;
    stiffness?: number;
    mass?: number;
  };
  /** Easing function for interpolation */
  easing?: (t: number) => number;
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
 * Transition - Provides scene transition effects
 *
 * This component wraps content and applies transition effects,
 * useful for scene changes and visual transitions in video compositions.
 *
 * @example
 * ```tsx
 * <Transition type="fade" duration={30} direction="in">
 *   <YourContent />
 * </Transition>
 * ```
 */
export const Transition: React.FC<TransitionProps> = ({
  type = 'fade',
  duration = 30,
  startFrame = 0,
  direction = 'in',
  color = '#000000',
  children,
  springConfig = DEFAULT_SPRING_CONFIG,
  easing,
}) => {
  const frame = useCurrentFrame();
  const { fps, width, height } = useVideoConfig();

  // Calculate progress (0 to 1)
  const rawProgress = interpolate(
    frame,
    [startFrame, startFrame + duration],
    [0, 1],
    {
      extrapolateLeft: 'clamp',
      extrapolateRight: 'clamp',
    }
  );

  // Apply easing if provided
  const progress = easing ? easing(rawProgress) : rawProgress;

  // Reverse progress for 'out' direction
  const effectiveProgress = direction === 'out' ? 1 - progress : progress;

  /**
   * Get opacity for fade-based transitions
   */
  const getOpacity = (): number => {
    switch (type) {
      case 'fade':
      case 'dissolve':
        return effectiveProgress;
      default:
        return 1;
    }
  };

  /**
   * Get transform for movement-based transitions
   */
  const getTransform = (): string => {
    switch (type) {
      case 'slideLeft':
        return `translateX(${interpolate(effectiveProgress, [0, 1], [100, 0])}%)`;
      case 'slideRight':
        return `translateX(${interpolate(effectiveProgress, [0, 1], [-100, 0])}%)`;
      case 'slideUp':
        return `translateY(${interpolate(effectiveProgress, [0, 1], [100, 0])}%)`;
      case 'slideDown':
        return `translateY(${interpolate(effectiveProgress, [0, 1], [-100, 0])}%)`;
      case 'zoom':
        const scale = interpolate(effectiveProgress, [0, 1], [0.5, 1]);
        return `scale(${scale})`;
      default:
        return 'none';
    }
  };

  /**
   * Get clip path for wipe transitions
   */
  const getClipPath = (): string | undefined => {
    switch (type) {
      case 'wipeLeft':
        return `inset(0 ${interpolate(effectiveProgress, [0, 1], [100, 0])}% 0 0)`;
      case 'wipeRight':
        return `inset(0 0 0 ${interpolate(effectiveProgress, [0, 1], [100, 0])}%)`;
      case 'wipeUp':
        return `inset(${interpolate(effectiveProgress, [0, 1], [100, 0])}% 0 0 0)`;
      case 'wipeDown':
        return `inset(0 0 ${interpolate(effectiveProgress, [0, 1], [100, 0])}% 0)`;
      case 'iris': {
        const radius = interpolate(effectiveProgress, [0, 1], [0, 150]);
        return `circle(${radius}% at 50% 50%)`;
      }
      case 'clock': {
        // Clock wipe using conic gradient simulation
        const angle = interpolate(effectiveProgress, [0, 1], [0, 360]);
        // Use polygon approximation for clock wipe
        const points = generateClockWipePoints(angle);
        return `polygon(${points})`;
      }
      default:
        return undefined;
    }
  };

  /**
   * Get filter for dissolve effect
   */
  const getFilter = (): string | undefined => {
    if (type === 'dissolve') {
      const blur = interpolate(effectiveProgress, [0, 0.5, 1], [10, 5, 0]);
      return `blur(${blur}px)`;
    }
    return undefined;
  };

  const opacity = getOpacity();
  const transform = getTransform();
  const clipPath = getClipPath();
  const filter = getFilter();

  return (
    <AbsoluteFill
      style={{
        opacity,
        transform,
        clipPath,
        filter,
        willChange: 'opacity, transform, clip-path, filter',
      }}
    >
      {children}
    </AbsoluteFill>
  );
};

/**
 * Generate polygon points for clock wipe effect
 */
function generateClockWipePoints(angle: number): string {
  if (angle <= 0) return '50% 50%, 50% 50%, 50% 50%';
  if (angle >= 360) return '50% 50%, 50% 0%, 100% 0%, 100% 100%, 0% 100%, 0% 0%, 50% 0%';

  const points: string[] = ['50% 50%', '50% 0%'];
  const corners = [
    { angle: 45, point: '100% 0%' },
    { angle: 135, point: '100% 100%' },
    { angle: 225, point: '0% 100%' },
    { angle: 315, point: '0% 0%' },
  ];

  for (const corner of corners) {
    if (angle >= corner.angle) {
      points.push(corner.point);
    }
  }

  // Calculate the final point on the edge
  const radians = ((angle - 90) * Math.PI) / 180;
  const edgePoint = getEdgePoint(radians);
  points.push(edgePoint);

  return points.join(', ');
}

/**
 * Get the point on the edge of the container for a given angle
 */
function getEdgePoint(radians: number): string {
  const tan = Math.tan(radians);
  const cos = Math.cos(radians);
  const sin = Math.sin(radians);

  let x: number, y: number;

  if (Math.abs(cos) > Math.abs(sin)) {
    // Intersects with left or right edge
    x = cos > 0 ? 100 : 0;
    y = 50 + (x - 50) * tan;
  } else {
    // Intersects with top or bottom edge
    y = sin > 0 ? 100 : 0;
    x = 50 + (y - 50) / tan;
  }

  // Clamp values
  x = Math.max(0, Math.min(100, x));
  y = Math.max(0, Math.min(100, y));

  return `${x}% ${y}%`;
}

/**
 * TransitionOverlay - A standalone transition overlay
 *
 * Use this for creating transition overlays that cover the entire screen,
 * useful for scene-to-scene transitions.
 */
export interface TransitionOverlayProps {
  /** Type of transition effect */
  type?: 'fade' | 'wipe' | 'iris';
  /** Duration of the transition in frames */
  duration?: number;
  /** Frame at which the transition starts */
  startFrame?: number;
  /** Direction: 'in' fades from color, 'out' fades to color */
  direction?: 'in' | 'out';
  /** Overlay color */
  color?: string;
}

export const TransitionOverlay: React.FC<TransitionOverlayProps> = ({
  type = 'fade',
  duration = 30,
  startFrame = 0,
  direction = 'in',
  color = '#000000',
}) => {
  const frame = useCurrentFrame();

  const progress = interpolate(
    frame,
    [startFrame, startFrame + duration],
    [0, 1],
    {
      extrapolateLeft: 'clamp',
      extrapolateRight: 'clamp',
    }
  );

  // For 'in' direction: overlay fades out (content appears)
  // For 'out' direction: overlay fades in (content disappears)
  const opacity = direction === 'in' ? 1 - progress : progress;

  if (opacity <= 0) return null;

  const getClipPath = (): string | undefined => {
    if (type === 'iris') {
      const radius = direction === 'in'
        ? interpolate(progress, [0, 1], [150, 0])
        : interpolate(progress, [0, 1], [0, 150]);
      return `circle(${radius}% at 50% 50%)`;
    }
    if (type === 'wipe') {
      const wipeProgress = direction === 'in' ? progress : 1 - progress;
      return `inset(0 ${wipeProgress * 100}% 0 0)`;
    }
    return undefined;
  };

  return (
    <AbsoluteFill
      style={{
        backgroundColor: color,
        opacity: type === 'fade' ? opacity : 1,
        clipPath: getClipPath(),
        zIndex: 1000,
      }}
    />
  );
};

export default Transition;
