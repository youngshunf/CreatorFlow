/**
 * BarItem Component
 *
 * Displays an animated bar chart item with spring animation.
 * Designed to be used within a Sequence for staggered animations.
 */
import React from 'react';
import { useCurrentFrame, useVideoConfig, spring, interpolate } from 'remotion';
import type { DataPoint } from '../compositions/DataVisualization';

export type BarItemProps = {
  /** Data point to display */
  point: DataPoint;
  /** Maximum value for scaling */
  maxValue: number;
  /** Chart height in pixels */
  chartHeight: number;
  /** Show label */
  showLabel?: boolean;
  /** Show value */
  showValue?: boolean;
  /** Primary color */
  primaryColor?: string;
};

/**
 * BarItem - Animated bar chart item
 *
 * This component handles the animation of a single bar in a bar chart.
 * It uses spring animation for smooth, natural motion.
 */
export const BarItem: React.FC<BarItemProps> = ({
  point,
  maxValue,
  chartHeight,
  showLabel = true,
  showValue = true,
  primaryColor = '#6366f1',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Bar animation with spring
  const barProgress = spring({
    frame,
    fps,
    config: {
      damping: 200,
      stiffness: 80,
      mass: 0.8,
    },
  });

  // Calculate bar height
  const barHeight = (point.value / maxValue) * chartHeight * barProgress;

  // Label fade-in
  const labelOpacity = interpolate(frame, [10, 30], [0, 1], {
    extrapolateRight: 'clamp',
  });

  // Value fade-in (delayed)
  const valueOpacity = interpolate(frame, [20, 40], [0, 1], {
    extrapolateRight: 'clamp',
  });

  return (
    <div
      style={{
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        gap: 10,
        flex: 1,
      }}
    >
      {/* Value */}
      {showValue && (
        <div
          style={{
            fontSize: 18,
            fontWeight: 'bold',
            color: point.color || primaryColor,
            opacity: valueOpacity,
            minHeight: 25,
          }}
        >
          {Math.round(point.value * barProgress)}
        </div>
      )}

      {/* Bar */}
      <div
        style={{
          width: '100%',
          height: chartHeight,
          display: 'flex',
          alignItems: 'flex-end',
          justifyContent: 'center',
        }}
      >
        <div
          style={{
            width: '80%',
            height: barHeight,
            backgroundColor: point.color || primaryColor,
            borderRadius: '8px 8px 0 0',
            transition: 'height 0.3s ease',
          }}
        />
      </div>

      {/* Label */}
      {showLabel && (
        <div
          style={{
            fontSize: 16,
            color: '#ffffff',
            opacity: labelOpacity,
            textAlign: 'center',
          }}
        >
          {point.label}
        </div>
      )}
    </div>
  );
};
