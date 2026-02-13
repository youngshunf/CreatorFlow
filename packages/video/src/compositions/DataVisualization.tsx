/**
 * DataVisualization Composition
 *
 * A video composition that displays animated data charts and graphs.
 * Supports bar charts, line charts, and pie charts with smooth animations.
 *
 * @requirements 3.3, 3.5, 3.6
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Img,
  staticFile,
} from 'remotion';
import type { TemplateProps } from '../templates/types';

/**
 * Data point interface
 */
export interface DataPoint {
  label: string;
  value: number;
  color?: string;
}

/**
 * Props for DataVisualization composition
 */
export type DataVisualizationProps = TemplateProps & {
  /** Chart title */
  title?: string;
  /** Data points to visualize */
  data?: DataPoint[];
  /** Chart type */
  chartType?: 'bar' | 'line' | 'pie' | 'donut';
  /** Show data labels */
  showLabels?: boolean;
  /** Show values on chart */
  showValues?: boolean;
  /** Animation delay between data points (in frames) */
  staggerDelay?: number;
};

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
 * Default chart colors
 */
const CHART_COLORS = [
  '#6366f1',
  '#8b5cf6',
  '#ec4899',
  '#f43f5e',
  '#f97316',
  '#eab308',
  '#22c55e',
  '#14b8a6',
];

/**
 * Default data for demo purposes
 */
const DEFAULT_DATA: DataPoint[] = [
  { label: 'Jan', value: 65 },
  { label: 'Feb', value: 85 },
  { label: 'Mar', value: 45 },
  { label: 'Apr', value: 95 },
  { label: 'May', value: 75 },
  { label: 'Jun', value: 55 },
];

/**
 * DataVisualization - Displays animated data charts
 *
 * This composition creates professional data visualizations with
 * smooth animations for presenting statistics and metrics.
 */
export const DataVisualization: React.FC<DataVisualizationProps> = ({
  title = 'Data Overview',
  subtitle,
  data = DEFAULT_DATA,
  colors = DEFAULT_COLORS,
  chartType = 'bar',
  showLabels = true,
  showValues = true,
  staggerDelay = 5,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Calculate max value for scaling
  const maxValue = Math.max(...data.map((d) => d.value));

  // Title animation
  const titleOpacity = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const titleTranslateY = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
    from: -30,
    to: 0,
  });

  // Fade out at the end
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 30, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  // Render bar chart
  const renderBarChart = () => {
    const barWidth = 60;
    const barGap = 30;
    const chartHeight = 400;
    const chartWidth = data.length * (barWidth + barGap);

    return (
      <div
        style={{
          display: 'flex',
          alignItems: 'flex-end',
          justifyContent: 'center',
          height: chartHeight,
          gap: barGap,
          marginTop: 60,
        }}
      >
        {data.map((point, index) => {
          const delay = index * staggerDelay;
          const barProgress = spring({
            frame: frame - delay - 30,
            fps,
            config: { damping: 200, stiffness: 80, mass: 0.8 },
          });

          const barHeight = (point.value / maxValue) * chartHeight * barProgress;
          const barColor = point.color || CHART_COLORS[index % CHART_COLORS.length];

          return (
            <div
              key={index}
              style={{
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                gap: 10,
              }}
            >
              {/* Value label */}
              {showValues && (
                <span
                  style={{
                    fontSize: 20,
                    fontWeight: 'bold',
                    color: colors.text,
                    opacity: barProgress,
                  }}
                >
                  {Math.round(point.value * barProgress)}
                </span>
              )}

              {/* Bar */}
              <div
                style={{
                  width: barWidth,
                  height: barHeight,
                  backgroundColor: barColor,
                  borderRadius: '8px 8px 0 0',
                  boxShadow: '0 4px 12px rgba(0, 0, 0, 0.3)',
                }}
              />

              {/* Label */}
              {showLabels && (
                <span
                  style={{
                    fontSize: 16,
                    color: 'rgba(255, 255, 255, 0.7)',
                    marginTop: 10,
                  }}
                >
                  {point.label}
                </span>
              )}
            </div>
          );
        })}
      </div>
    );
  };

  // Render line chart
  const renderLineChart = () => {
    const chartWidth = 800;
    const chartHeight = 400;
    const padding = 40;

    const points = data.map((point, index) => {
      const x = padding + (index / (data.length - 1)) * (chartWidth - padding * 2);
      const y = chartHeight - padding - (point.value / maxValue) * (chartHeight - padding * 2);
      return { x, y, ...point };
    });

    // Animate line drawing
    const lineProgress = spring({
      frame: frame - 30,
      fps,
      config: { damping: 200, stiffness: 60, mass: 1 },
    });

    // Create SVG path
    const pathData = points
      .map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x} ${p.y}`)
      .join(' ');

    return (
      <svg
        width={chartWidth}
        height={chartHeight}
        style={{ marginTop: 40 }}
      >
        {/* Grid lines */}
        {[0, 0.25, 0.5, 0.75, 1].map((ratio, i) => (
          <line
            key={i}
            x1={padding}
            y1={padding + ratio * (chartHeight - padding * 2)}
            x2={chartWidth - padding}
            y2={padding + ratio * (chartHeight - padding * 2)}
            stroke="rgba(255, 255, 255, 0.1)"
            strokeWidth={1}
          />
        ))}

        {/* Line */}
        <path
          d={pathData}
          fill="none"
          stroke={colors.primary}
          strokeWidth={4}
          strokeLinecap="round"
          strokeLinejoin="round"
          strokeDasharray={1000}
          strokeDashoffset={1000 * (1 - lineProgress)}
        />

        {/* Data points */}
        {points.map((point, index) => {
          const delay = index * staggerDelay;
          const pointProgress = spring({
            frame: frame - delay - 45,
            fps,
            config: { damping: 200, stiffness: 100, mass: 0.5 },
          });

          return (
            <g key={index}>
              <circle
                cx={point.x}
                cy={point.y}
                r={8 * pointProgress}
                fill={colors.primary}
                stroke={colors.background}
                strokeWidth={3}
              />
              {showLabels && (
                <text
                  x={point.x}
                  y={chartHeight - 10}
                  textAnchor="middle"
                  fill="rgba(255, 255, 255, 0.7)"
                  fontSize={14}
                  opacity={pointProgress}
                >
                  {point.label}
                </text>
              )}
              {showValues && (
                <text
                  x={point.x}
                  y={point.y - 20}
                  textAnchor="middle"
                  fill={colors.text}
                  fontSize={16}
                  fontWeight="bold"
                  opacity={pointProgress}
                >
                  {point.value}
                </text>
              )}
            </g>
          );
        })}
      </svg>
    );
  };

  // Render pie/donut chart
  const renderPieChart = () => {
    const size = 400;
    const centerX = size / 2;
    const centerY = size / 2;
    const radius = chartType === 'donut' ? 140 : 160;
    const innerRadius = chartType === 'donut' ? 80 : 0;

    const total = data.reduce((sum, d) => sum + d.value, 0);
    let currentAngle = -90; // Start from top

    return (
      <div style={{ display: 'flex', alignItems: 'center', gap: 60, marginTop: 40 }}>
        <svg width={size} height={size}>
          {data.map((point, index) => {
            const delay = index * staggerDelay;
            const sliceProgress = spring({
              frame: frame - delay - 30,
              fps,
              config: { damping: 200, stiffness: 80, mass: 0.8 },
            });

            const angle = (point.value / total) * 360 * sliceProgress;
            const startAngle = currentAngle;
            const endAngle = currentAngle + angle;
            currentAngle = endAngle;

            const startRad = (startAngle * Math.PI) / 180;
            const endRad = (endAngle * Math.PI) / 180;

            const x1 = centerX + radius * Math.cos(startRad);
            const y1 = centerY + radius * Math.sin(startRad);
            const x2 = centerX + radius * Math.cos(endRad);
            const y2 = centerY + radius * Math.sin(endRad);

            const largeArc = angle > 180 ? 1 : 0;
            const sliceColor = point.color || CHART_COLORS[index % CHART_COLORS.length];

            let pathD: string;
            if (chartType === 'donut') {
              const ix1 = centerX + innerRadius * Math.cos(startRad);
              const iy1 = centerY + innerRadius * Math.sin(startRad);
              const ix2 = centerX + innerRadius * Math.cos(endRad);
              const iy2 = centerY + innerRadius * Math.sin(endRad);

              pathD = `
                M ${x1} ${y1}
                A ${radius} ${radius} 0 ${largeArc} 1 ${x2} ${y2}
                L ${ix2} ${iy2}
                A ${innerRadius} ${innerRadius} 0 ${largeArc} 0 ${ix1} ${iy1}
                Z
              `;
            } else {
              pathD = `
                M ${centerX} ${centerY}
                L ${x1} ${y1}
                A ${radius} ${radius} 0 ${largeArc} 1 ${x2} ${y2}
                Z
              `;
            }

            return (
              <path
                key={index}
                d={pathD}
                fill={sliceColor}
                stroke={colors.background}
                strokeWidth={2}
              />
            );
          })}

          {/* Center text for donut */}
          {chartType === 'donut' && (
            <text
              x={centerX}
              y={centerY}
              textAnchor="middle"
              dominantBaseline="middle"
              fill={colors.text}
              fontSize={36}
              fontWeight="bold"
            >
              {total}
            </text>
          )}
        </svg>

        {/* Legend */}
        {showLabels && (
          <div style={{ display: 'flex', flexDirection: 'column', gap: 15 }}>
            {data.map((point, index) => {
              const delay = index * staggerDelay;
              const legendProgress = spring({
                frame: frame - delay - 30,
                fps,
                config: { damping: 200, stiffness: 100, mass: 0.5 },
              });

              const sliceColor = point.color || CHART_COLORS[index % CHART_COLORS.length];

              return (
                <div
                  key={index}
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 12,
                    opacity: legendProgress,
                    transform: `translateX(${20 * (1 - legendProgress)}px)`,
                  }}
                >
                  <div
                    style={{
                      width: 16,
                      height: 16,
                      borderRadius: 4,
                      backgroundColor: sliceColor,
                    }}
                  />
                  <span style={{ color: colors.text, fontSize: 18 }}>
                    {point.label}
                  </span>
                  {showValues && (
                    <span style={{ color: 'rgba(255, 255, 255, 0.6)', fontSize: 16 }}>
                      ({point.value})
                    </span>
                  )}
                </div>
              );
            })}
          </div>
        )}
      </div>
    );
  };

  // Render chart based on type
  const renderChart = () => {
    switch (chartType) {
      case 'bar':
        return renderBarChart();
      case 'line':
        return renderLineChart();
      case 'pie':
      case 'donut':
        return renderPieChart();
      default:
        return renderBarChart();
    }
  };

  return (
    <AbsoluteFill
      style={{
        backgroundColor: colors.background,
        fontFamily: 'system-ui, -apple-system, sans-serif',
        padding: 60,
        opacity: fadeOut,
      }}
    >
      {/* Header */}
      <div
        style={{
          opacity: titleOpacity,
          transform: `translateY(${titleTranslateY}px)`,
        }}
      >
        <h1
          style={{
            fontSize: 56,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            marginBottom: 10,
          }}
        >
          {title}
        </h1>
        {subtitle && (
          <p
            style={{
              fontSize: 24,
              color: colors.primary,
              margin: 0,
            }}
          >
            {subtitle}
          </p>
        )}
      </div>

      {/* Chart */}
      <div
        style={{
          flex: 1,
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center',
        }}
      >
        {renderChart()}
      </div>

      {/* Optional Logo */}
      {logo && (
        <Img
          src={staticFile(logo)}
          alt="Logo"
          style={{
            position: 'absolute',
            bottom: 40,
            right: 40,
            height: 40,
            opacity: titleOpacity,
          }}
        />
      )}
    </AbsoluteFill>
  );
};

export default DataVisualization;
