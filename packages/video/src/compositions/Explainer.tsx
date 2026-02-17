/**
 * Explainer Composition
 *
 * 概念讲解视频组合，简洁的讲解模板，适用于教育内容。
 * 将复杂主题拆解为易懂的要点。
 *
 * @requirements 5.3, 5.5
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
 * 讲解要点接口
 */
export interface ExplainerItem {
  title: string;
  description?: string;
  icon?: string;
}

/**
 * Props for Explainer composition
 */
export interface ExplainerProps extends TemplateProps {
  /** 主标题 */
  title?: string;
  /** 副标题 */
  subtitle?: string;
  /** 要点列表 */
  items?: ExplainerItem[];
  /** Logo 图片路径 */
  logo?: string;
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#6366f1',
  secondary: '#8b5cf6',
  background: '#ffffff',
  text: '#1f2937',
};

/**
 * Explainer - 概念讲解视频
 *
 * 清晰的设计用于解释概念
 */
export const Explainer: React.FC<ExplainerProps> = ({
  title = '概念解析',
  subtitle = '简单易懂的讲解',
  items = [
    { title: '核心要点一', description: '你需要理解的基础概念', icon: '💡' },
    { title: '核心要点二', description: '在基础上深入了解更多细节', icon: '📊' },
    { title: '核心要点三', description: '实际应用与关键收获', icon: '🎯' },
  ],
  colors = DEFAULT_COLORS,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 标题动画
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // 淡出
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 30, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  return (
    <AbsoluteFill
      style={{
        backgroundColor: colors.background,
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
        padding: 80,
      }}
    >
      {/* Logo */}
      {logo && (
        <Img
          src={staticFile(logo)}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 40,
            left: 60,
            height: 40,
            opacity: titleProgress,
          }}
        />
      )}

      {/* 头部 */}
      <div
        style={{
          marginBottom: 60,
          paddingTop: 40,
        }}
      >
        <h1
          style={{
            fontSize: 56,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            marginBottom: 16,
            opacity: titleProgress,
          }}
        >
          {title}
        </h1>
        <p
          style={{
            fontSize: 28,
            color: colors.primary,
            margin: 0,
            opacity: titleProgress,
          }}
        >
          {subtitle}
        </p>
      </div>

      {/* 内容网格 */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(3, 1fr)',
          gap: 40,
          flex: 1,
        }}
      >
        {items.map((item, index) => {
          const itemProgress = spring({
            frame: frame - 30 - index * 15,
            fps,
            config: { damping: 200, stiffness: 100, mass: 0.5 },
          });

          const itemTranslateY = interpolate(itemProgress, [0, 1], [30, 0], {
            extrapolateLeft: 'clamp',
            extrapolateRight: 'clamp',
          });

          return (
            <div
              key={index}
              style={{
                padding: 40,
                backgroundColor: `${colors.primary}08`,
                borderRadius: 20,
                border: `2px solid ${colors.primary}20`,
                opacity: Math.max(0, itemProgress),
                transform: `translateY(${itemTranslateY}px)`,
              }}
            >
              <span
                style={{
                  fontSize: 48,
                  display: 'block',
                  marginBottom: 20,
                }}
              >
                {item.icon || '📌'}
              </span>
              <h3
                style={{
                  fontSize: 28,
                  fontWeight: 'bold',
                  color: colors.text,
                  margin: 0,
                  marginBottom: 12,
                }}
              >
                {item.title}
              </h3>
              {item.description && (
                <p
                  style={{
                    fontSize: 18,
                    color: `${colors.text}99`,
                    margin: 0,
                    lineHeight: 1.5,
                  }}
                >
                  {item.description}
                </p>
              )}
            </div>
          );
        })}
      </div>
    </AbsoluteFill>
  );
};

/**
 * 默认 props
 */
export const defaultProps: ExplainerProps = {
  title: '概念解析',
  subtitle: '简单易懂的讲解',
  colors: DEFAULT_COLORS,
  items: [
    {
      title: '核心要点一',
      description: '你需要理解的基础概念',
      icon: '💡',
    },
    {
      title: '核心要点二',
      description: '在基础上深入了解更多细节',
      icon: '📊',
    },
    {
      title: '核心要点三',
      description: '实际应用与关键收获',
      icon: '🎯',
    },
  ],
};

export default Explainer;
