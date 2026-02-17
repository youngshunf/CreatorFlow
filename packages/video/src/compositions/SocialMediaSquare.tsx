/**
 * SocialMediaSquare Composition
 *
 * 方形社交媒体视频组合（1:1），适用于 Instagram 动态和 Facebook。
 * 简洁设计，内容居中展示。
 *
 * @requirements 5.1, 5.5
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
 * Props for SocialMediaSquare composition
 */
export interface SocialMediaSquareProps extends TemplateProps {
  /** 主标题 */
  title?: string;
  /** 副标题 */
  subtitle?: string;
  /** Logo 图片路径 */
  logo?: string;
  /** 行动号召按钮 */
  cta?: {
    text: string;
    url?: string;
  };
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#ec4899',
  secondary: '#8b5cf6',
  background: '#0f0f23',
  text: '#ffffff',
};

/**
 * SocialMediaSquare - 方形社交媒体视频
 *
 * 优化的方形视频组合，适合 Instagram 动态和 Facebook 帖子
 */
export const SocialMediaSquare: React.FC<SocialMediaSquareProps> = ({
  title = '你的精彩标题',
  subtitle = '双击点赞',
  colors = DEFAULT_COLORS,
  logo,
  cta,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 标题动画
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const titleScale = interpolate(titleProgress, [0, 1], [0.9, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  // 副标题动画
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // 淡出
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 20, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  return (
    <AbsoluteFill
      style={{
        background: `radial-gradient(circle at center, ${colors.primary}30 0%, ${colors.background} 70%)`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
        padding: 60,
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
            left: 40,
            height: 50,
            opacity: titleProgress,
          }}
        />
      )}

      {/* 内容 */}
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          textAlign: 'center',
          gap: 24,
        }}
      >
        <h1
          style={{
            fontSize: 56,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleProgress,
            transform: `scale(${titleScale})`,
            lineHeight: 1.2,
          }}
        >
          {title}
        </h1>

        <p
          style={{
            fontSize: 24,
            color: colors.primary,
            margin: 0,
            opacity: Math.max(0, subtitleProgress),
          }}
        >
          {subtitle}
        </p>

        {cta && (
          <div
            style={{
              marginTop: 30,
              padding: '16px 40px',
              backgroundColor: colors.primary,
              borderRadius: 30,
              fontSize: 20,
              fontWeight: 'bold',
              color: colors.text,
              opacity: Math.max(0, subtitleProgress),
            }}
          >
            {cta.text}
          </div>
        )}
      </div>
    </AbsoluteFill>
  );
};

/**
 * 默认 props
 */
export const defaultProps: SocialMediaSquareProps = {
  title: '你的精彩标题',
  subtitle: '双击点赞',
  colors: DEFAULT_COLORS,
  cta: {
    text: '立即购买',
  },
};

export default SocialMediaSquare;
