/**
 * SocialMediaVertical Composition
 *
 * 竖版社交媒体视频组合（9:16），适用于抖音、Instagram Reels 和 YouTube Shorts。
 * 包含吸睛动画和移动端优先设计。
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
 * Props for SocialMediaVertical composition
 */
export interface SocialMediaVerticalProps extends TemplateProps {
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
 * SocialMediaVertical - 竖版社交媒体视频
 *
 * 优化的竖版视频组合，使用 spring 动画，适合移动端观看
 */
export const SocialMediaVertical: React.FC<SocialMediaVerticalProps> = ({
  title = '你的精彩标题',
  subtitle = '上滑了解更多',
  colors = DEFAULT_COLORS,
  logo,
  cta,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 标题动画（spring）
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const titleTranslateY = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
    from: 50,
    to: 0,
  });

  // 副标题动画（延迟）
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // CTA 动画（延迟更多）
  const ctaProgress = spring({
    frame: frame - 30,
    fps,
    config: { damping: 200, stiffness: 80, mass: 1 },
  });

  const ctaScale = interpolate(ctaProgress, [0, 1], [0.8, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  // 结尾淡出
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 20, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  // 渐变背景动画
  const gradientRotation = interpolate(frame, [0, durationInFrames], [0, 360]);

  return (
    <AbsoluteFill
      style={{
        background: `linear-gradient(${gradientRotation}deg, ${colors.background} 0%, ${colors.primary}20 50%, ${colors.background} 100%)`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
        padding: 40,
      }}
    >
      {/* Logo 在顶部 */}
      {logo && (
        <Img
          src={staticFile(logo)}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 60,
            height: 60,
            opacity: titleProgress,
          }}
        />
      )}

      {/* 主内容 */}
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          textAlign: 'center',
          gap: 30,
        }}
      >
        {/* 标题 */}
        <h1
          style={{
            fontSize: 64,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleProgress,
            transform: `translateY(${titleTranslateY}px)`,
            lineHeight: 1.2,
            maxWidth: '90%',
          }}
        >
          {title}
        </h1>

        {/* 副标题 */}
        <p
          style={{
            fontSize: 28,
            color: colors.primary,
            margin: 0,
            opacity: Math.max(0, subtitleProgress),
            maxWidth: '80%',
          }}
        >
          {subtitle}
        </p>

        {/* CTA 按钮 */}
        {cta && (
          <div
            style={{
              marginTop: 40,
              padding: '20px 50px',
              backgroundColor: colors.primary,
              borderRadius: 50,
              fontSize: 24,
              fontWeight: 'bold',
              color: colors.text,
              opacity: Math.max(0, ctaProgress),
              transform: `scale(${ctaScale})`,
              boxShadow: `0 10px 30px ${colors.primary}60`,
            }}
          >
            {cta.text}
          </div>
        )}
      </div>

      {/* 装饰元素 */}
      <div
        style={{
          position: 'absolute',
          bottom: 100,
          width: 60,
          height: 4,
          backgroundColor: colors.primary,
          borderRadius: 2,
          opacity: ctaProgress,
        }}
      />
    </AbsoluteFill>
  );
};

/**
 * 默认 props
 */
export const defaultProps: SocialMediaVerticalProps = {
  title: '你的精彩标题',
  subtitle: '上滑了解更多',
  colors: DEFAULT_COLORS,
  cta: {
    text: '了解更多',
  },
};

export default SocialMediaVertical;
