/**
 * PromoAd Composition
 *
 * 促销广告视频组合，醒目的促销视频，适用于打折、特卖和限时优惠。
 * 大胆设计，吸引眼球。
 *
 * @requirements 5.2, 5.5
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
 * Props for PromoAd composition
 */
export interface PromoAdProps extends TemplateProps {
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
  /** 折扣标签 */
  discount?: string;
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#f43f5e',
  secondary: '#fb923c',
  background: '#1a1a2e',
  text: '#ffffff',
};

/**
 * PromoAd - 促销广告视频
 *
 * 大胆、吸引眼球的设计，适合广告
 */
export const PromoAd: React.FC<PromoAdProps> = ({
  title = '超级大促',
  subtitle = '限时优惠，不容错过！',
  colors = DEFAULT_COLORS,
  logo,
  cta = { text: '立即抢购' },
  discount = '50% OFF',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 脉冲动画（用于折扣标签）
  const pulse = Math.sin(frame * 0.15) * 0.05 + 1;

  // 标题动画
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 100, stiffness: 200, mass: 0.5 },
  });

  const titleScale = interpolate(titleProgress, [0, 1], [0.8, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  // 折扣动画
  const discountProgress = spring({
    frame: frame - 10,
    fps,
    config: { damping: 100, stiffness: 200, mass: 0.5 },
  });

  const discountScale = interpolate(discountProgress, [0, 1], [0.5, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  // CTA 动画
  const ctaProgress = spring({
    frame: frame - 25,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const ctaScale = interpolate(ctaProgress, [0, 1], [0.9, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
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
        background: `linear-gradient(135deg, ${colors.background} 0%, ${colors.primary}30 100%)`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
      }}
    >
      {/* 动画背景形状 */}
      <div
        style={{
          position: 'absolute',
          top: -100,
          right: -100,
          width: 400,
          height: 400,
          borderRadius: '50%',
          background: `radial-gradient(circle, ${colors.primary}40 0%, transparent 70%)`,
          transform: `scale(${pulse})`,
        }}
      />
      <div
        style={{
          position: 'absolute',
          bottom: -150,
          left: -150,
          width: 500,
          height: 500,
          borderRadius: '50%',
          background: `radial-gradient(circle, ${colors.secondary}30 0%, transparent 70%)`,
          transform: `scale(${pulse * 1.1})`,
        }}
      />

      {/* Logo */}
      {logo && (
        <Img
          src={staticFile(logo)}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 40,
            left: 60,
            height: 50,
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
          gap: 20,
        }}
      >
        {/* 折扣徽章 */}
        <div
          style={{
            padding: '12px 30px',
            backgroundColor: colors.primary,
            borderRadius: 50,
            fontSize: 28,
            fontWeight: 'bold',
            color: colors.text,
            opacity: discountProgress,
            transform: `scale(${discountScale * pulse})`,
            boxShadow: `0 10px 30px ${colors.primary}60`,
          }}
        >
          {discount}
        </div>

        {/* 标题 */}
        <h1
          style={{
            fontSize: 96,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleProgress,
            transform: `scale(${titleScale})`,
            textShadow: `0 4px 20px ${colors.primary}80`,
            letterSpacing: '0.05em',
          }}
        >
          {title}
        </h1>

        {/* 副标题 */}
        <p
          style={{
            fontSize: 32,
            color: 'rgba(255, 255, 255, 0.8)',
            margin: 0,
            opacity: titleProgress,
          }}
        >
          {subtitle}
        </p>

        {/* CTA */}
        <div
          style={{
            marginTop: 30,
            padding: '20px 60px',
            background: `linear-gradient(135deg, ${colors.primary} 0%, ${colors.secondary} 100%)`,
            borderRadius: 12,
            fontSize: 28,
            fontWeight: 'bold',
            color: colors.text,
            opacity: ctaProgress,
            transform: `scale(${ctaScale})`,
            boxShadow: `0 15px 40px ${colors.primary}50`,
          }}
        >
          {cta.text}
        </div>
      </div>
    </AbsoluteFill>
  );
};

/**
 * 默认 props
 */
export const defaultProps: PromoAdProps = {
  title: '超级大促',
  subtitle: '限时优惠，不容错过！',
  colors: DEFAULT_COLORS,
  cta: {
    text: '立即抢购',
  },
  discount: '50% OFF',
};

export default PromoAd;
