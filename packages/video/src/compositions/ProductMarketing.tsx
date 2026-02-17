/**
 * ProductMarketing Composition
 *
 * 产品营销视频组合，专业产品营销视频，突出功能亮点。
 * 适用于产品发布、演示和推广内容。
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
 * 产品功能接口
 */
export interface ProductFeature {
  title: string;
  description?: string;
  icon?: string;
}

/**
 * Props for ProductMarketing composition
 */
export interface ProductMarketingProps extends TemplateProps {
  /** 主标题 */
  title?: string;
  /** 副标题 */
  subtitle?: string;
  /** 功能列表 */
  items?: ProductFeature[];
  /** Logo 图片路径 */
  logo?: string;
  /** 行动号召按钮 */
  cta?: {
    text: string;
    url?: string;
  };
  /** 产品图片路径 */
  productImage?: string;
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#0ea5e9',
  secondary: '#0284c7',
  background: '#0f172a',
  text: '#f8fafc',
};

/**
 * ProductMarketing - 产品营销视频
 *
 * 专业产品展示，带有动画功能列表
 */
export const ProductMarketing: React.FC<ProductMarketingProps> = ({
  title = '全新产品发布',
  subtitle = '创新引领未来',
  items = [
    { title: '功能一', description: '惊人的能力', icon: '🚀' },
    { title: '功能二', description: '强大的性能', icon: '⚡' },
    { title: '功能三', description: '智能设计', icon: '🎯' },
  ],
  colors = DEFAULT_COLORS,
  logo,
  cta,
  productImage,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 标题动画
  const titleProgress = spring({
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

  // 副标题动画
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // 产品图片动画
  const productProgress = spring({
    frame: frame - 30,
    fps,
    config: { damping: 200, stiffness: 80, mass: 1 },
  });

  const productScale = interpolate(productProgress, [0, 1], [0.9, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
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
      }}
    >
      {/* 背景渐变 */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          background: `radial-gradient(ellipse at top right, ${colors.primary}20 0%, transparent 50%),
                       radial-gradient(ellipse at bottom left, ${colors.secondary}20 0%, transparent 50%)`,
        }}
      />

      {/* 内容容器 */}
      <div
        style={{
          display: 'flex',
          height: '100%',
          padding: 80,
        }}
      >
        {/* 左侧 - 文本内容 */}
        <div
          style={{
            flex: 1,
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center',
            paddingRight: 60,
          }}
        >
          {/* Logo */}
          {logo && (
            <Img
              src={staticFile(logo)}
              alt="Logo"
              style={{
                height: 50,
                marginBottom: 40,
                opacity: titleProgress,
              }}
            />
          )}

          {/* 标题 */}
          <h1
            style={{
              fontSize: 64,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 20,
              opacity: titleProgress,
              transform: `translateY(${titleTranslateY}px)`,
              lineHeight: 1.1,
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
              marginBottom: 40,
              opacity: Math.max(0, subtitleProgress),
            }}
          >
            {subtitle}
          </p>

          {/* 功能列表 */}
          <div
            style={{
              display: 'flex',
              flexDirection: 'column',
              gap: 20,
            }}
          >
            {items.slice(0, 3).map((item, index) => {
              const itemProgress = spring({
                frame: frame - 45 - index * 10,
                fps,
                config: { damping: 200, stiffness: 100, mass: 0.5 },
              });

              const itemTranslateX = interpolate(itemProgress, [0, 1], [-20, 0], {
                extrapolateLeft: 'clamp',
                extrapolateRight: 'clamp',
              });

              return (
                <div
                  key={index}
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 16,
                    opacity: Math.max(0, itemProgress),
                    transform: `translateX(${itemTranslateX}px)`,
                  }}
                >
                  <span style={{ fontSize: 32 }}>{item.icon || '✓'}</span>
                  <div>
                    <h3
                      style={{
                        fontSize: 22,
                        fontWeight: 'bold',
                        color: colors.text,
                        margin: 0,
                      }}
                    >
                      {item.title}
                    </h3>
                    {item.description && (
                      <p
                        style={{
                          fontSize: 16,
                          color: 'rgba(255, 255, 255, 0.7)',
                          margin: 0,
                        }}
                      >
                        {item.description}
                      </p>
                    )}
                  </div>
                </div>
              );
            })}
          </div>

          {/* CTA */}
          {cta && (
            <div
              style={{
                marginTop: 40,
                opacity: Math.max(
                  0,
                  spring({
                    frame: frame - 75,
                    fps,
                    config: { damping: 200, stiffness: 100, mass: 0.5 },
                  })
                ),
              }}
            >
              <div
                style={{
                  display: 'inline-block',
                  padding: '16px 40px',
                  backgroundColor: colors.primary,
                  borderRadius: 8,
                  fontSize: 20,
                  fontWeight: 'bold',
                  color: colors.text,
                  boxShadow: `0 8px 24px ${colors.primary}40`,
                }}
              >
                {cta.text}
              </div>
            </div>
          )}
        </div>

        {/* 右侧 - 产品图片 */}
        <div
          style={{
            flex: 1,
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
          }}
        >
          {productImage ? (
            <Img
              src={staticFile(productImage)}
              style={{
                maxWidth: '100%',
                maxHeight: '80%',
                borderRadius: 20,
                boxShadow: '0 30px 60px rgba(0, 0, 0, 0.4)',
                opacity: productProgress,
                transform: `scale(${productScale})`,
              }}
            />
          ) : (
            <div
              style={{
                width: 500,
                height: 400,
                borderRadius: 20,
                background: `linear-gradient(135deg, ${colors.primary} 0%, ${colors.secondary} 100%)`,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                boxShadow: '0 30px 60px rgba(0, 0, 0, 0.4)',
                opacity: productProgress,
                transform: `scale(${productScale})`,
              }}
            >
              <span
                style={{
                  fontSize: 100,
                  color: 'rgba(255, 255, 255, 0.3)',
                }}
              >
                📦
              </span>
            </div>
          )}
        </div>
      </div>
    </AbsoluteFill>
  );
};

/**
 * 默认 props
 */
export const defaultProps: ProductMarketingProps = {
  title: '全新产品发布',
  subtitle: '创新引领未来',
  colors: DEFAULT_COLORS,
  items: [
    {
      title: '极速性能',
      description: '满足你所有需求的极致性能',
      icon: '⚡',
    },
    {
      title: '智能设计',
      description: '直觉式交互，开箱即用',
      icon: '🎯',
    },
    {
      title: '安全可靠',
      description: '企业级安全防护，内置保障',
      icon: '🔒',
    },
  ],
  cta: {
    text: '立即体验',
  },
};

export default ProductMarketing;
