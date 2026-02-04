/**
 * ProductShowcase Composition
 *
 * A video composition that displays product features with animated
 * showcase effects. Perfect for marketing and promotional videos.
 *
 * @requirements 3.4, 3.5, 3.6
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Img,
  Sequence,
} from 'remotion';
import type { TemplateProps } from '../templates/types';

/**
 * Feature item interface
 */
export interface FeatureItem {
  title: string;
  description?: string;
  icon?: string;
  image?: string;
}

/**
 * Props for ProductShowcase composition
 */
export interface ProductShowcaseProps extends TemplateProps {
  /** Product name */
  title?: string;
  /** Product tagline */
  subtitle?: string;
  /** Product features */
  items?: FeatureItem[];
  /** Product image URL */
  productImage?: string;
  /** Layout style */
  layout?: 'centered' | 'split' | 'features-grid';
  /** Animation style */
  animationStyle?: 'fade' | 'slide' | 'zoom' | 'spring';
}

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
 * Default features for demo purposes
 */
const DEFAULT_FEATURES: FeatureItem[] = [
  {
    title: 'Feature One',
    description: 'Amazing feature that helps you achieve more',
    icon: '🚀',
  },
  {
    title: 'Feature Two',
    description: 'Powerful capability for better results',
    icon: '⚡',
  },
  {
    title: 'Feature Three',
    description: 'Smart solution for modern challenges',
    icon: '🎯',
  },
  {
    title: 'Feature Four',
    description: 'Innovative approach to common problems',
    icon: '💡',
  },
];

/**
 * ProductShowcase - Displays product features with animations
 *
 * This composition creates professional product showcase videos
 * with customizable layouts and animation effects.
 */
export const ProductShowcase: React.FC<ProductShowcaseProps> = ({
  title = 'Amazing Product',
  subtitle = 'The future of productivity',
  items = DEFAULT_FEATURES,
  colors = DEFAULT_COLORS,
  productImage,
  logo,
  cta,
  layout = 'centered',
  animationStyle = 'spring',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Calculate animation phases
  const introPhase = 60; // 2 seconds for intro
  const featurePhase = 60; // 2 seconds per feature
  const ctaPhase = 60; // 2 seconds for CTA

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const titleTranslateY = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
    from: -40,
    to: 0,
  });

  // Subtitle animation (delayed)
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // Product image animation
  const productImageProgress = spring({
    frame: frame - 30,
    fps,
    config: { damping: 200, stiffness: 80, mass: 1 },
  });

  const productImageScale = interpolate(productImageProgress, [0, 1], [0.8, 1]);

  // Fade out at the end
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 30, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  // Render feature item
  const renderFeature = (feature: FeatureItem, index: number) => {
    const delay = introPhase + index * 15;
    const featureProgress = spring({
      frame: frame - delay,
      fps,
      config: { damping: 200, stiffness: 100, mass: 0.5 },
    });

    const translateX =
      animationStyle === 'slide'
        ? interpolate(featureProgress, [0, 1], [50, 0])
        : 0;

    const translateY =
      animationStyle === 'spring'
        ? spring({
            frame: frame - delay,
            fps,
            config: { damping: 200, stiffness: 100, mass: 0.5 },
            from: 30,
            to: 0,
          })
        : 0;

    const scale =
      animationStyle === 'zoom'
        ? interpolate(featureProgress, [0, 1], [0.8, 1])
        : 1;

    return (
      <div
        key={index}
        style={{
          display: 'flex',
          alignItems: 'flex-start',
          gap: 20,
          padding: 20,
          backgroundColor: 'rgba(255, 255, 255, 0.05)',
          borderRadius: 16,
          opacity: featureProgress,
          transform: `translateX(${translateX}px) translateY(${translateY}px) scale(${scale})`,
        }}
      >
        {/* Icon or Image */}
        {feature.icon ? (
          <span style={{ fontSize: 40 }}>{feature.icon}</span>
        ) : feature.image ? (
          <Img
            src={feature.image}
            style={{
              width: 60,
              height: 60,
              borderRadius: 12,
              objectFit: 'cover',
            }}
          />
        ) : (
          <div
            style={{
              width: 60,
              height: 60,
              borderRadius: 12,
              backgroundColor: colors.primary,
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              fontSize: 24,
              color: colors.text,
              fontWeight: 'bold',
            }}
          >
            {index + 1}
          </div>
        )}

        {/* Text content */}
        <div style={{ flex: 1 }}>
          <h3
            style={{
              fontSize: 24,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 8,
            }}
          >
            {feature.title}
          </h3>
          {feature.description && (
            <p
              style={{
                fontSize: 16,
                color: 'rgba(255, 255, 255, 0.7)',
                margin: 0,
                lineHeight: 1.5,
              }}
            >
              {feature.description}
            </p>
          )}
        </div>
      </div>
    );
  };

  // Render CTA button
  const renderCTA = () => {
    if (!cta) return null;

    const ctaDelay = introPhase + items.length * 15 + 30;
    const ctaProgress = spring({
      frame: frame - ctaDelay,
      fps,
      config: { damping: 200, stiffness: 100, mass: 0.5 },
    });

    const ctaScale = interpolate(ctaProgress, [0, 1], [0.9, 1]);

    return (
      <div
        style={{
          marginTop: 40,
          opacity: ctaProgress,
          transform: `scale(${ctaScale})`,
        }}
      >
        <div
          style={{
            display: 'inline-block',
            padding: '16px 40px',
            backgroundColor: colors.primary,
            borderRadius: 12,
            fontSize: 24,
            fontWeight: 'bold',
            color: colors.text,
            boxShadow: `0 8px 24px ${colors.primary}40`,
          }}
        >
          {cta.text}
        </div>
      </div>
    );
  };

  // Render centered layout
  const renderCenteredLayout = () => (
    <AbsoluteFill
      style={{
        justifyContent: 'center',
        alignItems: 'center',
        padding: 80,
      }}
    >
      {/* Header */}
      <div
        style={{
          textAlign: 'center',
          marginBottom: 60,
        }}
      >
        <h1
          style={{
            fontSize: 72,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            marginBottom: 20,
            opacity: titleProgress,
            transform: `translateY(${titleTranslateY}px)`,
          }}
        >
          {title}
        </h1>
        <p
          style={{
            fontSize: 32,
            color: colors.primary,
            margin: 0,
            opacity: Math.max(0, subtitleProgress),
          }}
        >
          {subtitle}
        </p>
      </div>

      {/* Product Image */}
      {productImage && (
        <div
          style={{
            marginBottom: 60,
            opacity: productImageProgress,
            transform: `scale(${productImageScale})`,
          }}
        >
          <Img
            src={productImage}
            style={{
              maxWidth: 600,
              maxHeight: 400,
              borderRadius: 20,
              boxShadow: '0 20px 60px rgba(0, 0, 0, 0.4)',
            }}
          />
        </div>
      )}

      {/* Features Grid */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(2, 1fr)',
          gap: 20,
          maxWidth: 900,
        }}
      >
        {items.map((feature, index) => renderFeature(feature, index))}
      </div>

      {/* CTA */}
      {renderCTA()}
    </AbsoluteFill>
  );

  // Render split layout
  const renderSplitLayout = () => (
    <AbsoluteFill
      style={{
        flexDirection: 'row',
      }}
    >
      {/* Left side - Product Image */}
      <div
        style={{
          flex: 1,
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center',
          backgroundColor: `${colors.primary}20`,
          padding: 60,
        }}
      >
        {productImage ? (
          <Img
            src={productImage}
            style={{
              maxWidth: '80%',
              maxHeight: '80%',
              borderRadius: 20,
              boxShadow: '0 20px 60px rgba(0, 0, 0, 0.4)',
              opacity: productImageProgress,
              transform: `scale(${productImageScale})`,
            }}
          />
        ) : (
          <div
            style={{
              width: 400,
              height: 400,
              borderRadius: 20,
              background: `linear-gradient(135deg, ${colors.primary} 0%, ${colors.secondary} 100%)`,
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              opacity: productImageProgress,
              transform: `scale(${productImageScale})`,
            }}
          >
            <span
              style={{
                fontSize: 120,
                color: 'rgba(255, 255, 255, 0.3)',
              }}
            >
              📦
            </span>
          </div>
        )}
      </div>

      {/* Right side - Content */}
      <div
        style={{
          flex: 1,
          display: 'flex',
          flexDirection: 'column',
          justifyContent: 'center',
          padding: 60,
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
            transform: `translateY(${titleTranslateY}px)`,
          }}
        >
          {title}
        </h1>
        <p
          style={{
            fontSize: 24,
            color: colors.primary,
            margin: 0,
            marginBottom: 40,
            opacity: Math.max(0, subtitleProgress),
          }}
        >
          {subtitle}
        </p>

        {/* Features List */}
        <div
          style={{
            display: 'flex',
            flexDirection: 'column',
            gap: 16,
          }}
        >
          {items.slice(0, 4).map((feature, index) => renderFeature(feature, index))}
        </div>

        {/* CTA */}
        {renderCTA()}
      </div>
    </AbsoluteFill>
  );

  // Render features grid layout
  const renderFeaturesGridLayout = () => (
    <AbsoluteFill
      style={{
        padding: 80,
      }}
    >
      {/* Header */}
      <div
        style={{
          textAlign: 'center',
          marginBottom: 60,
        }}
      >
        <h1
          style={{
            fontSize: 64,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            marginBottom: 16,
            opacity: titleProgress,
            transform: `translateY(${titleTranslateY}px)`,
          }}
        >
          {title}
        </h1>
        <p
          style={{
            fontSize: 28,
            color: colors.primary,
            margin: 0,
            opacity: Math.max(0, subtitleProgress),
          }}
        >
          {subtitle}
        </p>
      </div>

      {/* Features Grid */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(2, 1fr)',
          gap: 30,
          flex: 1,
        }}
      >
        {items.map((feature, index) => {
          const delay = introPhase + index * 15;
          const featureProgress = spring({
            frame: frame - delay,
            fps,
            config: { damping: 200, stiffness: 100, mass: 0.5 },
          });

          return (
            <div
              key={index}
              style={{
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                justifyContent: 'center',
                padding: 40,
                backgroundColor: 'rgba(255, 255, 255, 0.05)',
                borderRadius: 20,
                textAlign: 'center',
                opacity: featureProgress,
                transform: `scale(${interpolate(featureProgress, [0, 1], [0.9, 1])})`,
              }}
            >
              {feature.icon ? (
                <span style={{ fontSize: 60, marginBottom: 20 }}>
                  {feature.icon}
                </span>
              ) : (
                <div
                  style={{
                    width: 80,
                    height: 80,
                    borderRadius: 20,
                    backgroundColor: colors.primary,
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    fontSize: 32,
                    color: colors.text,
                    fontWeight: 'bold',
                    marginBottom: 20,
                  }}
                >
                  {index + 1}
                </div>
              )}
              <h3
                style={{
                  fontSize: 28,
                  fontWeight: 'bold',
                  color: colors.text,
                  margin: 0,
                  marginBottom: 12,
                }}
              >
                {feature.title}
              </h3>
              {feature.description && (
                <p
                  style={{
                    fontSize: 18,
                    color: 'rgba(255, 255, 255, 0.7)',
                    margin: 0,
                    lineHeight: 1.5,
                  }}
                >
                  {feature.description}
                </p>
              )}
            </div>
          );
        })}
      </div>

      {/* CTA */}
      <div style={{ textAlign: 'center' }}>{renderCTA()}</div>
    </AbsoluteFill>
  );

  // Render layout based on type
  const renderLayout = () => {
    switch (layout) {
      case 'split':
        return renderSplitLayout();
      case 'features-grid':
        return renderFeaturesGridLayout();
      case 'centered':
      default:
        return renderCenteredLayout();
    }
  };

  return (
    <AbsoluteFill
      style={{
        backgroundColor: colors.background,
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
      }}
    >
      {renderLayout()}

      {/* Optional Logo */}
      {logo && (
        <img
          src={logo}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 40,
            right: 40,
            height: 50,
            opacity: titleProgress,
          }}
        />
      )}
    </AbsoluteFill>
  );
};

export default ProductShowcase;
