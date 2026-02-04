/**
 * Marketing Templates
 *
 * Pre-configured templates for marketing and promotional videos.
 * Optimized for 16:9 horizontal format suitable for YouTube, presentations, and ads.
 *
 * @requirements 5.2, 5.5
 */

import type {
  VideoTemplate,
  VideoTemplateWithVariants,
} from './types';
import {
  DEFAULT_COLOR_SCHEMES,
  ASPECT_RATIOS,
} from './types';

/**
 * Composition code for product marketing video
 * Features product showcase with animated features list
 */
const PRODUCT_MARKETING_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Img,
} from 'remotion';

interface ProductMarketingProps {
  title?: string;
  subtitle?: string;
  items?: Array<{
    title: string;
    description?: string;
    icon?: string;
  }>;
  colors?: {
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;
  cta?: {
    text: string;
    url?: string;
  };
  productImage?: string;
}

export const ProductMarketing: React.FC<ProductMarketingProps> = ({
  title = 'Introducing Our Product',
  subtitle = 'The future of innovation',
  items = [
    { title: 'Feature One', description: 'Amazing capability', icon: '🚀' },
    { title: 'Feature Two', description: 'Powerful performance', icon: '⚡' },
    { title: 'Feature Three', description: 'Smart design', icon: '🎯' },
  ],
  colors = {
    primary: '#0ea5e9',
    secondary: '#0284c7',
    background: '#0f172a',
    text: '#f8fafc',
  },
  logo,
  cta,
  productImage,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

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
    from: -30,
    to: 0,
  });

  // Subtitle animation
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // Product image animation
  const productProgress = spring({
    frame: frame - 30,
    fps,
    config: { damping: 200, stiffness: 80, mass: 1 },
  });

  // Fade out
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
      {/* Background gradient */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          background: \`radial-gradient(ellipse at top right, \${colors.primary}20 0%, transparent 50%),
                       radial-gradient(ellipse at bottom left, \${colors.secondary}20 0%, transparent 50%)\`,
        }}
      />

      {/* Content container */}
      <div
        style={{
          display: 'flex',
          height: '100%',
          padding: 80,
        }}
      >
        {/* Left side - Text content */}
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
            <img
              src={logo}
              alt="Logo"
              style={{
                height: 50,
                marginBottom: 40,
                opacity: titleProgress,
              }}
            />
          )}

          {/* Title */}
          <h1
            style={{
              fontSize: 64,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 20,
              opacity: titleProgress,
              transform: \`translateY(\${titleTranslateY}px)\`,
              lineHeight: 1.1,
            }}
          >
            {title}
          </h1>

          {/* Subtitle */}
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

          {/* Features */}
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

              return (
                <div
                  key={index}
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 16,
                    opacity: Math.max(0, itemProgress),
                    transform: \`translateX(\${interpolate(itemProgress, [0, 1], [-20, 0])}px)\`,
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
                opacity: Math.max(0, spring({
                  frame: frame - 75,
                  fps,
                  config: { damping: 200, stiffness: 100, mass: 0.5 },
                })),
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
                  boxShadow: \`0 8px 24px \${colors.primary}40\`,
                }}
              >
                {cta.text}
              </div>
            </div>
          )}
        </div>

        {/* Right side - Product image */}
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
              src={productImage}
              style={{
                maxWidth: '100%',
                maxHeight: '80%',
                borderRadius: 20,
                boxShadow: '0 30px 60px rgba(0, 0, 0, 0.4)',
                opacity: productProgress,
                transform: \`scale(\${interpolate(productProgress, [0, 1], [0.9, 1])})\`,
              }}
            />
          ) : (
            <div
              style={{
                width: 500,
                height: 400,
                borderRadius: 20,
                background: \`linear-gradient(135deg, \${colors.primary} 0%, \${colors.secondary} 100%)\`,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                boxShadow: '0 30px 60px rgba(0, 0, 0, 0.4)',
                opacity: productProgress,
                transform: \`scale(\${interpolate(productProgress, [0, 1], [0.9, 1])})\`,
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

export default ProductMarketing;
`;

/**
 * Composition code for promotional/ad video
 * Bold, attention-grabbing design for advertisements
 */
const PROMO_AD_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

interface PromoAdProps {
  title?: string;
  subtitle?: string;
  colors?: {
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;
  cta?: {
    text: string;
    url?: string;
  };
  discount?: string;
}

export const PromoAd: React.FC<PromoAdProps> = ({
  title = 'MEGA SALE',
  subtitle = 'Limited time offer',
  colors = {
    primary: '#f43f5e',
    secondary: '#fb923c',
    background: '#1a1a2e',
    text: '#ffffff',
  },
  logo,
  cta = { text: 'Shop Now' },
  discount = '50% OFF',
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Pulsing animation for discount
  const pulse = Math.sin(frame * 0.15) * 0.05 + 1;

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 100, stiffness: 200, mass: 0.5 },
  });

  // Discount animation
  const discountProgress = spring({
    frame: frame - 10,
    fps,
    config: { damping: 100, stiffness: 200, mass: 0.5 },
  });

  // CTA animation
  const ctaProgress = spring({
    frame: frame - 25,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // Fade out
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 20, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  return (
    <AbsoluteFill
      style={{
        background: \`linear-gradient(135deg, \${colors.background} 0%, \${colors.primary}30 100%)\`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
      }}
    >
      {/* Animated background shapes */}
      <div
        style={{
          position: 'absolute',
          top: -100,
          right: -100,
          width: 400,
          height: 400,
          borderRadius: '50%',
          background: \`radial-gradient(circle, \${colors.primary}40 0%, transparent 70%)\`,
          transform: \`scale(\${pulse})\`,
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
          background: \`radial-gradient(circle, \${colors.secondary}30 0%, transparent 70%)\`,
          transform: \`scale(\${pulse * 1.1})\`,
        }}
      />

      {/* Logo */}
      {logo && (
        <img
          src={logo}
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

      {/* Main content */}
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          textAlign: 'center',
          gap: 20,
        }}
      >
        {/* Discount badge */}
        <div
          style={{
            padding: '12px 30px',
            backgroundColor: colors.primary,
            borderRadius: 50,
            fontSize: 28,
            fontWeight: 'bold',
            color: colors.text,
            opacity: discountProgress,
            transform: \`scale(\${interpolate(discountProgress, [0, 1], [0.5, 1]) * pulse})\`,
            boxShadow: \`0 10px 30px \${colors.primary}60\`,
          }}
        >
          {discount}
        </div>

        {/* Title */}
        <h1
          style={{
            fontSize: 96,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleProgress,
            transform: \`scale(\${interpolate(titleProgress, [0, 1], [0.8, 1])})\`,
            textShadow: \`0 4px 20px \${colors.primary}80\`,
            letterSpacing: '0.05em',
          }}
        >
          {title}
        </h1>

        {/* Subtitle */}
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
            background: \`linear-gradient(135deg, \${colors.primary} 0%, \${colors.secondary} 100%)\`,
            borderRadius: 12,
            fontSize: 28,
            fontWeight: 'bold',
            color: colors.text,
            opacity: ctaProgress,
            transform: \`scale(\${interpolate(ctaProgress, [0, 1], [0.9, 1])})\`,
            boxShadow: \`0 15px 40px \${colors.primary}50\`,
          }}
        >
          {cta.text}
        </div>
      </div>
    </AbsoluteFill>
  );
};

export default PromoAd;
`;

/**
 * Product Marketing Template (16:9)
 * Professional product showcase with features
 */
export const productMarketingTemplate: VideoTemplate = {
  id: 'marketing-product',
  name: 'Product Marketing',
  description:
    'Professional product marketing video with feature highlights. Perfect for product launches, demos, and promotional content.',
  category: 'marketing',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 300, // 10 seconds
  },
  defaultProps: {
    title: 'Introducing Our Product',
    subtitle: 'The future of innovation',
    colors: DEFAULT_COLOR_SCHEMES.corporate,
    items: [
      {
        title: 'Lightning Fast',
        description: 'Blazing performance for all your needs',
        icon: '⚡',
      },
      {
        title: 'Smart Design',
        description: 'Intuitive interface that just works',
        icon: '🎯',
      },
      {
        title: 'Secure & Reliable',
        description: 'Enterprise-grade security built-in',
        icon: '🔒',
      },
    ],
    cta: {
      text: 'Get Started',
    },
    animationStyle: 'spring',
  },
  compositionCode: PRODUCT_MARKETING_CODE,
  aspectRatio: '16:9',
  tags: ['product', 'marketing', 'professional', 'business', 'launch'],
  useCases: [
    'Product launches',
    'Feature demonstrations',
    'Company presentations',
    'YouTube ads',
    'Website hero videos',
  ],
};

/**
 * Promotional Ad Template (16:9)
 * Bold, attention-grabbing design for sales and promotions
 */
export const promoAdTemplate: VideoTemplate = {
  id: 'marketing-promo',
  name: 'Promotional Ad',
  description:
    'Eye-catching promotional video for sales, discounts, and special offers. Bold design that grabs attention.',
  category: 'marketing',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 150, // 5 seconds
  },
  defaultProps: {
    title: 'MEGA SALE',
    subtitle: 'Limited time offer - Don\'t miss out!',
    colors: DEFAULT_COLOR_SCHEMES.playful,
    cta: {
      text: 'Shop Now',
    },
    animationStyle: 'spring',
  },
  compositionCode: PROMO_AD_CODE,
  aspectRatio: '16:9',
  tags: ['promo', 'sale', 'discount', 'ad', 'marketing', 'bold'],
  useCases: [
    'Sales promotions',
    'Discount announcements',
    'Flash sales',
    'Holiday campaigns',
    'Social media ads',
  ],
};

/**
 * Marketing Template with variants
 */
export const marketingTemplate: VideoTemplateWithVariants = {
  ...productMarketingTemplate,
  id: 'marketing',
  name: 'Marketing',
  description:
    'Professional marketing templates for product showcases, promotions, and business presentations.',
  variants: [
    {
      id: 'product',
      name: 'Product Showcase',
      aspectRatio: 'HORIZONTAL',
      config: productMarketingTemplate.defaultConfig,
    },
    {
      id: 'promo',
      name: 'Promotional Ad',
      aspectRatio: 'HORIZONTAL',
      config: promoAdTemplate.defaultConfig,
    },
  ],
};

/**
 * Get marketing template by variant
 */
export function getMarketingTemplate(
  variant: 'product' | 'promo' = 'product'
): VideoTemplate {
  switch (variant) {
    case 'promo':
      return promoAdTemplate;
    case 'product':
    default:
      return productMarketingTemplate;
  }
}

/**
 * All marketing templates
 */
export const MARKETING_TEMPLATES = [
  productMarketingTemplate,
  promoAdTemplate,
] as const;
