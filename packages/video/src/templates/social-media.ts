/**
 * Social Media Templates
 *
 * Pre-configured templates optimized for social media platforms.
 * Includes 9:16 vertical (TikTok, Reels, Shorts) and 1:1 square (Instagram, Facebook) formats.
 *
 * @requirements 5.1, 5.5
 */

import type {
  VideoTemplate,
  VideoTemplateWithVariants,
  TemplateVariant,
} from './types';
import {
  DEFAULT_COLOR_SCHEMES,
  ASPECT_RATIOS,
} from './types';

/**
 * Composition code for vertical social media video (9:16)
 * Uses TitleAnimation with spring animations optimized for mobile viewing
 */
const VERTICAL_COMPOSITION_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

interface SocialMediaVerticalProps {
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
}

export const SocialMediaVertical: React.FC<SocialMediaVerticalProps> = ({
  title = 'Your Title Here',
  subtitle = 'Swipe up to learn more',
  colors = {
    primary: '#ec4899',
    secondary: '#8b5cf6',
    background: '#0f0f23',
    text: '#ffffff',
  },
  logo,
  cta,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Title animation with spring
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

  // Subtitle animation (delayed)
  const subtitleProgress = spring({
    frame: frame - 15,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  // CTA animation (delayed more)
  const ctaProgress = spring({
    frame: frame - 30,
    fps,
    config: { damping: 200, stiffness: 80, mass: 1 },
  });

  // Fade out at the end
  const fadeOut = interpolate(
    frame,
    [durationInFrames - 20, durationInFrames],
    [1, 0],
    { extrapolateLeft: 'clamp', extrapolateRight: 'clamp' }
  );

  // Gradient background animation
  const gradientRotation = interpolate(frame, [0, durationInFrames], [0, 360]);

  return (
    <AbsoluteFill
      style={{
        background: \`linear-gradient(\${gradientRotation}deg, \${colors.background} 0%, \${colors.primary}20 50%, \${colors.background} 100%)\`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
        padding: 40,
      }}
    >
      {/* Logo at top */}
      {logo && (
        <img
          src={logo}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 60,
            height: 60,
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
          gap: 30,
        }}
      >
        {/* Title */}
        <h1
          style={{
            fontSize: 64,
            fontWeight: 'bold',
            color: colors.text,
            margin: 0,
            opacity: titleProgress,
            transform: \`translateY(\${titleTranslateY}px)\`,
            lineHeight: 1.2,
            maxWidth: '90%',
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
            opacity: Math.max(0, subtitleProgress),
            maxWidth: '80%',
          }}
        >
          {subtitle}
        </p>

        {/* CTA Button */}
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
              transform: \`scale(\${interpolate(ctaProgress, [0, 1], [0.8, 1])})\`,
              boxShadow: \`0 10px 30px \${colors.primary}60\`,
            }}
          >
            {cta.text}
          </div>
        )}
      </div>

      {/* Decorative elements */}
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

export default SocialMediaVertical;
`;

/**
 * Composition code for square social media video (1:1)
 * Optimized for Instagram feed and Facebook posts
 */
const SQUARE_COMPOSITION_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

interface SocialMediaSquareProps {
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
}

export const SocialMediaSquare: React.FC<SocialMediaSquareProps> = ({
  title = 'Your Title Here',
  subtitle = 'Double tap to like',
  colors = {
    primary: '#ec4899',
    secondary: '#8b5cf6',
    background: '#0f0f23',
    text: '#ffffff',
  },
  logo,
  cta,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
  });

  const titleScale = interpolate(titleProgress, [0, 1], [0.9, 1]);

  // Subtitle animation
  const subtitleProgress = spring({
    frame: frame - 15,
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
        background: \`radial-gradient(circle at center, \${colors.primary}30 0%, \${colors.background} 70%)\`,
        justifyContent: 'center',
        alignItems: 'center',
        fontFamily: 'system-ui, -apple-system, sans-serif',
        opacity: fadeOut,
        padding: 60,
      }}
    >
      {/* Logo */}
      {logo && (
        <img
          src={logo}
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

      {/* Content */}
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
            transform: \`scale(\${titleScale})\`,
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

export default SocialMediaSquare;
`;

/**
 * Social Media Vertical Template (9:16)
 * Optimized for TikTok, Instagram Reels, YouTube Shorts
 */
export const socialMediaVerticalTemplate: VideoTemplate = {
  id: 'social-media-vertical',
  name: '社交媒体竖版',
  description:
    '竖版视频模板，适用于抖音、Instagram Reels 和 YouTube Shorts。包含吸睛动画和移动端优先设计。',
  category: 'social-media',
  defaultConfig: {
    width: ASPECT_RATIOS.VERTICAL.width,
    height: ASPECT_RATIOS.VERTICAL.height,
    fps: 30,
    durationInFrames: 150, // 5 seconds
  },
  defaultProps: {
    title: '你的精彩标题',
    subtitle: '上滑了解更多',
    colors: DEFAULT_COLOR_SCHEMES.vibrant,
    animationStyle: 'spring',
    cta: {
      text: '了解更多',
    },
  },
  compositionCode: VERTICAL_COMPOSITION_CODE,
  aspectRatio: '9:16',
  tags: ['抖音', '快手', '短视频', '竖屏', '移动端', '社交'],
  useCases: [
    '抖音短视频',
    '快手短视频',
    '微信视频号',
    'B站竖屏视频',
    '移动端优先内容',
  ],
};

/**
 * Social Media Square Template (1:1)
 * Optimized for Instagram Feed, Facebook Posts
 */
export const socialMediaSquareTemplate: VideoTemplate = {
  id: 'social-media-square',
  name: '社交媒体方形',
  description:
    '方形视频模板，适用于 Instagram 动态和 Facebook。简洁设计，内容居中展示。',
  category: 'social-media',
  defaultConfig: {
    width: ASPECT_RATIOS.SQUARE.width,
    height: ASPECT_RATIOS.SQUARE.height,
    fps: 30,
    durationInFrames: 150, // 5 seconds
  },
  defaultProps: {
    title: '你的精彩标题',
    subtitle: '双击点赞',
    colors: DEFAULT_COLOR_SCHEMES.vibrant,
    animationStyle: 'spring',
    cta: {
      text: '立即购买',
    },
  },
  compositionCode: SQUARE_COMPOSITION_CODE,
  aspectRatio: '1:1',
  tags: ['小红书', '微信', '方形', '动态', '社交'],
  useCases: [
    '小红书图文动态',
    '微信朋友圈',
    '微博动态',
    'B站动态',
  ],
};

/**
 * Social Media Template with all variants
 * Provides both vertical and square options
 */
export const socialMediaTemplate: VideoTemplateWithVariants = {
  ...socialMediaVerticalTemplate,
  id: 'social-media',
  name: '社交媒体',
  description:
    '多比例社交媒体模板，适用于跨平台内容创作。',
  variants: [
    {
      id: 'vertical',
      name: 'Vertical (9:16)',
      aspectRatio: 'VERTICAL',
      config: socialMediaVerticalTemplate.defaultConfig,
    },
    {
      id: 'square',
      name: 'Square (1:1)',
      aspectRatio: 'SQUARE',
      config: socialMediaSquareTemplate.defaultConfig,
    },
    {
      id: 'portrait',
      name: 'Portrait (4:5)',
      aspectRatio: 'PORTRAIT',
      config: {
        width: ASPECT_RATIOS.PORTRAIT.width,
        height: ASPECT_RATIOS.PORTRAIT.height,
        fps: 30,
        durationInFrames: 150,
      },
    },
  ],
};

/**
 * Get social media template by variant
 */
export function getSocialMediaTemplate(
  variant: 'vertical' | 'square' = 'vertical'
): VideoTemplate {
  switch (variant) {
    case 'square':
      return socialMediaSquareTemplate;
    case 'vertical':
    default:
      return socialMediaVerticalTemplate;
  }
}

/**
 * All social media templates
 */
export const SOCIAL_MEDIA_TEMPLATES = [
  socialMediaVerticalTemplate,
  socialMediaSquareTemplate,
] as const;
