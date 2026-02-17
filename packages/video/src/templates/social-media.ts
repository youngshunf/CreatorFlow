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
  compositionId: 'SocialMediaVertical',
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
  compositionId: 'SocialMediaSquare',
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
