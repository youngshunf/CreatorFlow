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
 * Product Marketing Template (16:9)
 * Professional product showcase with feature highlights
 */
export const productMarketingTemplate: VideoTemplate = {
  id: 'marketing-product',
  name: '产品营销',
  description:
    '专业的产品展示模板，适用于产品发布和功能介绍。横版设计，适合 YouTube 和演示文稿。',
  category: 'marketing',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 300, // 10 seconds
  },
  defaultProps: {
    title: '全新产品发布',
    subtitle: '创新引领未来',
    colors: DEFAULT_COLOR_SCHEMES.corporate,
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
    animationStyle: 'spring',
  },
  compositionId: 'ProductMarketing',
  aspectRatio: '16:9',
  tags: ['产品', '营销', '专业', '商务', '发布'],
  useCases: [
    '产品发布',
    '功能演示',
    '企业宣传',
    '信息流广告',
    '官网首页视频',
  ],
};

/**
 * Promotional Ad Template (16:9)
 * Bold, attention-grabbing design for sales and promotions
 */
export const promoAdTemplate: VideoTemplate = {
  id: 'marketing-promo',
  name: '促销广告',
  description:
    '醒目的促销视频，适用于打折、特卖和限时优惠。大胆设计，吸引眼球。',
  category: 'marketing',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 180, // 6 seconds
  },
  defaultProps: {
    title: '限时特惠',
    subtitle: '全场5折起',
    colors: DEFAULT_COLOR_SCHEMES.playful,
    cta: {
      text: '立即抢购',
    },
    animationStyle: 'spring',
  },
  compositionId: 'PromoAd',
  aspectRatio: '16:9',
  tags: ['促销', '广告', '特卖', '优惠', '营销'],
  useCases: [
    '促销活动',
    '限时优惠',
    '节日营销',
    '清仓甩卖',
    '新品上市',
  ],
};

/**
 * Marketing Template with all variants
 */
export const marketingTemplate: VideoTemplateWithVariants = {
  ...productMarketingTemplate,
  id: 'marketing',
  name: '营销推广',
  description:
    '多场景营销模板，适用于产品推广和促销活动。',
  variants: [
    {
      id: 'product',
      name: 'Product Marketing',
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
