/**
 * Tutorial Templates
 *
 * Pre-configured templates for educational and tutorial videos.
 * Optimized for clear presentation of step-by-step instructions.
 *
 * @requirements 5.3, 5.5
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
 * Step-by-Step Tutorial Template (16:9)
 * Clear progression through numbered steps
 */
export const stepByStepTemplate: VideoTemplate = {
  id: 'tutorial-step-by-step',
  name: '分步教程',
  description:
    '分步教程模板，适用于操作指南和说明。清晰的步骤展示，易于跟随。',
  category: 'tutorial',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 360, // 12 seconds
  },
  defaultProps: {
    title: '操作指南',
    subtitle: '跟着这些简单步骤操作',
    colors: DEFAULT_COLOR_SCHEMES.nature,
    items: [
      {
        title: '准备工作环境',
        description: '收集所有必要材料并搭建环境',
        icon: '1️⃣',
      },
      {
        title: '按步骤执行',
        description: '仔细且有条理地执行每个步骤',
        icon: '2️⃣',
      },
      {
        title: '检查与优化',
        description: '检查成果并进行必要的调整',
        icon: '3️⃣',
      },
      {
        title: '完成与分享',
        description: '完善项目并分享你的成果',
        icon: '4️⃣',
      },
    ],
    animationStyle: 'spring',
  },
  compositionId: 'StepByStepTutorial',
  aspectRatio: '16:9',
  tags: ['教程', '操作指南', '步骤', '指导', '教学', '说明'],
  useCases: [
    '操作指南',
    '软件教程',
    'DIY 手工教学',
    '流程演示',
    '培训视频',
  ],
};

/**
 * Explainer Template (16:9)
 * Clean design for explaining concepts
 */
export const explainerTemplate: VideoTemplate = {
  id: 'tutorial-explainer',
  name: '概念讲解',
  description:
    '简洁的讲解模板，适用于教育内容。将复杂主题拆解为易懂的要点。',
  category: 'tutorial',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 300, // 10 seconds
  },
  defaultProps: {
    title: '理解核心概念',
    subtitle: '简单易懂的讲解',
    colors: DEFAULT_COLOR_SCHEMES.minimal,
    items: [
      {
        title: '核心要点一',
        description: '基础概念与定义',
        icon: '💡',
      },
      {
        title: '核心要点二',
        description: '工作原理与机制',
        icon: '⚙️',
      },
      {
        title: '核心要点三',
        description: '实际应用与关键收获',
        icon: '🎯',
      },
    ],
    animationStyle: 'spring',
  },
  compositionId: 'Explainer',
  aspectRatio: '16:9',
  tags: ['讲解', '教育', '概念', '学习', '演示'],
  useCases: [
    '概念讲解',
    '教育内容',
    '课程素材',
    '演示文稿',
    '知识分享',
  ],
};

/**
 * Tips/Listicle Template (16:9)
 * Quick tips format with numbered items
 */
export const tipsTemplate: VideoTemplate = {
  id: 'tutorial-tips',
  name: '技巧清单',
  description:
    '快速技巧模板，适用于清单式内容。以引人入胜的格式分享多个见解。',
  category: 'tutorial',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 240, // 8 seconds
  },
  defaultProps: {
    title: '5 个实用技巧',
    subtitle: '今天就提升你的效率',
    colors: DEFAULT_COLOR_SCHEMES.cinematic,
    items: [
      { title: '早起行动', description: '带着目标开始新的一天' },
      { title: '保持专注', description: '排除干扰，集中精力' },
      { title: '适当休息', description: '休息是为了保持充沛精力' },
      { title: '回顾进展', description: '追踪你的成就和进步' },
      { title: '持续学习', description: '永远不要停止成长' },
    ],
    animationStyle: 'spring',
  },
  compositionId: 'Tips',
  aspectRatio: '16:9',
  tags: ['技巧', '清单', '建议', '快速', '实用'],
  useCases: [
    '技巧分享',
    '清单内容',
    '快速建议',
    '生活窍门',
    '效率提升',
  ],
};

/**
 * Tutorial Template with all variants
 */
export const tutorialTemplate: VideoTemplateWithVariants = {
  ...stepByStepTemplate,
  id: 'tutorial',
  name: '教程',
  description:
    '多场景教程模板，适用于教育和培训内容。',
  variants: [
    {
      id: 'step-by-step',
      name: 'Step-by-Step',
      aspectRatio: 'HORIZONTAL',
      config: stepByStepTemplate.defaultConfig,
    },
    {
      id: 'explainer',
      name: 'Explainer',
      aspectRatio: 'HORIZONTAL',
      config: explainerTemplate.defaultConfig,
    },
    {
      id: 'tips',
      name: 'Tips/Listicle',
      aspectRatio: 'HORIZONTAL',
      config: tipsTemplate.defaultConfig,
    },
  ],
};

/**
 * Get tutorial template by variant
 */
export function getTutorialTemplate(
  variant: 'step-by-step' | 'explainer' | 'tips' = 'step-by-step'
): VideoTemplate {
  switch (variant) {
    case 'explainer':
      return explainerTemplate;
    case 'tips':
      return tipsTemplate;
    case 'step-by-step':
    default:
      return stepByStepTemplate;
  }
}

/**
 * All tutorial templates
 */
export const TUTORIAL_TEMPLATES = [
  stepByStepTemplate,
  explainerTemplate,
  tipsTemplate,
] as const;
