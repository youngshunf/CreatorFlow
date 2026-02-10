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
 * Composition code for step-by-step tutorial video
 * Features numbered steps with clear progression
 */
const STEP_BY_STEP_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Sequence,
} from 'remotion';

interface TutorialStep {
  title: string;
  description?: string;
  icon?: string;
}

interface StepByStepTutorialProps {
  title?: string;
  subtitle?: string;
  items?: TutorialStep[];
  colors?: {
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;
}

export const StepByStepTutorial: React.FC<StepByStepTutorialProps> = ({
  title = 'How To Guide',
  subtitle = 'Follow these simple steps',
  items = [
    { title: 'Step 1', description: 'First, do this', icon: '1️⃣' },
    { title: 'Step 2', description: 'Then, do that', icon: '2️⃣' },
    { title: 'Step 3', description: 'Finally, complete', icon: '3️⃣' },
  ],
  colors = {
    primary: '#10b981',
    secondary: '#059669',
    background: '#0f172a',
    text: '#f8fafc',
  },
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Calculate timing for each step
  const introFrames = 60; // 2 seconds intro
  const stepFrames = Math.floor((durationInFrames - introFrames - 30) / items.length);
  
  // Current step based on frame
  const currentStepIndex = Math.min(
    Math.floor((frame - introFrames) / stepFrames),
    items.length - 1
  );
  const isIntro = frame < introFrames;

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
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
      {/* Header */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          padding: '40px 60px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          borderBottom: \`1px solid \${colors.primary}30\`,
        }}
      >
        {/* Logo */}
        {logo ? (
          <img
            src={logo}
            alt="Logo"
            style={{
              height: 40,
              opacity: titleProgress,
            }}
          />
        ) : (
          <div style={{ width: 40 }} />
        )}

        {/* Progress indicator */}
        <div
          style={{
            display: 'flex',
            gap: 8,
          }}
        >
          {items.map((_, index) => (
            <div
              key={index}
              style={{
                width: 40,
                height: 4,
                borderRadius: 2,
                backgroundColor:
                  index <= currentStepIndex && !isIntro
                    ? colors.primary
                    : \`\${colors.primary}30\`,
                transition: 'background-color 0.3s',
              }}
            />
          ))}
        </div>
      </div>

      {/* Main content */}
      <div
        style={{
          flex: 1,
          display: 'flex',
          flexDirection: 'column',
          justifyContent: 'center',
          alignItems: 'center',
          padding: 80,
          paddingTop: 120,
        }}
      >
        {/* Intro section */}
        {isIntro && (
          <div
            style={{
              textAlign: 'center',
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
              }}
            >
              {title}
            </h1>
            <p
              style={{
                fontSize: 32,
                color: colors.primary,
                margin: 0,
                opacity: titleProgress,
              }}
            >
              {subtitle}
            </p>
          </div>
        )}

        {/* Step content */}
        {!isIntro && items[currentStepIndex] && (
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: 60,
              maxWidth: 1200,
            }}
          >
            {/* Step number */}
            <div
              style={{
                width: 150,
                height: 150,
                borderRadius: '50%',
                backgroundColor: colors.primary,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                fontSize: 64,
                fontWeight: 'bold',
                color: colors.text,
                boxShadow: \`0 20px 40px \${colors.primary}40\`,
                opacity: spring({
                  frame: frame - introFrames - currentStepIndex * stepFrames,
                  fps,
                  config: { damping: 200, stiffness: 100, mass: 0.5 },
                }),
                transform: \`scale(\${spring({
                  frame: frame - introFrames - currentStepIndex * stepFrames,
                  fps,
                  config: { damping: 200, stiffness: 100, mass: 0.5 },
                  from: 0.5,
                  to: 1,
                })})\`,
              }}
            >
              {currentStepIndex + 1}
            </div>

            {/* Step content */}
            <div
              style={{
                flex: 1,
              }}
            >
              <h2
                style={{
                  fontSize: 48,
                  fontWeight: 'bold',
                  color: colors.text,
                  margin: 0,
                  marginBottom: 16,
                  opacity: spring({
                    frame: frame - introFrames - currentStepIndex * stepFrames - 5,
                    fps,
                    config: { damping: 200, stiffness: 100, mass: 0.5 },
                  }),
                }}
              >
                {items[currentStepIndex].title}
              </h2>
              {items[currentStepIndex].description && (
                <p
                  style={{
                    fontSize: 28,
                    color: 'rgba(255, 255, 255, 0.7)',
                    margin: 0,
                    lineHeight: 1.5,
                    opacity: spring({
                      frame: frame - introFrames - currentStepIndex * stepFrames - 10,
                      fps,
                      config: { damping: 200, stiffness: 100, mass: 0.5 },
                    }),
                  }}
                >
                  {items[currentStepIndex].description}
                </p>
              )}
            </div>
          </div>
        )}
      </div>
    </AbsoluteFill>
  );
};

export default StepByStepTutorial;
`;

/**
 * Composition code for explainer/educational video
 * Clean design for explaining concepts
 */
const EXPLAINER_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

interface ExplainerProps {
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
}

export const Explainer: React.FC<ExplainerProps> = ({
  title = 'Understanding the Concept',
  subtitle = 'A simple explanation',
  items = [
    { title: 'Point One', description: 'Key insight here', icon: '💡' },
    { title: 'Point Two', description: 'Another important fact', icon: '📊' },
    { title: 'Point Three', description: 'Final takeaway', icon: '🎯' },
  ],
  colors = {
    primary: '#6366f1',
    secondary: '#8b5cf6',
    background: '#ffffff',
    text: '#1f2937',
  },
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
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
        padding: 80,
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
            left: 60,
            height: 40,
            opacity: titleProgress,
          }}
        />
      )}

      {/* Header */}
      <div
        style={{
          marginBottom: 60,
          paddingTop: 40,
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
          }}
        >
          {title}
        </h1>
        <p
          style={{
            fontSize: 28,
            color: colors.primary,
            margin: 0,
            opacity: titleProgress,
          }}
        >
          {subtitle}
        </p>
      </div>

      {/* Content grid */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(3, 1fr)',
          gap: 40,
          flex: 1,
        }}
      >
        {items.map((item, index) => {
          const itemProgress = spring({
            frame: frame - 30 - index * 15,
            fps,
            config: { damping: 200, stiffness: 100, mass: 0.5 },
          });

          return (
            <div
              key={index}
              style={{
                padding: 40,
                backgroundColor: \`\${colors.primary}08\`,
                borderRadius: 20,
                border: \`2px solid \${colors.primary}20\`,
                opacity: Math.max(0, itemProgress),
                transform: \`translateY(\${interpolate(itemProgress, [0, 1], [30, 0])}px)\`,
              }}
            >
              <span
                style={{
                  fontSize: 48,
                  display: 'block',
                  marginBottom: 20,
                }}
              >
                {item.icon || '📌'}
              </span>
              <h3
                style={{
                  fontSize: 28,
                  fontWeight: 'bold',
                  color: colors.text,
                  margin: 0,
                  marginBottom: 12,
                }}
              >
                {item.title}
              </h3>
              {item.description && (
                <p
                  style={{
                    fontSize: 18,
                    color: \`\${colors.text}99\`,
                    margin: 0,
                    lineHeight: 1.5,
                  }}
                >
                  {item.description}
                </p>
              )}
            </div>
          );
        })}
      </div>
    </AbsoluteFill>
  );
};

export default Explainer;
`;

/**
 * Composition code for tips/listicle video
 * Quick tips format with numbered items
 */
const TIPS_CODE = `
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
} from 'remotion';

interface TipsProps {
  title?: string;
  subtitle?: string;
  items?: Array<{
    title: string;
    description?: string;
  }>;
  colors?: {
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;
}

export const Tips: React.FC<TipsProps> = ({
  title = '5 Quick Tips',
  subtitle = 'Boost your productivity',
  items = [
    { title: 'Tip 1', description: 'Start with the basics' },
    { title: 'Tip 2', description: 'Practice regularly' },
    { title: 'Tip 3', description: 'Learn from mistakes' },
    { title: 'Tip 4', description: 'Stay consistent' },
    { title: 'Tip 5', description: 'Never stop learning' },
  ],
  colors = {
    primary: '#f59e0b',
    secondary: '#d97706',
    background: '#1f2937',
    text: '#ffffff',
  },
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Title animation
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
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
      {/* Background accent */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          right: 0,
          width: '40%',
          height: '100%',
          background: \`linear-gradient(180deg, \${colors.primary}20 0%, transparent 100%)\`,
        }}
      />

      {/* Content */}
      <div
        style={{
          display: 'flex',
          height: '100%',
          padding: 80,
        }}
      >
        {/* Left side - Title */}
        <div
          style={{
            flex: 1,
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center',
          }}
        >
          {logo && (
            <img
              src={logo}
              alt="Logo"
              style={{
                height: 40,
                marginBottom: 40,
                opacity: titleProgress,
              }}
            />
          )}
          <h1
            style={{
              fontSize: 64,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 16,
              opacity: titleProgress,
            }}
          >
            {title}
          </h1>
          <p
            style={{
              fontSize: 28,
              color: colors.primary,
              margin: 0,
              opacity: titleProgress,
            }}
          >
            {subtitle}
          </p>
        </div>

        {/* Right side - Tips list */}
        <div
          style={{
            flex: 1,
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center',
            gap: 20,
          }}
        >
          {items.slice(0, 5).map((item, index) => {
            const itemProgress = spring({
              frame: frame - 20 - index * 12,
              fps,
              config: { damping: 200, stiffness: 100, mass: 0.5 },
            });

            return (
              <div
                key={index}
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: 20,
                  opacity: Math.max(0, itemProgress),
                  transform: \`translateX(\${interpolate(itemProgress, [0, 1], [30, 0])}px)\`,
                }}
              >
                <div
                  style={{
                    width: 50,
                    height: 50,
                    borderRadius: 12,
                    backgroundColor: colors.primary,
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    fontSize: 24,
                    fontWeight: 'bold',
                    color: colors.background,
                    flexShrink: 0,
                  }}
                >
                  {index + 1}
                </div>
                <div>
                  <h3
                    style={{
                      fontSize: 24,
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
                        color: 'rgba(255, 255, 255, 0.6)',
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
      </div>
    </AbsoluteFill>
  );
};

export default Tips;
`;

/**
 * Step-by-Step Tutorial Template (16:9)
 * Clear progression through numbered steps
 */
export const stepByStepTemplate: VideoTemplate = {
  id: 'tutorial-steps',
  name: '分步教程',
  description:
    '带有清晰步骤进度的教程模板。适用于操作指南、演练和教学内容。',
  category: 'tutorial',
  defaultConfig: {
    width: ASPECT_RATIOS.HORIZONTAL.width,
    height: ASPECT_RATIOS.HORIZONTAL.height,
    fps: 30,
    durationInFrames: 450, // 15 seconds
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
  compositionCode: STEP_BY_STEP_CODE,
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
    title: '概念解析',
    subtitle: '简单易懂的讲解',
    colors: DEFAULT_COLOR_SCHEMES.minimal,
    items: [
      {
        title: '核心要点一',
        description: '你需要理解的基础概念',
        icon: '💡',
      },
      {
        title: '核心要点二',
        description: '在基础上深入了解更多细节',
        icon: '📊',
      },
      {
        title: '核心要点三',
        description: '实际应用与关键收获',
        icon: '🎯',
      },
    ],
    animationStyle: 'spring',
  },
  compositionCode: EXPLAINER_CODE,
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
  compositionCode: TIPS_CODE,
  aspectRatio: '16:9',
  tags: ['技巧', '清单', '快速', '建议', '效率'],
  useCases: [
    '快速技巧视频',
    '清单类内容',
    '经验汇总',
    '最佳实践',
    '生活妙招',
  ],
};

/**
 * Tutorial Template with variants
 */
export const tutorialTemplate: VideoTemplateWithVariants = {
  ...stepByStepTemplate,
  id: 'tutorial',
  name: '教程',
  description:
    '教育类模板，适用于教程、讲解和教学内容。清晰专业的设计，提升学习效果。',
  variants: [
    {
      id: 'steps',
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
      name: 'Tips & Listicle',
      aspectRatio: 'HORIZONTAL',
      config: tipsTemplate.defaultConfig,
    },
  ],
};

/**
 * Get tutorial template by variant
 */
export function getTutorialTemplate(
  variant: 'steps' | 'explainer' | 'tips' = 'steps'
): VideoTemplate {
  switch (variant) {
    case 'explainer':
      return explainerTemplate;
    case 'tips':
      return tipsTemplate;
    case 'steps':
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
