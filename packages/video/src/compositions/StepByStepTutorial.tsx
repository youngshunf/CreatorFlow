/**
 * StepByStepTutorial Composition
 *
 * 分步教程视频组合，带有清晰步骤进度。
 * 适用于操作指南、演练和教学内容。
 *
 * @requirements 5.3, 5.5
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
 * 教程步骤接口
 */
export interface TutorialStep {
  title: string;
  description?: string;
  icon?: string;
}

/**
 * Props for StepByStepTutorial composition
 */
export interface StepByStepTutorialProps extends TemplateProps {
  /** 主标题 */
  title?: string;
  /** 副标题 */
  subtitle?: string;
  /** 步骤列表 */
  items?: TutorialStep[];
  /** Logo 图片路径 */
  logo?: string;
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#10b981',
  secondary: '#059669',
  background: '#0f172a',
  text: '#f8fafc',
};

/**
 * StepByStepTutorial - 分步教程视频
 *
 * 清晰的步骤进度展示，适合教学内容
 */
export const StepByStepTutorial: React.FC<StepByStepTutorialProps> = ({
  title = '操作指南',
  subtitle = '跟着这些简单步骤操作',
  items = [
    { title: '步骤 1', description: '首先，做这个', icon: '1️⃣' },
    { title: '步骤 2', description: '然后，做那个', icon: '2️⃣' },
    { title: '步骤 3', description: '最后，完成', icon: '3️⃣' },
  ],
  colors = DEFAULT_COLORS,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // 计算每个步骤的时间
  const introFrames = 60; // 2 秒介绍
  const stepFrames = Math.floor((durationInFrames - introFrames - 30) / items.length);

  // 基于帧数计算当前步骤
  const currentStepIndex = Math.min(
    Math.floor((frame - introFrames) / stepFrames),
    items.length - 1
  );
  const isIntro = frame < introFrames;

  // 标题动画
  const titleProgress = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100, mass: 0.5 },
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
      {/* 头部 */}
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
          borderBottom: `1px solid ${colors.primary}30`,
        }}
      >
        {/* Logo */}
        {logo ? (
          <Img
            src={staticFile(logo)}
            alt="Logo"
            style={{
              height: 40,
              opacity: titleProgress,
            }}
          />
        ) : (
          <div style={{ width: 40 }} />
        )}

        {/* 进度指示器 */}
        <div
          style={{
            display: 'flex',
            gap: 8,
          }}
        >
          {items.map((_, index) => {
            const isActive = index <= currentStepIndex && !isIntro;
            const barOpacity = isActive ? 1 : 0.3;

            return (
              <div
                key={index}
                style={{
                  width: 40,
                  height: 4,
                  borderRadius: 2,
                  backgroundColor: colors.primary,
                  opacity: barOpacity,
                }}
              />
            );
          })}
        </div>
      </div>

      {/* 主内容 */}
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
        {/* 介绍部分 */}
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

        {/* 步骤内容 */}
        {!isIntro && items[currentStepIndex] && (
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: 60,
              maxWidth: 1200,
            }}
          >
            {/* 步骤编号 */}
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
                boxShadow: `0 20px 40px ${colors.primary}40`,
                opacity: spring({
                  frame: frame - introFrames - currentStepIndex * stepFrames,
                  fps,
                  config: { damping: 200, stiffness: 100, mass: 0.5 },
                }),
                transform: `scale(${spring({
                  frame: frame - introFrames - currentStepIndex * stepFrames,
                  fps,
                  config: { damping: 200, stiffness: 100, mass: 0.5 },
                  from: 0.5,
                  to: 1,
                })})`,
              }}
            >
              {currentStepIndex + 1}
            </div>

            {/* 步骤内容 */}
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

/**
 * 默认 props
 */
export const defaultProps: StepByStepTutorialProps = {
  title: '操作指南',
  subtitle: '跟着这些简单步骤操作',
  colors: DEFAULT_COLORS,
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
};

export default StepByStepTutorial;
