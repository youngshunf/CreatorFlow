/**
 * Tips Composition
 *
 * 技巧清单视频组合，快速技巧模板，适用于清单式内容。
 * 以引人入胜的格式分享多个见解。
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
 * 技巧项接口
 */
export interface TipItem {
  title: string;
  description?: string;
}

/**
 * Props for Tips composition
 */
export interface TipsProps extends TemplateProps {
  /** 主标题 */
  title?: string;
  /** 副标题 */
  subtitle?: string;
  /** 技巧列表 */
  items?: TipItem[];
  /** Logo 图片路径 */
  logo?: string;
}

/**
 * 默认颜色配置
 */
const DEFAULT_COLORS = {
  primary: '#f59e0b',
  secondary: '#d97706',
  background: '#1f2937',
  text: '#ffffff',
};

/**
 * Tips - 技巧清单视频
 *
 * 快速技巧格式，带编号项目
 */
export const Tips: React.FC<TipsProps> = ({
  title = '5 个实用技巧',
  subtitle = '今天就提升你的效率',
  items = [
    { title: '技巧 1', description: '从基础开始' },
    { title: '技巧 2', description: '定期练习' },
    { title: '技巧 3', description: '从错误中学习' },
    { title: '技巧 4', description: '保持一致' },
    { title: '技巧 5', description: '永不停止学习' },
  ],
  colors = DEFAULT_COLORS,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

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
      {/* 背景装饰 */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          right: 0,
          width: '40%',
          height: '100%',
          background: `linear-gradient(180deg, ${colors.primary}20 0%, transparent 100%)`,
        }}
      />

      {/* 内容 */}
      <div
        style={{
          display: 'flex',
          height: '100%',
          padding: 80,
        }}
      >
        {/* 左侧 - 标题 */}
        <div
          style={{
            flex: 1,
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center',
          }}
        >
          {logo && (
            <Img
              src={staticFile(logo)}
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

        {/* 右侧 - 技巧列表 */}
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

            const itemTranslateX = interpolate(itemProgress, [0, 1], [30, 0], {
              extrapolateLeft: 'clamp',
              extrapolateRight: 'clamp',
            });

            return (
              <div
                key={index}
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: 20,
                  opacity: Math.max(0, itemProgress),
                  transform: `translateX(${itemTranslateX}px)`,
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

/**
 * 默认 props
 */
export const defaultProps: TipsProps = {
  title: '5 个实用技巧',
  subtitle: '今天就提升你的效率',
  colors: DEFAULT_COLORS,
  items: [
    { title: '早起行动', description: '带着目标开始新的一天' },
    { title: '保持专注', description: '排除干扰，集中精力' },
    { title: '适当休息', description: '休息是为了保持充沛精力' },
    { title: '回顾进展', description: '追踪你的成就和进步' },
    { title: '持续学习', description: '永远不要停止成长' },
  ],
};

export default Tips;
