/**
 * SceneComposer - 场景编排器核心组件
 *
 * 使用 TransitionSeries 编排多个场景片段，支持 5 种过渡效果。
 * 这是视频创作的核心组合，Agent 通过配置 scenes + transitions 来编排视频。
 */
import React from "react";
import { AbsoluteFill, useVideoConfig } from "remotion";
import { TransitionSeries, linearTiming } from "@remotion/transitions";
import { fade } from "@remotion/transitions/fade";
import { slide } from "@remotion/transitions/slide";
import { wipe } from "@remotion/transitions/wipe";
import { flip } from "@remotion/transitions/flip";
import { clockWipe } from "@remotion/transitions/clock-wipe";
import { z } from "zod";
import type { CalculateMetadataFunction } from "remotion";
import type { TransitionPresentation } from "@remotion/transitions";

// 导入所有组合组件
import { TitleAnimation } from "./TitleAnimation";
import { Slideshow } from "./Slideshow";
import { DataVisualization } from "./DataVisualization";
import { ProductShowcase } from "./ProductShowcase";
import { SocialMediaVertical } from "./SocialMediaVertical";
import { SocialMediaSquare } from "./SocialMediaSquare";
import { StepByStepTutorial } from "./StepByStepTutorial";
import { Explainer } from "./Explainer";
import { Tips } from "./Tips";
import { ProductMarketing } from "./ProductMarketing";
import { PromoAd } from "./PromoAd";

// ============================================================================
// 组合映射表
// ============================================================================

/**
 * compositionId → React 组件的映射
 * 所有 11 个内置组合都在此注册
 */
const COMPOSITION_MAP: Record<string, React.FC<any>> = {
  TitleAnimation,
  Slideshow,
  DataVisualization,
  ProductShowcase,
  SocialMediaVertical,
  SocialMediaSquare,
  StepByStepTutorial,
  Explainer,
  Tips,
  ProductMarketing,
  PromoAd,
};

// ============================================================================
// Zod Schema
// ============================================================================

const TransitionDirectionSchema = z.enum([
  "from-left",
  "from-right",
  "from-top",
  "from-bottom",
]);

const TransitionConfigSchema = z.object({
  type: z.enum(["fade", "slide", "wipe", "flip", "clock-wipe", "none"]),
  durationInFrames: z.number().int().positive(),
  direction: TransitionDirectionSchema.optional(),
});

const SceneItemSchema = z.object({
  id: z.string(),
  name: z.string(),
  compositionId: z.string(),
  durationInFrames: z.number().int().positive(),
  props: z.record(z.string(), z.any()).default({}),
});

export const SceneComposerSchema = z.object({
  scenes: z.array(SceneItemSchema),
  transitions: z.array(TransitionConfigSchema),
});

// ============================================================================
// 类型
// ============================================================================

type TransitionConfig = z.infer<typeof TransitionConfigSchema>;
export type SceneComposerProps = z.infer<typeof SceneComposerSchema>;

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * 根据过渡配置获取 Remotion presentation 对象
 */
const getPresentation = (
  transition: TransitionConfig,
  videoSize: { width: number; height: number },
): TransitionPresentation<Record<string, unknown>> => {
  switch (transition.type) {
    case "fade":
      return fade() as TransitionPresentation<Record<string, unknown>>;
    case "slide":
      return slide({
        direction: transition.direction ?? "from-right",
      }) as TransitionPresentation<Record<string, unknown>>;
    case "wipe":
      return wipe({
        direction: transition.direction ?? "from-left",
      }) as TransitionPresentation<Record<string, unknown>>;
    case "flip":
      return flip({
        direction: transition.direction ?? "from-left",
      }) as TransitionPresentation<Record<string, unknown>>;
    case "clock-wipe":
      return clockWipe({
        width: videoSize.width,
        height: videoSize.height,
      }) as unknown as TransitionPresentation<Record<string, unknown>>;
    default:
      return fade() as TransitionPresentation<Record<string, unknown>>;
  }
};

// ============================================================================
// SceneRenderer - 根据 compositionId 渲染对应组件
// ============================================================================

interface SceneRendererProps {
  compositionId: string;
  props: Record<string, unknown>;
}

const SceneRenderer: React.FC<SceneRendererProps> = ({
  compositionId,
  props,
}) => {
  const Component = COMPOSITION_MAP[compositionId];

  if (!Component) {
    return (
      <AbsoluteFill
        style={{
          backgroundColor: "#1a1a2e",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          color: "#ef4444",
          fontSize: 24,
          fontFamily: "sans-serif",
        }}
      >
        未知组合: {compositionId}
      </AbsoluteFill>
    );
  }

  return <Component {...props} />;
};

// ============================================================================
// SceneComposer 主组件
// ============================================================================

/**
 * SceneComposer - 使用 TransitionSeries 编排多个场景
 *
 * 接收 scenes 和 transitions 数组，按顺序渲染每个场景，
 * 并在场景之间插入过渡效果。
 */
export const SceneComposer: React.FC<SceneComposerProps> = ({
  scenes,
  transitions,
}) => {
  const { width, height } = useVideoConfig();

  if (!scenes || scenes.length === 0) {
    return (
      <AbsoluteFill
        style={{
          backgroundColor: "#1a1a2e",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          color: "#94a3b8",
          fontSize: 24,
          fontFamily: "sans-serif",
        }}
      >
        暂无场景
      </AbsoluteFill>
    );
  }

  return (
    <TransitionSeries>
      {scenes.map((scene, i) => (
        <React.Fragment key={scene.id}>
          <TransitionSeries.Sequence durationInFrames={scene.durationInFrames}>
            <SceneRenderer
              compositionId={scene.compositionId}
              props={scene.props}
            />
          </TransitionSeries.Sequence>
          {transitions[i] && transitions[i].type !== "none" && (
            <TransitionSeries.Transition
              presentation={getPresentation(transitions[i], { width, height })}
              timing={linearTiming({
                durationInFrames: transitions[i].durationInFrames,
              })}
            />
          )}
        </React.Fragment>
      ))}
    </TransitionSeries>
  );
};

// ============================================================================
// calculateMetadata - 动态计算总时长
// ============================================================================

/**
 * 动态计算 SceneComposer 的总时长
 * 总时长 = sum(scenes.durationInFrames) - sum(transitions.durationInFrames)
 * 因为过渡效果会让相邻场景重叠
 */
export const calculateSceneComposerMetadata: CalculateMetadataFunction<
  SceneComposerProps
> = async ({ props }) => {
  const scenesTotal = props.scenes.reduce(
    (sum, s) => sum + s.durationInFrames,
    0,
  );
  const transitionsTotal = props.transitions
    .filter((t) => t.type !== "none")
    .reduce((sum, t) => sum + t.durationInFrames, 0);

  return {
    durationInFrames: Math.max(1, scenesTotal - transitionsTotal),
  };
};
