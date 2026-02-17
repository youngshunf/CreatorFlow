/**
 * Video Compositions
 *
 * This module exports all available video composition components.
 * Compositions are the main building blocks for creating videos in Remotion.
 */

// 场景编排器（核心组合）
export {
  SceneComposer,
  SceneComposerSchema,
  calculateSceneComposerMetadata,
} from "./SceneComposer";
export type { SceneComposerProps } from "./SceneComposer";

// Export all composition components
export { TitleAnimation } from "./TitleAnimation";
export type { TitleAnimationProps } from "./TitleAnimation";

export { Slideshow } from "./Slideshow";
export type { SlideshowProps } from "./Slideshow";

export { DataVisualization } from "./DataVisualization";
export type { DataVisualizationProps, DataPoint } from "./DataVisualization";

export { ProductShowcase } from "./ProductShowcase";
export type { ProductShowcaseProps, FeatureItem } from "./ProductShowcase";

// 社交媒体组合
export { SocialMediaVertical } from "./SocialMediaVertical";
export type { SocialMediaVerticalProps } from "./SocialMediaVertical";

export { SocialMediaSquare } from "./SocialMediaSquare";
export type { SocialMediaSquareProps } from "./SocialMediaSquare";

// 教程组合
export { StepByStepTutorial } from "./StepByStepTutorial";
export type {
  StepByStepTutorialProps,
  TutorialStep,
} from "./StepByStepTutorial";

export { Explainer } from "./Explainer";
export type { ExplainerProps, ExplainerItem } from "./Explainer";

export { Tips } from "./Tips";
export type { TipsProps, TipItem } from "./Tips";

// 营销组合
export { ProductMarketing } from "./ProductMarketing";
export type { ProductMarketingProps, ProductFeature } from "./ProductMarketing";

export { PromoAd } from "./PromoAd";
export type { PromoAdProps } from "./PromoAd";

/**
 * List of all available composition IDs
 */
export const COMPOSITION_IDS = [
  "TitleAnimation",
  "Slideshow",
  "DataVisualization",
  "ProductShowcase",
  "SocialMediaVertical",
  "SocialMediaSquare",
  "StepByStepTutorial",
  "Explainer",
  "Tips",
  "ProductMarketing",
  "PromoAd",
] as const;

export type CompositionId = (typeof COMPOSITION_IDS)[number];
