/**
 * @sprouty-ai/video
 *
 * Video creation engine for Sprouty AI - Remotion integration,
 * compositions, templates, and rendering capabilities.
 */

// Re-export types
export * from "./types";

// Re-export Root component
export { RemotionRoot } from "./Root";

// Re-export compositions
export * from "./compositions";

// Re-export templates
export * from "./templates";

// Re-export reusable components (Transition renamed to avoid conflict with Transition type)
export {
  AnimatedText,
  type AnimatedTextProps,
  type AnimationType,
  Transition as TransitionComponent,
  TransitionOverlay,
  type TransitionProps,
  type TransitionOverlayProps,
  type TransitionEffectType,
  Background,
  GradientBackground,
  ImageBackground,
  type BackgroundProps,
  type GradientBackgroundProps,
  type ImageBackgroundProps,
  type BackgroundType,
  type GradientDirection,
  type BackgroundAnimation,
  Logo,
  AnimatedLogo,
  type LogoProps,
  type AnimatedLogoProps,
  type LogoAnimation,
  type LogoPosition,
  TitleText,
  type TitleTextProps,
  SubtitleText,
  type SubtitleTextProps,
  SlideItem,
  type SlideItemProps,
  BarItem,
  type BarItemProps,
} from "./components";

// Re-export hooks
export * from "./hooks";

// Re-export skills for AI Agent
export * from "./skills";

// Re-export utility functions
export * from "./utils";
