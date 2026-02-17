/**
 * Reusable Video Components
 *
 * This module exports reusable UI components for video compositions.
 * These components provide common animation effects and visual elements.
 *
 * @requirements 4.1, 4.2, 4.3, 4.4, 4.5
 */

// AnimatedText - Text with various animation effects
export {
  AnimatedText,
  type AnimatedTextProps,
  type AnimationType,
} from "./AnimatedText";

// Transition - Scene transition effects
export {
  Transition,
  TransitionOverlay,
  type TransitionProps,
  type TransitionOverlayProps,
  type TransitionEffectType,
} from "./Transition";

// Background - Various background types with animations
export {
  Background,
  GradientBackground,
  ImageBackground,
  type BackgroundProps,
  type GradientBackgroundProps,
  type ImageBackgroundProps,
  type BackgroundType,
  type GradientDirection,
  type BackgroundAnimation,
} from "./Background";

// Logo - Logo display with animation effects
export {
  Logo,
  AnimatedLogo,
  type LogoProps,
  type AnimatedLogoProps,
  type LogoAnimation,
  type LogoPosition,
} from "./Logo";

// TitleText - Animated title component for Sequence-based compositions
export { TitleText, type TitleTextProps } from "./TitleText";

// SubtitleText - Animated subtitle component for Sequence-based compositions
export { SubtitleText, type SubtitleTextProps } from "./SubtitleText";

// SlideItem - Animated slide component for slideshow compositions
export { SlideItem, type SlideItemProps } from "./SlideItem";

// BarItem - Animated bar chart item for data visualization
export { BarItem, type BarItemProps } from "./BarItem";
