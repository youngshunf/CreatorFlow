/**
 * Template Types
 *
 * Type definitions for video templates and their props.
 * These interfaces define the standard props that all compositions
 * and templates should accept.
 *
 * @requirements 5.4
 */

/**
 * Color scheme interface for templates
 */
export interface ColorSchemeConfig {
  primary: string;
  secondary: string;
  background: string;
  text: string;
}

/**
 * Standard template props interface
 *
 * This interface defines the common props that all video compositions
 * and templates should accept for consistency and reusability.
 */
export interface TemplateProps {
  /** Main title text */
  title?: string;
  /** Subtitle or tagline */
  subtitle?: string;
  /** Array of items (slides, features, data points, etc.) */
  items?: Array<{
    title: string;
    description?: string;
    /** Image path relative to public/ directory (e.g., 'images/slide1.jpg') */
    image?: string;
    icon?: string;
  }>;
  /** Color scheme */
  colors?: ColorSchemeConfig;
  /** Logo image path relative to public/ directory (e.g., 'images/logo.png') */
  logo?: string;
  /** Call-to-action configuration */
  cta?: {
    text: string;
    url?: string;
  };
  /** Font family override */
  fontFamily?: string;
  /** Animation style */
  animationStyle?: 'fade' | 'slide' | 'scale' | 'spring' | 'zoom';
}

/**
 * Video configuration for templates
 */
export interface VideoTemplateConfig {
  /** Video width in pixels */
  width: number;
  /** Video height in pixels */
  height: number;
  /** Frames per second */
  fps: number;
  /** Total duration in frames */
  durationInFrames: number;
}

/**
 * Video template definition
 *
 * Defines a complete video template including metadata,
 * default configuration, and composition code.
 */
export interface VideoTemplate {
  /** Unique template identifier */
  id: string;
  /** Display name of the template */
  name: string;
  /** Template description */
  description: string;
  /** Template category */
  category: 'social-media' | 'marketing' | 'tutorial';
  /** Default video configuration */
  defaultConfig: VideoTemplateConfig;
  /** Default props for the template */
  defaultProps: TemplateProps;
  /** Composition component code or reference */
  compositionCode: string;
  /** Preview thumbnail URL (optional) */
  thumbnail?: string;
  /** Tags for searchability */
  tags?: string[];
  /** Aspect ratio label (e.g., '9:16', '16:9', '1:1') */
  aspectRatio?: string;
  /** Recommended use cases */
  useCases?: string[];
}

/**
 * Template category type
 */
export type TemplateCategory = VideoTemplate['category'];

/**
 * Default color schemes for templates
 */
export const DEFAULT_COLOR_SCHEMES = {
  modern: {
    primary: '#6366f1',
    secondary: '#8b5cf6',
    background: '#1a1a2e',
    text: '#ffffff',
  },
  minimal: {
    primary: '#000000',
    secondary: '#666666',
    background: '#ffffff',
    text: '#000000',
  },
  playful: {
    primary: '#f43f5e',
    secondary: '#fb923c',
    background: '#fef3c7',
    text: '#1f2937',
  },
  corporate: {
    primary: '#0ea5e9',
    secondary: '#0284c7',
    background: '#0f172a',
    text: '#f8fafc',
  },
  cinematic: {
    primary: '#fbbf24',
    secondary: '#f59e0b',
    background: '#000000',
    text: '#ffffff',
  },
  vibrant: {
    primary: '#ec4899',
    secondary: '#8b5cf6',
    background: '#0f0f23',
    text: '#ffffff',
  },
  nature: {
    primary: '#10b981',
    secondary: '#059669',
    background: '#064e3b',
    text: '#ecfdf5',
  },
} as const;

/**
 * Color scheme type
 */
export type ColorScheme = keyof typeof DEFAULT_COLOR_SCHEMES;

/**
 * Get color scheme by name
 */
export function getColorScheme(name: ColorScheme): ColorSchemeConfig {
  return DEFAULT_COLOR_SCHEMES[name];
}

/**
 * Common aspect ratios for video templates
 */
export const ASPECT_RATIOS = {
  /** 9:16 vertical (TikTok, Instagram Reels, YouTube Shorts) */
  VERTICAL: { width: 1080, height: 1920, label: '9:16' },
  /** 1:1 square (Instagram Feed, Facebook) */
  SQUARE: { width: 1080, height: 1080, label: '1:1' },
  /** 16:9 horizontal (YouTube, presentations) */
  HORIZONTAL: { width: 1920, height: 1080, label: '16:9' },
  /** 4:5 portrait (Instagram Feed optimal) */
  PORTRAIT: { width: 1080, height: 1350, label: '4:5' },
} as const;

/**
 * Aspect ratio type
 */
export type AspectRatio = keyof typeof ASPECT_RATIOS;

/**
 * Get aspect ratio configuration
 */
export function getAspectRatio(ratio: AspectRatio) {
  return ASPECT_RATIOS[ratio];
}

/**
 * Template variant for different aspect ratios
 */
export interface TemplateVariant {
  /** Variant identifier */
  id: string;
  /** Variant name */
  name: string;
  /** Aspect ratio */
  aspectRatio: AspectRatio;
  /** Video configuration */
  config: VideoTemplateConfig;
}

/**
 * Extended template with variants
 */
export interface VideoTemplateWithVariants extends VideoTemplate {
  /** Available variants for different aspect ratios */
  variants?: TemplateVariant[];
}
