/**
 * Video Compositions
 *
 * This module exports all available video composition components.
 * Compositions are the main building blocks for creating videos in Remotion.
 */

// Export all composition components
export { TitleAnimation } from './TitleAnimation';
export type { TitleAnimationProps } from './TitleAnimation';

export { Slideshow } from './Slideshow';
export type { SlideshowProps } from './Slideshow';

export { DataVisualization } from './DataVisualization';
export type { DataVisualizationProps, DataPoint } from './DataVisualization';

export { ProductShowcase } from './ProductShowcase';
export type { ProductShowcaseProps, FeatureItem } from './ProductShowcase';

/**
 * List of all available composition IDs
 */
export const COMPOSITION_IDS = [
  'TitleAnimation',
  'Slideshow',
  'DataVisualization',
  'ProductShowcase',
] as const;

export type CompositionId = (typeof COMPOSITION_IDS)[number];
