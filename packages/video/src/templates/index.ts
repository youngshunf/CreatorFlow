/**
 * Video Templates
 *
 * This module exports all available video templates.
 * Templates provide pre-configured compositions for common use cases.
 *
 * @requirements 5.1, 5.2, 5.3, 5.4, 5.5
 */

// Export template types
export * from './types';

// Export social media templates
export {
  socialMediaTemplate,
  socialMediaVerticalTemplate,
  socialMediaSquareTemplate,
  getSocialMediaTemplate,
  SOCIAL_MEDIA_TEMPLATES,
} from './social-media';

// Export marketing templates
export {
  marketingTemplate,
  productMarketingTemplate,
  promoAdTemplate,
  getMarketingTemplate,
  MARKETING_TEMPLATES,
} from './marketing';

// Export tutorial templates
export {
  tutorialTemplate,
  stepByStepTemplate,
  explainerTemplate,
  tipsTemplate,
  getTutorialTemplate,
  TUTORIAL_TEMPLATES,
} from './tutorial';

import type { VideoTemplate, TemplateCategory } from './types';
import { SOCIAL_MEDIA_TEMPLATES } from './social-media';
import { MARKETING_TEMPLATES } from './marketing';
import { TUTORIAL_TEMPLATES } from './tutorial';

/**
 * All available templates
 */
export const ALL_TEMPLATES: readonly VideoTemplate[] = [
  ...SOCIAL_MEDIA_TEMPLATES,
  ...MARKETING_TEMPLATES,
  ...TUTORIAL_TEMPLATES,
];

/**
 * Templates grouped by category
 */
export const TEMPLATES_BY_CATEGORY: Record<TemplateCategory, readonly VideoTemplate[]> = {
  'social-media': SOCIAL_MEDIA_TEMPLATES,
  'marketing': MARKETING_TEMPLATES,
  'tutorial': TUTORIAL_TEMPLATES,
};

/**
 * Get template by ID
 * @param id Template ID
 * @returns VideoTemplate or undefined if not found
 */
export function getTemplateById(id: string): VideoTemplate | undefined {
  return ALL_TEMPLATES.find((template) => template.id === id);
}

/**
 * Get templates by category
 * @param category Template category
 * @returns Array of templates in the category
 */
export function getTemplatesByCategory(category: TemplateCategory): readonly VideoTemplate[] {
  return TEMPLATES_BY_CATEGORY[category] || [];
}

/**
 * Get template by ID (throws if not found)
 * @param id Template ID
 * @returns VideoTemplate
 * @throws Error if template not found
 */
export function getTemplate(id: string): VideoTemplate {
  const template = getTemplateById(id);
  if (!template) {
    throw new Error(`Template not found: ${id}`);
  }
  return template;
}

/**
 * Search templates by tags
 * @param tags Tags to search for
 * @returns Array of matching templates
 */
export function searchTemplatesByTags(tags: string[]): VideoTemplate[] {
  const lowerTags = tags.map((t) => t.toLowerCase());
  return ALL_TEMPLATES.filter((template) =>
    template.tags?.some((tag) => lowerTags.includes(tag.toLowerCase()))
  );
}

/**
 * Get all template IDs
 */
export const TEMPLATE_IDS = ALL_TEMPLATES.map((t) => t.id);

/**
 * Template ID type
 */
export type TemplateId = (typeof TEMPLATE_IDS)[number];
