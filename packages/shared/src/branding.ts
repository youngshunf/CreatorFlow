/**
 * Centralized branding assets for CreatorFlow
 * Used by OAuth callback pages
 */

export const CREATOR_FLOW_LOGO = [
  ' █████  ██████  ██████ ██████ ████████ ██████  ██████ ██      ██████  ██  ██  ██',
  '██      ██   ██ ██     ██  ██    ██ ██   ██ ██   ██ ██      ██   ██ ██  ██  ██',
  '██      ██████  ████   █████  ██   ██ ██████  ████   ██      ██   ██ ██ ████ ██',
  '██      ██   ██ ██     ██  ██    ██ ██   ██ ██     ██   ██ ██   ██ ███  ████',
  ' █████  ██   ██ ██████ ██  ██    ██ ██   ██ ██     ████  ██████  ██    ██',
] as const;

/** Logo as a single string for HTML templates */
export const CREATOR_FLOW_LOGO_HTML = CREATOR_FLOW_LOGO.map((line) => line.trimEnd()).join('\n');

/** Viewer URL for session sharing */
export const VIEWER_URL = 'https://viewer.creatorflow.app';

// Legacy export for backward compatibility
export const CRAFT_LOGO = CREATOR_FLOW_LOGO;
export const CRAFT_LOGO_HTML = CREATOR_FLOW_LOGO_HTML;
