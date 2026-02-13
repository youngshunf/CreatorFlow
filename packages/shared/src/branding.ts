/**
 * Centralized branding assets for Sprouty AI
 * Used by OAuth callback pages
 */

export const SPROUTY_LOGO = [
  '███████ ██████  ██████   ██████  ██    ██ ████████ ██    ██',
  '██      ██   ██ ██   ██ ██    ██ ██    ██    ██     ██  ██ ',
  '███████ ██████  ██████  ██    ██ ██    ██    ██      ████  ',
  '     ██ ██      ██   ██ ██    ██ ██    ██    ██       ██   ',
  '███████ ██      ██   ██  ██████   ██████     ██       ██   ',
] as const;

/** Logo as a single string for HTML templates */
export const SPROUTY_LOGO_HTML = SPROUTY_LOGO.map((line) => line.trimEnd()).join('\n');

/** Viewer URL for session sharing */
export const VIEWER_URL = 'https://viewer.creatorflow.app';

// Legacy export for backward compatibility
export const CRAFT_LOGO = SPROUTY_LOGO;
export const CRAFT_LOGO_HTML = SPROUTY_LOGO_HTML;
export const CREATOR_FLOW_LOGO = SPROUTY_LOGO;
export const CREATOR_FLOW_LOGO_HTML = SPROUTY_LOGO_HTML;
