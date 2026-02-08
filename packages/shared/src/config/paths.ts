/**
 * Centralized path configuration for Sprouty AI (智小芽).
 *
 * Supports multi-instance development via SPROUTY_CONFIG_DIR environment variable.
 * When running from a numbered folder (e.g., sprouty-ai-1), the detect-instance.sh
 * script sets SPROUTY_CONFIG_DIR to ~/.sprouty-ai-1, allowing multiple instances to run
 * simultaneously with separate configurations.
 *
 * Default (non-numbered folders): ~/.sprouty-ai/
 * Instance 1 (-1 suffix): ~/.sprouty-ai-1/
 * Instance 2 (-2 suffix): ~/.sprouty-ai-2/
 */

import { homedir } from 'os';
import { join } from 'path';
import { existsSync, renameSync } from 'fs';

const DEFAULT_CONFIG_DIR = join(homedir(), '.sprouty-ai');
const LEGACY_CONFIG_DIR = join(homedir(), '.creator-flow');

/**
 * Auto-migrate legacy config directory (~/.creator-flow/ → ~/.sprouty-ai/).
 * Only migrates if the old directory exists and the new one does not.
 */
function migrateConfigDir(): string {
  if (!existsSync(DEFAULT_CONFIG_DIR) && existsSync(LEGACY_CONFIG_DIR)) {
    try {
      renameSync(LEGACY_CONFIG_DIR, DEFAULT_CONFIG_DIR);
      console.log(`[sprouty] Migrated config: ${LEGACY_CONFIG_DIR} → ${DEFAULT_CONFIG_DIR}`);
    } catch (err) {
      console.error(`[sprouty] Failed to migrate config dir:`, err);
      // Fall back to legacy dir if migration fails
      return LEGACY_CONFIG_DIR;
    }
  }
  return DEFAULT_CONFIG_DIR;
}

// Allow override via environment variable for multi-instance dev
// Falls back to default ~/.sprouty-ai/ for production and non-numbered dev folders
// Supports legacy CREATOR_FLOW_CONFIG_DIR for backward compatibility
export const CONFIG_DIR =
  process.env.SPROUTY_CONFIG_DIR ||
  process.env.CREATOR_FLOW_CONFIG_DIR ||
  migrateConfigDir();
