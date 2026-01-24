/**
 * Centralized path configuration for CreatorFlow.
 *
 * Supports multi-instance development via CREATOR_FLOW_CONFIG_DIR environment variable.
 * When running from a numbered folder (e.g., creator-flow-1), the detect-instance.sh
 * script sets CREATOR_FLOW_CONFIG_DIR to ~/.creator-flow-1, allowing multiple instances to run
 * simultaneously with separate configurations.
 *
 * Default (non-numbered folders): ~/.creator-flow/
 * Instance 1 (-1 suffix): ~/.creator-flow-1/
 * Instance 2 (-2 suffix): ~/.creator-flow-2/
 */

import { homedir } from 'os';
import { join } from 'path';

// Allow override via environment variable for multi-instance dev
// Falls back to default ~/.creator-flow/ for production and non-numbered dev folders
export const CONFIG_DIR = process.env.CREATOR_FLOW_CONFIG_DIR || join(homedir(), '.creator-flow');
