/**
 * URL validation for MCP servers
 *
 * Simple URL validation without AI - checks basic URL format.
 */

import { debug } from '../utils/debug.ts';
import type { AgentError } from '../agent/errors.ts';

export interface UrlValidationResult {
  valid: boolean;
  /** Simple error message for validation failures */
  error?: string;
  /** Typed error for API/billing failures - display as ErrorBanner */
  typedError?: AgentError;
}

/**
 * Validate a MCP server URL
 * Simple validation without AI - checks basic URL format
 */
export async function validateMcpUrl(
  url: string,
  _apiKey?: string,
  _oauthToken?: string,
): Promise<UrlValidationResult> {
  debug('[url-validator] Validating URL:', url);

  try {
    // Basic URL validation
    const parsed = new URL(url);

    // Check protocol
    if (parsed.protocol !== 'https:' && parsed.protocol !== 'http:') {
      return { valid: false, error: 'URL must use http:// or https:// protocol' };
    }

    // Check for credentials in URL (security risk)
    if (parsed.username || parsed.password) {
      return { valid: false, error: 'URL should not contain credentials' };
    }

    // URL is valid
    return { valid: true };
  } catch (err) {
    debug('[url-validator] Error:', err);
    return { valid: false, error: 'Invalid URL format' };
  }
}
