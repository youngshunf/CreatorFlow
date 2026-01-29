/**
 * Error diagnostics - runs quick checks to identify the specific cause
 * of a generic "process exited" error from the SDK.
 */

import Anthropic from '@anthropic-ai/sdk';
import { getLastApiError } from '../network-interceptor.ts';
import { getAnthropicApiKey, getAnthropicBaseUrl, getClaudeOAuthToken, type AuthType } from '../config/storage.ts';

export type DiagnosticCode =
  | 'billing_error'         // HTTP 402 from Anthropic API
  | 'token_expired'
  | 'invalid_credentials'
  | 'rate_limited'          // HTTP 429 from Anthropic API
  | 'mcp_unreachable'
  | 'service_unavailable'
  | 'unknown_error';

export interface DiagnosticResult {
  code: DiagnosticCode;
  title: string;
  message: string;
  /** Diagnostic check results for debugging */
  details: string[];
}

interface DiagnosticConfig {
  authType?: AuthType;
  workspaceId?: string;
  rawError: string;
}

interface CheckResult {
  ok: boolean;
  detail: string;
  failCode?: DiagnosticCode;
  failTitle?: string;
  failMessage?: string;
}

/** Run a check with a timeout, returns default result if times out */
async function withTimeout<T>(promise: Promise<T>, timeoutMs: number, defaultValue: T): Promise<T> {
  const timeoutPromise = new Promise<T>((resolve) => setTimeout(() => resolve(defaultValue), timeoutMs));
  return Promise.race([promise, timeoutPromise]);
}

/**
 * Check if a recent API error was captured during the failed request.
 * This is the most accurate source of truth for API failures since it
 * captures the actual HTTP status code before the SDK wraps it.
 */
async function checkCapturedApiError(): Promise<CheckResult> {
  const apiError = getLastApiError();

  if (!apiError) {
    return { ok: true, detail: '✓ API error: None captured' };
  }

  // HTTP 402 - Payment Required
  if (apiError.status === 402) {
    return {
      ok: false,
      detail: `✗ API 错误: 402 ${apiError.message}`,
      failCode: 'billing_error',
      failTitle: '需要付款',
      failMessage: apiError.message || '您的 Anthropic API 账户存在计费问题。',
    };
  }

  // HTTP 401 - Unauthorized / Invalid Credentials
  if (apiError.status === 401) {
    return {
      ok: false,
      detail: `✗ API 错误: 401 ${apiError.message}`,
      failCode: 'invalid_credentials',
      failTitle: '凭证无效',
      failMessage: apiError.message || '您的 API 凭证无效或已过期。',
    };
  }

  // HTTP 429 - Rate Limited
  if (apiError.status === 429) {
    return {
      ok: false,
      detail: `✗ API 错误: 429 ${apiError.message}`,
      failCode: 'rate_limited',
      failTitle: '请求频率限制',
      failMessage: '请求过于频繁，请稍后再试。',
    };
  }

  // HTTP 5xx - Service Error
  if (apiError.status >= 500) {
    return {
      ok: false,
      detail: `✗ API 错误: ${apiError.status} ${apiError.message}`,
      failCode: 'service_unavailable',
      failTitle: 'Anthropic 服务错误',
      failMessage: `Anthropic API 返回错误 (${apiError.status})，这通常是临时问题。`,
    };
  }

  // Other 4xx errors - report but don't fail (might be expected)
  // Include the message so users can see what actually went wrong
  return { ok: true, detail: `✓ API error: ${apiError.status} - ${apiError.message}` };
}

/**
 * Derive a user-facing label from the configured API base URL.
 * Used in diagnostics messages so errors reference the correct provider.
 */
function getProviderLabel(baseUrl: string): string {
  if (baseUrl.includes('openrouter')) return 'OpenRouter';
  if (baseUrl.includes('anthropic')) return 'Anthropic';
  return 'API endpoint';
}

/**
 * Check if the configured API endpoint is reachable.
 * Uses a simple HEAD request to check connectivity without authentication.
 * Respects ANTHROPIC_BASE_URL so diagnostics are meaningful for all providers.
 */
async function checkApiAvailability(): Promise<CheckResult> {
  // Use the same base URL resolution as network-interceptor.ts
  const baseUrl = process.env.ANTHROPIC_BASE_URL?.trim() || 'https://api.anthropic.com';
  const label = getProviderLabel(baseUrl);

  try {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 3000);

    try {
      // HEAD request doesn't require auth and checks if service is up
      const response = await fetch(`${baseUrl}/v1/models`, {
        method: 'HEAD',
        signal: controller.signal,
      });
      clearTimeout(timeoutId);

      // Any response means the service is reachable
      // 401/403 = reachable but auth required (expected without key)
      // 5xx = service issues
      if (response.status >= 500) {
        return {
          ok: false,
          detail: `✗ ${label}: 服务错误 (${response.status})`,
          failCode: 'service_unavailable',
          failTitle: `${label} 服务错误`,
          failMessage: `${label} 出现问题，请稍后再试。`,
        };
      }

      return { ok: true, detail: `✓ ${label}: 可达 (${response.status})` };
    } catch (fetchError) {
      clearTimeout(timeoutId);

      if (fetchError instanceof Error && fetchError.name === 'AbortError') {
        return {
          ok: false,
          detail: `✗ ${label}: 超时`,
          failCode: 'service_unavailable',
          failTitle: `${label} 不可达`,
          failMessage: `无法连接到 ${label}，请检查网络连接。`,
        };
      }

      const msg = fetchError instanceof Error ? fetchError.message : String(fetchError);
      if (msg.includes('ECONNREFUSED') || msg.includes('ENOTFOUND') || msg.includes('fetch failed')) {
        return {
          ok: false,
          detail: `✗ ${label}: 不可达 (${msg})`,
          failCode: 'service_unavailable',
          failTitle: `${label} 不可达`,
          failMessage: `无法连接到 ${label}，请检查网络连接。`,
        };
      }

      return { ok: true, detail: `✓ ${label}: 未知 (${msg})` };
    }
  } catch (error) {
    const msg = error instanceof Error ? error.message : String(error);
    return { ok: true, detail: `✓ ${label}: 检查失败 (${msg})` };
  }
}

/** Check workspace token expiry - placeholder, always returns valid */
async function checkWorkspaceToken(_workspaceId: string): Promise<CheckResult> {
  // Token expiry checking was removed in a refactoring
  // For now, just assume tokens are valid - the actual API call will fail if expired
  return { ok: true, detail: '✓ Workspace token: Present' };
}

/**
 * Validate an API key by making a test request to Anthropic.
 * Uses models.list() which is lightweight and doesn't incur AI costs.
 */
async function validateApiKeyWithAnthropic(apiKey: string, baseUrl?: string | null): Promise<CheckResult> {
  try {
    const client = new Anthropic({
      apiKey,
      ...(baseUrl ? { baseURL: baseUrl } : {})
    });
    const result = await client.models.list();
    const modelCount = result.data?.length ?? 0;
    return {
      ok: true,
      detail: `✓ API key: Valid (${modelCount} models available)`,
    };
  } catch (error) {
    const msg = error instanceof Error ? error.message : String(error);

    // 401 = Invalid key
    if (msg.includes('401') || msg.includes('invalid') || msg.includes('Unauthorized') || msg.includes('authentication')) {
      return {
        ok: false,
        detail: '✗ API 密钥: 无效或已过期',
        failCode: 'invalid_credentials',
        failTitle: 'API 密钥无效',
        failMessage: '您的 Anthropic API 密钥无效或已过期，请在设置中更新。',
      };
    }

    // 403 = Key valid but no permission
    if (msg.includes('403') || msg.includes('permission') || msg.includes('Forbidden')) {
      return {
        ok: false,
        detail: '✗ API 密钥: 权限不足',
        failCode: 'invalid_credentials',
        failTitle: 'API 密钥权限错误',
        failMessage: '您的 API 密钥没有访问 API 的权限，请检查 Anthropic 控制台。',
      };
    }

    // Network/other errors - don't fail on these, just note them
    return {
      ok: true,
      detail: `✓ API 密钥: 跳过验证 (${msg.slice(0, 50)})`,
    };
  }
}

/** Check API key presence and validity */
async function checkApiKey(): Promise<CheckResult> {
  try {
    const apiKey = await getAnthropicApiKey();
    const baseUrl = getAnthropicBaseUrl();

    if (!apiKey) {
      return {
        ok: false,
        detail: '✗ API 密钥: 未找到',
        failCode: 'invalid_credentials',
        failTitle: 'API 密钥缺失',
        failMessage: '您的 Anthropic API 密钥缺失，请在设置中添加。',
      };
    }

    // Actually validate the key works
    return await validateApiKeyWithAnthropic(apiKey, baseUrl);
  } catch (error) {
    const msg = error instanceof Error ? error.message : String(error);
    return { ok: true, detail: `✓ API 密钥: 检查失败 (${msg})` };
  }
}

/** Check OAuth token presence */
async function checkOAuthToken(): Promise<CheckResult> {
  try {
    const token = await getClaudeOAuthToken();
    if (!token) {
      return {
        ok: false,
        detail: '✗ OAuth 令牌: 未找到',
        failCode: 'invalid_credentials',
        failTitle: 'OAuth 令牌缺失',
        failMessage: '您的 Claude Max OAuth 令牌缺失，请重新认证。',
      };
    }
    return { ok: true, detail: '✓ OAuth 令牌: 已存在' };
  } catch (error) {
    const msg = error instanceof Error ? error.message : String(error);
    return { ok: true, detail: `✓ OAuth 令牌: 检查失败 (${msg})` };
  }
}

/** Check MCP server connectivity with a quick HEAD request */
async function checkMcpConnectivity(mcpUrl: string): Promise<CheckResult> {
  try {
    // Parse the URL to get just the base server
    const url = new URL(mcpUrl);
    const baseUrl = `${url.protocol}//${url.host}`;

    // Quick HEAD request with short timeout
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 3000);

    try {
      const response = await fetch(baseUrl, {
        method: 'HEAD',
        signal: controller.signal,
      });
      clearTimeout(timeoutId);

      // Any response (even 4xx) means the server is reachable
      return { ok: true, detail: `✓ MCP 服务器: 可达 (${response.status})` };
    } catch (fetchError) {
      clearTimeout(timeoutId);
      if (fetchError instanceof Error && fetchError.name === 'AbortError') {
        return {
          ok: false,
          detail: '✗ MCP 服务器: 超时',
          failCode: 'mcp_unreachable',
          failTitle: 'MCP 服务器不可达',
          failMessage: '无法连接到 Craft MCP 服务器（超时），请检查网络连接。',
        };
      }
      const msg = fetchError instanceof Error ? fetchError.message : String(fetchError);
      // Check for common network errors
      if (msg.includes('ECONNREFUSED') || msg.includes('ENOTFOUND') || msg.includes('fetch failed')) {
        return {
          ok: false,
          detail: `✗ MCP 服务器: 不可达 (${msg})`,
          failCode: 'mcp_unreachable',
          failTitle: 'MCP 服务器不可达',
          failMessage: '无法连接到 Craft MCP 服务器，请检查网络连接。',
        };
      }
      return { ok: true, detail: `✓ MCP 服务器: 未知 (${msg})` };
    }
  } catch (error) {
    const msg = error instanceof Error ? error.message : String(error);
    return { ok: true, detail: `✓ MCP 服务器: 检查失败 (${msg})` };
  }
}

/**
 * Run error diagnostics to identify the specific cause of a failure.
 * All checks run in parallel with 5s timeouts.
 */
export async function runErrorDiagnostics(config: DiagnosticConfig): Promise<DiagnosticResult> {
  const { authType, workspaceId, rawError } = config;
  const details: string[] = [];
  const defaultResult: CheckResult = { ok: true, detail: '? Check: Timeout' };

  // Build list of checks to run based on config
  const checks: Promise<CheckResult>[] = [];

  // 0. FIRST: Check captured API error (most accurate source of truth)
  // This captures the actual HTTP status code from the failed request
  checks.push(withTimeout(checkCapturedApiError(), 1000, defaultResult));

  // 1. API endpoint availability check (uses configured base URL)
  checks.push(withTimeout(checkApiAvailability(), 4000, defaultResult));

  // 2. API key check with validation (only for api_key auth)
  if (authType === 'api_key') {
    checks.push(withTimeout(checkApiKey(), 5000, defaultResult));
  }

  // 3. OAuth token check (only for oauth_token auth)
  if (authType === 'oauth_token') {
    checks.push(withTimeout(checkOAuthToken(), 5000, defaultResult));
  }

  // Run all checks in parallel
  const results = await Promise.all(checks);

  // Collect details and find first failure
  let firstFailure: CheckResult | null = null;
  for (const result of results) {
    details.push(result.detail);
    if (!result.ok && !firstFailure) {
      firstFailure = result;
    }
  }

  // Add raw error to details
  details.push(`Raw error: ${rawError.slice(0, 200)}${rawError.length > 200 ? '...' : ''}`);

  // Return specific issue if found
  if (firstFailure && firstFailure.failCode && firstFailure.failTitle && firstFailure.failMessage) {
    return {
      code: firstFailure.failCode,
      title: firstFailure.failTitle,
      message: firstFailure.failMessage,
      details,
    };
  }

  // All checks passed but still failed - likely Anthropic service issue
  return {
    code: 'service_unavailable',
    title: '服务不可用',
    message: 'AI 服务出现问题，所有凭证均有效。请稍后再试。',
    details,
  };
}
