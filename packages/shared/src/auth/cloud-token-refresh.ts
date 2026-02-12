/**
 * 云端令牌刷新服务
 *
 * 处理 access_token 的自动刷新逻辑。
 * 云端 API 要求 refresh_token 通过 Cookie 传递。
 */

import { debug } from '../utils/debug.ts';

/**
 * 刷新云端 access_token
 *
 * @param apiBaseUrl - 云端 API 基础 URL（如 https://api.example.com/api/v1）
 * @param refreshToken - refresh_token
 * @returns 新的 accessToken 和过期时间
 */
export async function refreshCloudToken(
  apiBaseUrl: string,
  refreshToken: string
): Promise<{ accessToken: string; expiresAt: number }> {
  // 从 gatewayUrl 推导 API base URL
  // gatewayUrl 格式: https://xxx/api/v1/llm/proxy → 需要 https://xxx/api/v1/auth/refresh
  const baseUrl = apiBaseUrl.replace(/\/llm\/proxy\/?$/, '')

  const url = `${baseUrl}/auth/refresh`
  debug(`[CloudTokenRefresh] Refreshing token at ${url}`)

  const response = await fetch(url, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Cookie': `fba_refresh_token=${refreshToken}`,
    },
  })

  if (!response.ok) {
    const text = await response.text().catch(() => '')
    throw new Error(`Cloud token refresh failed (${response.status}): ${text}`)
  }

  const data = await response.json() as {
    access_token?: string;
    access_token_expire_time?: string;
    data?: {
      access_token?: string;
      access_token_expire_time?: string;
    };
  }

  // 兼容两种响应格式：直接字段或嵌套在 data 中
  const result = data.data || data
  const accessToken = result.access_token
  const expireTime = result.access_token_expire_time

  if (!accessToken) {
    throw new Error('Cloud token refresh response missing access_token')
  }

  const expiresAt = expireTime
    ? new Date(expireTime).getTime()
    : Date.now() + 24 * 60 * 60 * 1000 // 默认 24 小时

  debug(`[CloudTokenRefresh] Token refreshed, expires at ${new Date(expiresAt).toISOString()}`)

  return { accessToken, expiresAt }
}

/**
 * 检查 access_token 是否即将过期（默认 5 分钟缓冲）
 */
export function isTokenExpiringSoon(expiresAt: number, bufferMs = 5 * 60 * 1000): boolean {
  return Date.now() > expiresAt - bufferMs
}

/**
 * 检查 refresh_token 是否已过期
 */
export function isRefreshTokenExpired(refreshExpiresAt: number): boolean {
  return Date.now() > refreshExpiresAt
}
