import { verify } from 'hono/jwt'
import type { Context, Next } from 'hono'
import type { Env } from '../app'

export interface JWTPayload {
  userId: string
  email: string
  iat?: number
  exp?: number
}

/**
 * Verify JWT token and return payload
 */
export async function verifyToken(token: string): Promise<JWTPayload | null> {
  try {
    const payload = await verify(token, process.env.JWT_SECRET!, 'HS256')
    return payload as JWTPayload
  } catch (err) {
    console.error('[JWT] Verification failed:', err)
    return null
  }
}

/**
 * Auth middleware - verifies JWT token from Authorization header
 */
export const authMiddleware = async (c: Context<Env>, next: Next) => {
  const authHeader = c.req.header('Authorization')
  
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return c.json({ error: 'Missing or invalid authorization header' }, 401)
  }
  
  const token = authHeader.slice(7)  // Remove 'Bearer '
  
  const payload = await verifyToken(token)
  
  if (!payload) {
    return c.json({ error: 'Invalid or expired token' }, 401)
  }
  
  // Set user info in context
  c.set('userId', payload.userId)
  c.set('userEmail', payload.email)
  
  await next()
}

/**
 * Optional auth middleware - sets user info if token present, but doesn't require it
 */
export const optionalAuthMiddleware = async (c: Context<Env>, next: Next) => {
  const authHeader = c.req.header('Authorization')
  
  if (authHeader && authHeader.startsWith('Bearer ')) {
    const token = authHeader.slice(7)
    const payload = await verifyToken(token)
    
    if (payload) {
      c.set('userId', payload.userId)
      c.set('userEmail', payload.email)
    }
  }
  
  await next()
}

/**
 * WebSocket auth - verify token from message payload
 * Returns userId if valid, null otherwise
 */
export async function wsAuthMiddleware(token: string): Promise<string | null> {
  const payload = await verifyToken(token)
  return payload?.userId ?? null
}

/**
 * API key auth middleware - for programmatic access
 */
export const apiKeyMiddleware = async (c: Context<Env>, next: Next) => {
  const apiKey = c.req.header('X-API-Key')
  
  if (!apiKey) {
    return c.json({ error: 'Missing API key' }, 401)
  }
  
  // Import here to avoid circular dependency
  const { UserService } = await import('../services/user')
  const user = await UserService.getByApiKey(apiKey)
  
  if (!user) {
    return c.json({ error: 'Invalid API key' }, 401)
  }
  
  c.set('userId', user.id)
  c.set('userEmail', user.email)
  
  await next()
}
