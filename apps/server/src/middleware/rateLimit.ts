import type { Context, Next } from 'hono'
import type { Env } from '../app'

// Simple in-memory rate limiter
// In production, use Redis for distributed rate limiting
interface RateLimitEntry {
  count: number
  resetAt: number
}

const rateLimitStore = new Map<string, RateLimitEntry>()

// Clean up expired entries periodically
setInterval(() => {
  const now = Date.now()
  for (const [key, entry] of rateLimitStore.entries()) {
    if (entry.resetAt < now) {
      rateLimitStore.delete(key)
    }
  }
}, 60000)  // Clean up every minute

interface RateLimitOptions {
  windowMs: number    // Time window in milliseconds
  max: number         // Max requests per window
  keyGenerator?: (c: Context<Env>) => string
}

/**
 * Create a rate limiting middleware
 */
export function createRateLimiter(options: RateLimitOptions) {
  const {
    windowMs,
    max,
    keyGenerator = (c) => c.get('userId') || c.req.header('x-forwarded-for') || 'anonymous',
  } = options
  
  return async (c: Context<Env>, next: Next) => {
    const key = `rate:${keyGenerator(c)}`
    const now = Date.now()
    
    let entry = rateLimitStore.get(key)
    
    if (!entry || entry.resetAt < now) {
      entry = {
        count: 0,
        resetAt: now + windowMs,
      }
    }
    
    entry.count++
    rateLimitStore.set(key, entry)
    
    // Set rate limit headers
    c.header('X-RateLimit-Limit', String(max))
    c.header('X-RateLimit-Remaining', String(Math.max(0, max - entry.count)))
    c.header('X-RateLimit-Reset', String(Math.ceil(entry.resetAt / 1000)))
    
    if (entry.count > max) {
      c.header('Retry-After', String(Math.ceil((entry.resetAt - now) / 1000)))
      return c.json({
        error: 'Too many requests',
        code: 'RATE_LIMITED',
        retryAfter: Math.ceil((entry.resetAt - now) / 1000),
      }, 429)
    }
    
    await next()
  }
}

// Pre-configured rate limiters
export const authRateLimit = createRateLimiter({
  windowMs: 15 * 60 * 1000,  // 15 minutes
  max: 10,                    // 10 attempts per window
  keyGenerator: (c) => `auth:${c.req.header('x-forwarded-for') || 'unknown'}`,
})

export const apiRateLimit = createRateLimiter({
  windowMs: 60 * 1000,       // 1 minute
  max: 60,                    // 60 requests per minute
})

export const agentRateLimit = createRateLimiter({
  windowMs: 60 * 1000,       // 1 minute
  max: 10,                    // 10 agent requests per minute
})
