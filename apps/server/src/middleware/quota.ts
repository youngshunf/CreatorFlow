import type { Context, Next } from 'hono'
import { SubscriptionService } from '../services/subscription'
import { PLANS, isWithinLimit } from '../config/plans'
import type { Env } from '../app'

/**
 * Quota middleware - checks if user has remaining quota before executing agent
 */
export const quotaMiddleware = async (c: Context<Env>, next: Next) => {
  const userId = c.get('userId')
  
  if (!userId) {
    return c.json({ error: 'Unauthorized' }, 401)
  }
  
  // Get user's subscription and usage
  const [subscription, usage] = await Promise.all([
    SubscriptionService.getByUserId(userId),
    SubscriptionService.getUsage(userId),
  ])
  
  const plan = subscription?.plan || 'free'
  const limits = PLANS[plan as keyof typeof PLANS]?.limits || PLANS.free!.limits
  
  // Check daily message limit
  if (!isWithinLimit(usage.messagesToday, limits.messagesPerDay)) {
    return c.json({
      error: 'Daily message limit exceeded',
      code: 'QUOTA_EXCEEDED',
      details: {
        limit: limits.messagesPerDay,
        used: usage.messagesToday,
        resetAt: getEndOfDay(),
      },
    }, 429)
  }
  
  // Check monthly token limit
  if (!isWithinLimit(usage.tokensThisMonth, limits.tokensPerMonth)) {
    return c.json({
      error: 'Monthly token limit exceeded',
      code: 'QUOTA_EXCEEDED',
      details: {
        limit: limits.tokensPerMonth,
        used: usage.tokensThisMonth,
        resetAt: getEndOfMonth(),
      },
    }, 429)
  }
  
  await next()
}

/**
 * Get end of current day (UTC)
 */
function getEndOfDay(): string {
  const now = new Date()
  const endOfDay = new Date(Date.UTC(
    now.getUTCFullYear(),
    now.getUTCMonth(),
    now.getUTCDate() + 1,
    0, 0, 0, 0
  ))
  return endOfDay.toISOString()
}

/**
 * Get end of current month (UTC)
 */
function getEndOfMonth(): string {
  const now = new Date()
  const endOfMonth = new Date(Date.UTC(
    now.getUTCFullYear(),
    now.getUTCMonth() + 1,
    1, 0, 0, 0, 0
  ))
  return endOfMonth.toISOString()
}
