import { Hono } from 'hono'
import { zValidator } from '@hono/zod-validator'
import { z } from 'zod'
import { eq, desc, and, gte } from 'drizzle-orm'
import { authMiddleware } from '../middleware/auth'
import { UserService } from '../services/user'
import { SubscriptionService } from '../services/subscription'
import { db } from '../db'
import { orders } from '../db/schema'
import type { Env } from '../app'

const admin = new Hono<Env>()

// Admin middleware - check if user is admin
const adminMiddleware = async (c: any, next: () => Promise<void>) => {
  const userId = c.get('userId')
  const user = await UserService.getById(userId)
  
  // Check if user is admin (you can implement your own logic)
  // For now, check against an env variable list of admin emails
  const adminEmails = (process.env.ADMIN_EMAILS || '').split(',').map(e => e.trim())
  
  if (!user || !adminEmails.includes(user.email)) {
    return c.json({ error: 'Forbidden' }, 403)
  }
  
  await next()
}

// Dashboard stats
admin.get('/stats', authMiddleware, adminMiddleware, async (c) => {
  const now = new Date()
  const startOfDay = new Date(Date.UTC(now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate()))
  
  const [userCount, subsByPlan, revenue, recentOrders] = await Promise.all([
    UserService.count(),
    SubscriptionService.countByPlan(),
    SubscriptionService.getRevenueStats(),
    db ? db.select().from(orders).orderBy(desc(orders.createdAt)).limit(10) : [],
  ])
  
  // 计算活跃订阅数
  const activeSubscriptions = Object.values(subsByPlan).reduce((a, b) => a + b, 0)
  
  return c.json({
    userCount,
    activeSubscriptions,
    revenue,
    recentOrders,
  })
})

// List users
admin.get(
  '/users',
  authMiddleware,
  adminMiddleware,
  zValidator('query', z.object({
    page: z.coerce.number().default(1),
    limit: z.coerce.number().default(20),
    search: z.string().optional(),
  })),
  async (c) => {
    const { page, limit, search } = c.req.valid('query')
    
    const result = await UserService.list({ page, limit, search })
    
    return c.json(result)
  }
)

// Get user details
admin.get('/users/:userId', authMiddleware, adminMiddleware, async (c) => {
  const userId = c.req.param('userId')
  
  const [user, subscription, usage] = await Promise.all([
    UserService.getById(userId),
    SubscriptionService.getByUserId(userId),
    SubscriptionService.getUsage(userId),
  ])
  
  if (!user) {
    return c.json({ error: 'User not found' }, 404)
  }
  
  return c.json({
    user,
    subscription,
    usage,
  })
})

// Update user
admin.patch(
  '/users/:userId',
  authMiddleware,
  adminMiddleware,
  zValidator('json', z.object({
    name: z.string().optional(),
    emailVerified: z.boolean().optional(),
  })),
  async (c) => {
    const userId = c.req.param('userId')
    const updates = c.req.valid('json')
    
    const user = await UserService.update(userId, updates)
    
    return c.json({ user })
  }
)

// List subscriptions
admin.get(
  '/subscriptions',
  authMiddleware,
  adminMiddleware,
  zValidator('query', z.object({
    page: z.coerce.number().default(1),
    limit: z.coerce.number().default(20),
    status: z.string().optional(),
    plan: z.string().optional(),
  })),
  async (c) => {
    const { page, limit, status, plan } = c.req.valid('query')
    
    const result = await SubscriptionService.list({ page, limit, status, plan })
    
    return c.json(result)
  }
)

// Manual subscription management (for support cases)
admin.post(
  '/subscriptions/:userId/grant',
  authMiddleware,
  adminMiddleware,
  zValidator('json', z.object({
    plan: z.enum(['pro', 'team']),
    durationDays: z.number().min(1).max(365),
    reason: z.string(),
  })),
  async (c) => {
    const userId = c.req.param('userId')
    const { plan, durationDays, reason } = c.req.valid('json')
    
    const subscription = await SubscriptionService.grantManual({
      userId,
      plan,
      durationDays,
      reason,
      grantedBy: c.get('userId'),
    })
    
    return c.json({ subscription })
  }
)

// Revoke subscription
admin.post(
  '/subscriptions/:userId/revoke',
  authMiddleware,
  adminMiddleware,
  zValidator('json', z.object({
    reason: z.string(),
  })),
  async (c) => {
    const userId = c.req.param('userId')
    const { reason } = c.req.valid('json')
    
    await SubscriptionService.revoke({
      userId,
      reason,
      revokedBy: c.get('userId'),
    })
    
    return c.json({ success: true })
  }
)

// List orders
admin.get(
  '/orders',
  authMiddleware,
  adminMiddleware,
  zValidator('query', z.object({
    page: z.coerce.number().default(1),
    limit: z.coerce.number().default(20),
    status: z.string().optional(),
    channel: z.string().optional(),
  })),
  async (c) => {
    const { page, limit, status, channel } = c.req.valid('query')
    const offset = (page - 1) * limit
    
    if (!db) {
      return c.json({ orders: [], total: 0, page, limit })
    }
    
    const conditions = []
    if (status) conditions.push(eq(orders.status, status))
    if (channel) conditions.push(eq(orders.channel, channel))
    
    const whereClause = conditions.length > 0 ? and(...conditions) : undefined
    
    const [results, countResult] = await Promise.all([
      db.select().from(orders).where(whereClause).orderBy(desc(orders.createdAt)).limit(limit).offset(offset),
      db.select().from(orders).where(whereClause),
    ])
    
    return c.json({
      orders: results,
      total: countResult.length,
      page,
      limit,
    })
  }
)

// Grant subscription (alternative endpoint for admin frontend)
admin.post(
  '/subscriptions/grant',
  authMiddleware,
  adminMiddleware,
  zValidator('json', z.object({
    userId: z.string(),
    plan: z.enum(['pro', 'team']),
    months: z.number().min(1).max(24),
    reason: z.string(),
  })),
  async (c) => {
    const { userId, plan, months, reason } = c.req.valid('json')
    
    const subscription = await SubscriptionService.grantManual({
      userId,
      plan,
      durationDays: months * 30,
      reason,
      grantedBy: c.get('userId'),
    })
    
    return c.json({ subscription })
  }
)

export { admin }
