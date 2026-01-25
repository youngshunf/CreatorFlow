import { eq, and, gte, sql, count } from 'drizzle-orm'
import { nanoid } from 'nanoid'
import { db, subscriptions, usageRecords, auditLogs, type Subscription, type NewSubscription } from '../db'

export interface UsageStats {
  messagesToday: number
  tokensThisMonth: number
}

export class SubscriptionService {
  /**
   * Get subscription by user ID
   */
  static async getByUserId(userId: string): Promise<Subscription | null> {
    if (!db) return null
    
    const sub = await db.query.subscriptions.findFirst({
      where: and(
        eq(subscriptions.userId, userId),
        eq(subscriptions.status, 'active')
      ),
    })
    
    return sub ?? null
  }
  
  /**
   * Get Stripe customer ID for user
   */
  static async getStripeCustomerId(userId: string): Promise<string | null> {
    if (!db) return null
    
    const sub = await db.query.subscriptions.findFirst({
      where: eq(subscriptions.userId, userId),
    })
    
    return sub?.stripeCustomerId ?? null
  }
  
  /**
   * Set Stripe customer ID for user
   */
  static async setStripeCustomerId(userId: string, customerId: string): Promise<void> {
    if (!db) return
    
    const existing = await db.query.subscriptions.findFirst({
      where: eq(subscriptions.userId, userId),
    })
    
    if (existing) {
      await db.update(subscriptions)
        .set({ stripeCustomerId: customerId, updatedAt: new Date() })
        .where(eq(subscriptions.id, existing.id))
    } else {
      await db.insert(subscriptions).values({
        id: nanoid(),
        userId,
        plan: 'free',
        status: 'active',
        stripeCustomerId: customerId,
      })
    }
  }
  
  /**
   * Activate subscription after checkout
   */
  static async activate(data: {
    userId: string
    plan: string
    stripeSubscriptionId: string
    stripeCustomerId: string
  }): Promise<Subscription> {
    if (!db) throw new Error('Database not configured')
    
    // Deactivate any existing subscription
    await db.update(subscriptions)
      .set({ status: 'canceled', updatedAt: new Date() })
      .where(and(
        eq(subscriptions.userId, data.userId),
        eq(subscriptions.status, 'active')
      ))
    
    // Create new subscription
    const newSub: NewSubscription = {
      id: nanoid(),
      userId: data.userId,
      plan: data.plan,
      status: 'active',
      stripeCustomerId: data.stripeCustomerId,
      stripeSubscriptionId: data.stripeSubscriptionId,
      currentPeriodStart: new Date(),
      currentPeriodEnd: new Date(Date.now() + 30 * 24 * 60 * 60 * 1000),  // +30 days
    }
    
    const [sub] = await db.insert(subscriptions).values(newSub).returning()
    return sub!
  }
  
  /**
   * Update subscription period
   */
  static async updatePeriod(stripeSubscriptionId: string, start: Date, end: Date): Promise<void> {
    if (!db) return
    
    await db.update(subscriptions)
      .set({
        currentPeriodStart: start,
        currentPeriodEnd: end,
        status: 'active',
        updatedAt: new Date(),
      })
      .where(eq(subscriptions.stripeSubscriptionId, stripeSubscriptionId))
  }
  
  /**
   * Mark subscription as past due
   */
  static async markPastDue(stripeSubscriptionId: string): Promise<void> {
    if (!db) return
    
    await db.update(subscriptions)
      .set({ status: 'past_due', updatedAt: new Date() })
      .where(eq(subscriptions.stripeSubscriptionId, stripeSubscriptionId))
  }
  
  /**
   * Mark subscription as canceling (will cancel at period end)
   */
  static async markCanceling(stripeSubscriptionId: string): Promise<void> {
    if (!db) return
    
    await db.update(subscriptions)
      .set({ 
        status: 'canceling',
        cancelAtPeriodEnd: true,
        updatedAt: new Date(),
      })
      .where(eq(subscriptions.stripeSubscriptionId, stripeSubscriptionId))
  }
  
  /**
   * Cancel subscription
   */
  static async cancel(stripeSubscriptionId: string): Promise<void> {
    if (!db) return
    
    await db.update(subscriptions)
      .set({ status: 'canceled', updatedAt: new Date() })
      .where(eq(subscriptions.stripeSubscriptionId, stripeSubscriptionId))
  }
  
  /**
   * Get user's usage stats
   */
  static async getUsage(userId: string): Promise<UsageStats> {
    if (!db) {
      return { messagesToday: 0, tokensThisMonth: 0 }
    }
    
    const now = new Date()
    const startOfDay = new Date(Date.UTC(now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate()))
    const startOfMonth = new Date(Date.UTC(now.getUTCFullYear(), now.getUTCMonth(), 1))
    
    const [messagesResult, tokensResult] = await Promise.all([
      // Count messages today
      db.select({ count: count() })
        .from(usageRecords)
        .where(and(
          eq(usageRecords.userId, userId),
          eq(usageRecords.type, 'message'),
          gte(usageRecords.createdAt, startOfDay)
        )),
      // Sum tokens this month
      db.select({ sum: sql<number>`COALESCE(SUM(amount), 0)` })
        .from(usageRecords)
        .where(and(
          eq(usageRecords.userId, userId),
          eq(usageRecords.type, 'tokens'),
          gte(usageRecords.createdAt, startOfMonth)
        )),
    ])
    
    return {
      messagesToday: messagesResult[0]?.count ?? 0,
      tokensThisMonth: Number(tokensResult[0]?.sum ?? 0),
    }
  }
  
  /**
   * Check if user has quota available
   */
  static async checkQuota(userId: string): Promise<boolean> {
    // Import here to avoid circular dependency
    const { PLANS, isWithinLimit } = await import('../config/plans')
    
    const [subscription, usage] = await Promise.all([
      this.getByUserId(userId),
      this.getUsage(userId),
    ])
    
    const plan = subscription?.plan || 'free'
    const limits = PLANS[plan as keyof typeof PLANS]?.limits || PLANS.free!.limits
    
    return (
      isWithinLimit(usage.messagesToday, limits.messagesPerDay) &&
      isWithinLimit(usage.tokensThisMonth, limits.tokensPerMonth)
    )
  }
  
  /**
   * Record usage
   */
  static async recordUsage(userId: string, usage: { inputTokens?: number; outputTokens?: number }): Promise<void> {
    if (!db) return
    
    const totalTokens = (usage.inputTokens || 0) + (usage.outputTokens || 0)
    
    // Record message
    await db.insert(usageRecords).values({
      id: nanoid(),
      userId,
      type: 'message',
      amount: 1,
    })
    
    // Record tokens
    if (totalTokens > 0) {
      await db.insert(usageRecords).values({
        id: nanoid(),
        userId,
        type: 'tokens',
        amount: totalTokens,
        metadata: { inputTokens: usage.inputTokens, outputTokens: usage.outputTokens },
      })
    }
  }
  
  /**
   * Count subscriptions by plan (admin)
   */
  static async countByPlan(): Promise<Record<string, number>> {
    if (!db) return {}
    
    const results = await db.select({
      plan: subscriptions.plan,
      count: count(),
    })
      .from(subscriptions)
      .where(eq(subscriptions.status, 'active'))
      .groupBy(subscriptions.plan)
    
    return Object.fromEntries(results.map(r => [r.plan, r.count]))
  }
  
  /**
   * Get revenue stats (admin)
   */
  static async getRevenueStats(): Promise<{ mrr: number }> {
    // This would normally query Stripe or calculate from local data
    // For now, return placeholder
    return { mrr: 0 }
  }
  
  /**
   * List subscriptions (admin)
   */
  static async list(options: { page: number; limit: number; status?: string; plan?: string }) {
    if (!db) return { subscriptions: [], total: 0 }
    
    const { page, limit, status, plan } = options
    const offset = (page - 1) * limit
    
    const conditions = []
    if (status) conditions.push(eq(subscriptions.status, status))
    if (plan) conditions.push(eq(subscriptions.plan, plan))
    
    const whereClause = conditions.length > 0 ? and(...conditions) : undefined
    
    const [results, totalResult] = await Promise.all([
      db.query.subscriptions.findMany({
        where: whereClause,
        limit,
        offset,
        orderBy: (subscriptions, { desc }) => [desc(subscriptions.createdAt)],
      }),
      db.select({ count: count() }).from(subscriptions).where(whereClause),
    ])
    
    return {
      subscriptions: results,
      total: totalResult[0]?.count ?? 0,
      page,
      limit,
    }
  }
  
  /**
   * Grant manual subscription (admin)
   */
  static async grantManual(data: {
    userId: string
    plan: string
    durationDays: number
    reason: string
    grantedBy: string
  }): Promise<Subscription> {
    if (!db) throw new Error('Database not configured')
    
    const endDate = new Date(Date.now() + data.durationDays * 24 * 60 * 60 * 1000)
    
    const newSub: NewSubscription = {
      id: nanoid(),
      userId: data.userId,
      plan: data.plan,
      status: 'active',
      currentPeriodStart: new Date(),
      currentPeriodEnd: endDate,
      grantedBy: data.grantedBy,
      grantReason: data.reason,
    }
    
    const [sub] = await db.insert(subscriptions).values(newSub).returning()
    
    // Audit log
    await db.insert(auditLogs).values({
      id: nanoid(),
      actorId: data.grantedBy,
      action: 'subscription.grant',
      targetId: data.userId,
      targetType: 'user',
      details: { plan: data.plan, durationDays: data.durationDays, reason: data.reason },
    })
    
    return sub!
  }
  
  /**
   * Revoke subscription (admin)
   */
  static async revoke(data: {
    userId: string
    reason: string
    revokedBy: string
  }): Promise<void> {
    if (!db) return
    
    await db.update(subscriptions)
      .set({ status: 'canceled', updatedAt: new Date() })
      .where(and(
        eq(subscriptions.userId, data.userId),
        eq(subscriptions.status, 'active')
      ))
    
    // Audit log
    await db.insert(auditLogs).values({
      id: nanoid(),
      actorId: data.revokedBy,
      action: 'subscription.revoke',
      targetId: data.userId,
      targetType: 'user',
      details: { reason: data.reason },
    })
  }
}
