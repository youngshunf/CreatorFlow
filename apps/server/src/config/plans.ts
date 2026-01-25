export interface PlanLimits {
  messagesPerDay: number    // -1 for unlimited
  tokensPerMonth: number    // -1 for unlimited
  maxSessions: number       // -1 for unlimited
}

export interface Plan {
  name: string
  price: number             // Monthly price in USD
  priceCNY: number          // Monthly price in CNY (分)
  yearlyPriceCNY: number    // Yearly price in CNY (分)
  limits: PlanLimits
}

export const PLANS: Record<string, Plan> = {
  free: {
    name: 'Free',
    price: 0,
    priceCNY: 0,
    yearlyPriceCNY: 0,
    limits: {
      messagesPerDay: 20,
      tokensPerMonth: 100_000,
      maxSessions: 10,
    },
  },
  pro: {
    name: 'Pro',
    price: 19,
    priceCNY: 12900,        // 129元/月
    yearlyPriceCNY: 119900, // 1199元/年 (约85折)
    limits: {
      messagesPerDay: 500,
      tokensPerMonth: 2_000_000,
      maxSessions: 100,
    },
  },
  team: {
    name: 'Team',
    price: 49,
    priceCNY: 34900,        // 349元/月
    yearlyPriceCNY: 299900, // 2999元/年 (约72折)
    limits: {
      messagesPerDay: -1,  // unlimited
      tokensPerMonth: 10_000_000,
      maxSessions: -1,     // unlimited
    },
  },
}

export function getPlanPrice(planId: string, period: 'monthly' | 'yearly'): number {
  const plan = getPlan(planId)
  return period === 'monthly' ? plan.priceCNY : plan.yearlyPriceCNY
}

export function getPlan(planId: string): Plan {
  return PLANS[planId] ?? PLANS.free!
}

export function getPlanLimits(planId: string): PlanLimits {
  return getPlan(planId).limits
}

export function isWithinLimit(current: number, limit: number): boolean {
  return limit === -1 || current < limit
}
