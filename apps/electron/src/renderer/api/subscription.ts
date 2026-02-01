/**
 * 订阅与积分 API
 * @author Ysf
 */

import request from './request';

// ============================================
// Types
// ============================================

export interface CreditBalance {
  id: number;
  credit_type: string;
  original_amount: number;
  used_amount: number;
  remaining_amount: number;
  expires_at: string | null;
  granted_at: string;
  source_type: string;
  description: string | null;
}

export interface SubscriptionInfo {
  user_id: number;
  tier: string;
  tier_display_name: string | null;
  subscription_type: 'monthly' | 'yearly';  // 订阅类型
  monthly_credits: number;
  current_credits: number;
  used_credits: number;
  purchased_credits: number;
  monthly_remaining?: number;
  bonus_remaining?: number;
  billing_cycle_start: string;
  billing_cycle_end: string;
  subscription_start_date?: string | null;  // 订阅开始时间
  subscription_end_date?: string | null;    // 订阅结束时间
  next_grant_date?: string | null;          // 下次赠送时间（年度订阅）
  status: string;
  auto_renew?: boolean;
  balances?: CreditBalance[];
}

export interface SubscriptionTier {
  id: number;
  tier_name: string;
  display_name: string;
  monthly_credits: number;
  monthly_price: number;
  yearly_price?: number | null;     // 年费价格
  yearly_discount?: number | null;  // 年费折扣
  features: Record<string, unknown> | null;
}

export interface CreditPackage {
  id: number;
  package_name: string;
  credits: number;
  price: number;
  bonus_credits: number;
  description: string | null;
}

export interface PaymentResult {
  success: boolean;
  order_id: string;
  message: string;
  new_credits?: number;
  new_tier?: string;
}

export interface UpgradePriceResult {
  can_upgrade: boolean;
  message: string;
  target_tier: string;
  target_tier_display: string;
  subscription_type: string;
  original_price: number;
  remaining_value: number;
  final_price: number;
  remaining_days: number;
  current_tier: string;
  current_subscription_type: string;
}

// ============================================
// API
// ============================================

export const subscriptionApi = {
  /**
   * 获取当前用户订阅信息
   */
  getInfo: () => {
    return request.get<SubscriptionInfo>('/user_tier/my/subscription/info');
  },

  /**
   * 获取订阅等级列表
   */
  getTiers: () => {
    return request.get<SubscriptionTier[]>('/user_tier/my/subscription/tiers');
  },

  /**
   * 获取积分包列表
   */
  getPackages: () => {
    return request.get<CreditPackage[]>('/user_tier/my/subscription/packages');
  },

  /**
   * 计算升级价格
   * @param tierName 订阅等级
   * @param subscriptionType 订阅类型 (monthly/yearly)
   */
  calculateUpgradePrice: (tierName: string, subscriptionType: 'monthly' | 'yearly' = 'monthly') => {
    return request.post<UpgradePriceResult>('/user_tier/my/subscription/upgrade/calculate', {
      tier_name: tierName,
      subscription_type: subscriptionType,
    });
  },

  /**
   * 升级订阅（模拟支付）
   * @param tierName 订阅等级
   * @param subscriptionType 订阅类型 (monthly/yearly)
   */
  upgrade: (tierName: string, subscriptionType: 'monthly' | 'yearly' = 'monthly') => {
    return request.post<PaymentResult>('/user_tier/my/subscription/upgrade', {
      tier_name: tierName,
      subscription_type: subscriptionType,
    });
  },

  /**
   * 购买积分包（模拟支付）
   */
  purchase: (packageId: number) => {
    return request.post<PaymentResult>('/user_tier/my/subscription/purchase', {
      package_id: packageId,
    });
  },

  /**
   * 获取历史积分记录（已过期）
   */
  getBalanceHistory: () => {
    return request.get<CreditBalance[]>('/user_tier/my/subscription/balances/history');
  },
};
