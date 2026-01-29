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
  monthly_credits: number;
  current_credits: number;
  used_credits: number;
  purchased_credits: number;
  monthly_remaining?: number;
  bonus_remaining?: number;
  billing_cycle_start: string;
  billing_cycle_end: string;
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
   * 升级订阅（模拟支付）
   */
  upgrade: (tierName: string) => {
    return request.post<PaymentResult>('/user_tier/my/subscription/upgrade', {
      tier_name: tierName,
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
};
