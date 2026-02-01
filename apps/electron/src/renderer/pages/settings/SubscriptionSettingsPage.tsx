/**
 * SubscriptionSettingsPage
 *
 * 订阅与积分设置页面
 * - 显示当前订阅状态和积分余额
 * - 升级订阅
 * - 购买积分包
 *
 * @author Ysf
 */

import { useState, useEffect, useCallback } from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { useT } from '@/context/LocaleContext'
import { routes } from '@/lib/navigate'
import { Crown, Coins, TrendingUp, Calendar, Zap, RefreshCw, Check, X, Loader2, Gift, Sparkles, History, ChevronRight } from 'lucide-react'
import type { DetailsPageMeta } from '@/lib/navigation-registry'
import { cn } from '@/lib/utils'

import {
  SettingsSection,
  SettingsCard,
  SettingsRow,
} from '@/components/settings'

import {
  subscriptionApi,
  type SubscriptionInfo,
  type SubscriptionTier,
  type CreditPackage,
  type UpgradePriceResult,
  type CreditBalance,
} from '@/api/subscription'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'subscription',
}

// ============================================
// Constants
// ============================================

const TIER_COLORS: Record<string, { bg: string; text: string; accent: string; cardBg: string }> = {
  free: { 
    bg: 'bg-gradient-to-br from-slate-500 to-slate-700', 
    text: 'text-white', 
    accent: 'text-slate-200',
    cardBg: 'bg-white/20 backdrop-blur'
  },
  pro: { 
    bg: 'bg-gradient-to-br from-blue-500 to-indigo-600', 
    text: 'text-white', 
    accent: 'text-blue-100',
    cardBg: 'bg-white/20 backdrop-blur'
  },
  max: { 
    bg: 'bg-gradient-to-br from-purple-500 to-pink-600', 
    text: 'text-white', 
    accent: 'text-purple-100',
    cardBg: 'bg-white/20 backdrop-blur'
  },
  ultra: { 
    bg: 'bg-gradient-to-br from-amber-500 to-orange-600', 
    text: 'text-white', 
    accent: 'text-amber-100',
    cardBg: 'bg-white/20 backdrop-blur'
  },
}

// ============================================
// Main Component
// ============================================

export default function SubscriptionSettingsPage() {
  const t = useT()

  // State
  const [subscription, setSubscription] = useState<SubscriptionInfo | null>(null)
  const [tiers, setTiers] = useState<SubscriptionTier[]>([])
  const [packages, setPackages] = useState<CreditPackage[]>([])
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  // Dialog states
  const [showUpgradeDialog, setShowUpgradeDialog] = useState(false)
  const [showCreditDialog, setShowCreditDialog] = useState(false)
  const [showHistoryDialog, setShowHistoryDialog] = useState(false)
  const [historyBalances, setHistoryBalances] = useState<CreditBalance[]>([])
  const [isLoadingHistory, setIsLoadingHistory] = useState(false)
  const [selectedTier, setSelectedTier] = useState<string | null>(null)
  const [selectedSubscriptionType, setSelectedSubscriptionType] = useState<'monthly' | 'yearly'>('monthly')
  const [selectedPackage, setSelectedPackage] = useState<number | null>(null)
  const [isProcessing, setIsProcessing] = useState(false)
  const [result, setResult] = useState<{ success: boolean; message: string } | null>(null)
  const [priceInfo, setPriceInfo] = useState<UpgradePriceResult | null>(null)
  const [isCalculating, setIsCalculating] = useState(false)

  // Load data
  const loadData = useCallback(async () => {
    setIsLoading(true)
    setError(null)
    try {
      const [sub, tierList, pkgList] = await Promise.all([
        subscriptionApi.getInfo(),
        subscriptionApi.getTiers(),
        subscriptionApi.getPackages(),
      ])
      setSubscription(sub)
      setTiers(tierList)
      setPackages(pkgList)
    } catch (e) {
      setError(e instanceof Error ? e.message : '加载失败')
    } finally {
      setIsLoading(false)
    }
  }, [])

  useEffect(() => {
    loadData()
  }, [loadData])

  // Handlers
  const loadHistory = useCallback(async () => {
    setIsLoadingHistory(true)
    try {
      const data = await subscriptionApi.getBalanceHistory()
      setHistoryBalances(data)
    } catch (e) {
      console.error('Failed to load history:', e)
    } finally {
      setIsLoadingHistory(false)
    }
  }, [])

  const handleShowHistory = useCallback(() => {
    setShowHistoryDialog(true)
    loadHistory()
  }, [loadHistory])

  const handleSelectUpgrade = async (tierName: string, subType: 'monthly' | 'yearly') => {
    setSelectedTier(tierName)
    setSelectedSubscriptionType(subType)
    setIsCalculating(true)
    setPriceInfo(null)
    try {
      const result = await subscriptionApi.calculateUpgradePrice(tierName, subType)
      setPriceInfo(result)
    } catch (e) {
      console.error('Calculate price failed:', e)
    } finally {
      setIsCalculating(false)
    }
  }

  const handleUpgrade = async () => {
    if (!selectedTier || !priceInfo?.can_upgrade) return
    setIsProcessing(true)
    setResult(null)
    try {
      const res = await subscriptionApi.upgrade(selectedTier, selectedSubscriptionType)
      setResult({ success: res.success, message: res.message })
      if (res.success) {
        await loadData()
        setTimeout(() => {
          setShowUpgradeDialog(false)
          setResult(null)
          setSelectedTier(null)
          setSelectedSubscriptionType('monthly')
          setPriceInfo(null)
        }, 2000)
      }
    } catch (e) {
      setResult({ success: false, message: e instanceof Error ? e.message : '升级失败' })
    } finally {
      setIsProcessing(false)
    }
  }

  const handlePurchase = async () => {
    if (!selectedPackage) return
    setIsProcessing(true)
    setResult(null)
    try {
      const res = await subscriptionApi.purchase(selectedPackage)
      setResult({ success: res.success, message: res.message })
      if (res.success) {
        await loadData()
        setTimeout(() => {
          setShowCreditDialog(false)
          setResult(null)
          setSelectedPackage(null)
        }, 2000)
      }
    } catch (e) {
      setResult({ success: false, message: e instanceof Error ? e.message : '购买失败' })
    } finally {
      setIsProcessing(false)
    }
  }

  // Helpers
  const formatCredits = (credits: number | string) => {
    const num = Number(credits)
    if (isNaN(num)) return '0'
    return Math.round(num).toLocaleString()
  }

  const formatDate = (dateStr: string) => {
    return new Date(dateStr).toLocaleDateString('zh-CN', {
      year: 'numeric',
      month: 'short',
      day: 'numeric',
    })
  }

  const getUsagePercentage = () => {
    if (!subscription) return 0
    const total = Number(subscription.monthly_credits) + Number(subscription.purchased_credits)
    if (total === 0) return 0
    return Math.min(100, (Number(subscription.used_credits) / total) * 100)
  }

  const tierColors = TIER_COLORS[subscription?.tier || 'free'] || TIER_COLORS.free
  // 显示所有非 free 等级，根据价格判断是否可升级
  const currentTierPrice = tiers.find(t => t.tier_name === subscription?.tier)?.monthly_price ?? 0
  const upgradeTiers = tiers.filter(t => t.tier_name !== 'free')

  return (
    <div className="h-full flex flex-col">
      <PanelHeader
        title={t('订阅与积分')}
        actions={<HeaderMenu route={routes.view.settings('subscription')} />}
      />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-3xl mx-auto">
            {isLoading ? (
              <div className="flex items-center justify-center py-20">
                <RefreshCw className="w-6 h-6 animate-spin text-muted-foreground" />
              </div>
            ) : error ? (
              <div className="text-center py-20">
                <Crown className="w-12 h-12 mx-auto mb-3 opacity-50" />
                <p className="text-muted-foreground">{error}</p>
                <Button variant="outline" className="mt-4" onClick={loadData}>
                  重试
                </Button>
              </div>
            ) : !subscription ? (
              <div className="text-center py-20">
                <Crown className="w-12 h-12 mx-auto mb-3 opacity-50" />
                <p className="text-muted-foreground">请先配置 API Key</p>
              </div>
            ) : (
              <div className="space-y-8">
                {/* Current Subscription */}
                <SettingsSection title={t('当前订阅')}>
                  <SettingsCard className="overflow-hidden">
                    <div className={cn('p-6 rounded-xl', tierColors.bg)}>
                      {/* Header */}
                      <div className="flex items-center justify-between mb-6">
                        <div className="flex items-center gap-4">
                          <div className="w-14 h-14 rounded-2xl bg-white/20 backdrop-blur flex items-center justify-center shadow-lg">
                            <Crown className="w-7 h-7 text-white" />
                          </div>
                          <div>
                            <div className="flex items-center gap-2">
                              <h4 className="text-xl font-bold text-white">
                                {subscription.tier_display_name || subscription.tier}
                              </h4>
                              <span className="px-2 py-0.5 rounded-full bg-white/25 text-white text-xs font-medium">
                                {subscription.subscription_type === 'yearly' ? '年度' : '月度'}
                              </span>
                            </div>
                            <p className={cn('text-sm mt-0.5', tierColors.accent)}>
                              {subscription.status === 'active' ? '✓ 已激活' : '× 已过期'}
                            </p>
                          </div>
                        </div>
                        {/* 只有 ultra 年度订阅不显示升级按钮，其他都可以升级 */}
                        {!(subscription.tier === 'ultra' && subscription.subscription_type === 'yearly') && (
                          <Button 
                            onClick={() => setShowUpgradeDialog(true)}
                            className="bg-white/20 hover:bg-white/30 text-white border-0 backdrop-blur"
                          >
                            <TrendingUp className="w-4 h-4 mr-2" />
                            升级
                          </Button>
                        )}
                      </div>
                      
                      {/* Info Cards */}
                      <div className="grid grid-cols-2 gap-4">
                        <div className={cn('p-4 rounded-xl', tierColors.cardBg)}>
                          <div className="flex items-center gap-2 text-white/70 mb-2">
                            <Calendar className="w-4 h-4" />
                            <span className="text-xs font-medium">计费周期</span>
                          </div>
                          <p className="text-sm font-semibold text-white">
                            {formatDate(subscription.billing_cycle_start)} - {formatDate(subscription.billing_cycle_end)}
                          </p>
                        </div>
                        <div className={cn('p-4 rounded-xl', tierColors.cardBg)}>
                          <div className="flex items-center gap-2 text-white/70 mb-2">
                            <Sparkles className="w-4 h-4" />
                            <span className="text-xs font-medium">每月配额</span>
                          </div>
                          <p className="text-lg font-bold text-white">
                            {formatCredits(subscription.monthly_credits)} <span className="text-sm font-normal text-white/70">积分</span>
                          </p>
                        </div>
                      </div>
                    </div>
                  </SettingsCard>
                </SettingsSection>

                {/* Credits Balance */}
                <SettingsSection title={t('积分余额')}>
                  <SettingsCard>
                    <div className="p-6">
                      <div className="flex items-center justify-between mb-4">
                        <div className="flex items-center gap-3">
                          <div className="w-10 h-10 rounded-xl bg-amber-100 dark:bg-amber-900 flex items-center justify-center">
                            <Coins className="w-5 h-5 text-amber-600 dark:text-amber-400" />
                          </div>
                          <div>
                            <p className="text-sm text-muted-foreground">总可用积分</p>
                            <p className="text-2xl font-bold">{formatCredits(subscription.current_credits)}</p>
                          </div>
                        </div>
                        <Button variant="outline" onClick={() => setShowCreditDialog(true)}>
                          <Zap className="w-4 h-4 mr-2" />
                          购买积分
                        </Button>
                      </div>

                      {/* Balance records list */}
                      {subscription.balances && subscription.balances.length > 0 ? (
                        <div className="space-y-3">
                          {subscription.balances.map((balance) => {
                            const usagePercent = Number(balance.original_amount) > 0 
                              ? (Number(balance.used_amount) / Number(balance.original_amount)) * 100 
                              : 0
                            const isFullyUsed = Number(balance.remaining_amount) === 0
                            return (
                              <div key={balance.id} className={cn('p-3 rounded-lg', isFullyUsed ? 'bg-muted/30' : 'bg-muted/50')}>
                                {/* 第一行：类型 + 过期时间 */}
                                <div className="flex items-center justify-between mb-2">
                                  <div className="flex items-center gap-2">
                                    <span className={cn(
                                      'w-2 h-2 rounded-full',
                                      isFullyUsed ? 'bg-gray-400' :
                                      balance.credit_type === 'monthly' ? 'bg-blue-500' :
                                      balance.credit_type === 'purchased' ? 'bg-green-500' : 'bg-orange-500'
                                    )} />
                                    <span className={cn(
                                      'font-medium text-sm',
                                      isFullyUsed ? 'text-gray-500' :
                                      balance.credit_type === 'monthly' ? 'text-blue-600' :
                                      balance.credit_type === 'purchased' ? 'text-green-600' : 'text-orange-600'
                                    )}>
                                      {balance.credit_type === 'monthly' ? '月度积分' :
                                       balance.credit_type === 'purchased' ? '购买积分' : '赠送积分'}
                                    </span>
                                    {isFullyUsed && (
                                      <span className="text-xs text-gray-400">(已用完)</span>
                                    )}
                                    {balance.description && (
                                      <span className="text-xs text-muted-foreground">- {balance.description}</span>
                                    )}
                                  </div>
                                  {balance.expires_at ? (
                                    <span className="text-xs text-muted-foreground">
                                      {formatDate(balance.expires_at)} 过期
                                    </span>
                                  ) : (
                                    <span className="text-xs text-green-600">永久有效</span>
                                  )}
                                </div>
                                {/* 第二行：进度条 */}
                                <div className="flex items-center gap-3">
                                  <div className="flex-1 h-2 bg-secondary rounded-full overflow-hidden">
                                    <div
                                      className={cn(
                                        'h-full transition-all',
                                        isFullyUsed ? 'bg-gray-400' :
                                        balance.credit_type === 'monthly' ? 'bg-blue-500' :
                                        balance.credit_type === 'purchased' ? 'bg-green-500' : 'bg-orange-500'
                                      )}
                                      style={{ width: `${Math.min(100, usagePercent)}%` }}
                                    />
                                  </div>
                                  <span className="text-xs text-muted-foreground whitespace-nowrap">
                                    {formatCredits(balance.used_amount)} / {formatCredits(balance.original_amount)}
                                  </span>
                                </div>
                              </div>
                            )
                          })}
                        </div>
                      ) : (
                        <div className="text-center py-8 text-muted-foreground">
                          <Coins className="w-8 h-8 mx-auto mb-2 opacity-50" />
                          <p>暂无积分记录</p>
                        </div>
                      )}
                      
                      {/* 查看历史记录入口 */}
                      <button 
                        onClick={handleShowHistory}
                        className="w-full mt-4 p-3 rounded-lg border border-dashed border-muted-foreground/30 text-muted-foreground hover:bg-muted/50 hover:border-muted-foreground/50 transition-colors flex items-center justify-center gap-2"
                      >
                        <History className="w-4 h-4" />
                        <span className="text-sm">查看历史积分记录</span>
                        <ChevronRight className="w-4 h-4" />
                      </button>
                    </div>
                  </SettingsCard>
                </SettingsSection>
              </div>
            )}
          </div>
        </ScrollArea>
      </div>

      {/* Upgrade Dialog */}
      {showUpgradeDialog && (
        <div className="fixed inset-0 z-50 flex items-center justify-center">
          <div className="absolute inset-0 bg-black/50" onClick={() => !isProcessing && setShowUpgradeDialog(false)} />
          <div className="relative bg-background rounded-2xl shadow-xl w-full max-w-2xl max-h-[90vh] overflow-auto m-4">
            <div className="sticky top-0 bg-background border-b px-6 py-4 flex items-center justify-between z-10">
              <div className="flex items-center gap-3">
                <Crown className="w-5 h-5 text-primary" />
                <h2 className="text-lg font-semibold">升级订阅</h2>
              </div>
              <button onClick={() => !isProcessing && setShowUpgradeDialog(false)} className="p-2 rounded-lg hover:bg-foreground/5 transition-colors">
                <X className="w-5 h-5" />
              </button>
            </div>
            <div className="p-6">
              {result ? (
                <div className={cn('p-8 text-center rounded-xl', result.success ? 'bg-green-50 dark:bg-green-900/20' : 'bg-red-50 dark:bg-red-900/20')}>
                  <div className={cn('w-16 h-16 rounded-full mx-auto mb-4 flex items-center justify-center', result.success ? 'bg-green-100 dark:bg-green-900' : 'bg-red-100 dark:bg-red-900')}>
                    {result.success ? <Check className="w-8 h-8 text-green-600" /> : <X className="w-8 h-8 text-red-600" />}
                  </div>
                  <h3 className={cn('text-xl font-semibold mb-2', result.success ? 'text-green-700 dark:text-green-400' : 'text-red-700 dark:text-red-400')}>
                    {result.success ? '升级成功！' : '升级失败'}
                  </h3>
                  <p className="text-muted-foreground">{result.message}</p>
                </div>
              ) : (
                <div className="space-y-6">
                  {/* 说明文字 */}
                  {subscription?.subscription_type === 'yearly' && (
                    <div className="p-3 rounded-lg bg-amber-50 dark:bg-amber-950/30 text-amber-700 dark:text-amber-400 text-sm">
                      👑 您已是年度订阅用户，可以升级到更高等级
                    </div>
                  )}
                  
                  {upgradeTiers.map((tier) => {
                    const tierPrice = Number(tier.monthly_price)
                    const isCurrentTier = tier.tier_name === subscription?.tier
                    const isLowerTier = tierPrice < Number(currentTierPrice)
                    const isYearlyUser = subscription?.subscription_type === 'yearly'
                    
                    // 年度订阅用户不能选月度，只能选更高等级的年度
                    // 低等级不能选择
                    const canSelectMonthly = !isYearlyUser && !isLowerTier && !isCurrentTier
                    const canSelectYearly = tier.yearly_price && !isLowerTier && !(isCurrentTier && isYearlyUser)
                    
                    const discount = tier.yearly_discount
                      ? Math.round((1 - Number(tier.yearly_discount)) * 100)
                      : null
                    
                    return (
                      <div key={tier.id} className="border rounded-xl overflow-hidden">
                        {/* 等级标题 */}
                        <div className="px-5 py-3 bg-muted/50 flex items-center justify-between">
                          <div className="flex items-center gap-2">
                            <h4 className="font-semibold text-base">{tier.display_name}</h4>
                            {isCurrentTier && (
                              <span className="px-2 py-0.5 rounded-full bg-blue-500 text-white text-xs">当前</span>
                            )}
                            {isLowerTier && (
                              <span className="px-2 py-0.5 rounded-full bg-gray-400 text-white text-xs">低于当前</span>
                            )}
                          </div>
                          <span className="text-sm text-muted-foreground">
                            {formatCredits(tier.monthly_credits)} 积分/月
                          </span>
                        </div>
                        
                        {/* 订阅选项 */}
                        <div className="p-4 grid grid-cols-2 gap-3">
                          {/* 月度订阅 */}
                          <button
                            onClick={() => canSelectMonthly && handleSelectUpgrade(tier.tier_name, 'monthly')}
                            disabled={!canSelectMonthly}
                            className={cn(
                              'p-4 rounded-lg border-2 text-left transition-all relative',
                              !canSelectMonthly
                                ? 'opacity-40 cursor-not-allowed border-border bg-muted/30'
                                : 'cursor-pointer hover:border-blue-300',
                              selectedTier === tier.tier_name && selectedSubscriptionType === 'monthly'
                                ? 'ring-2 ring-blue-400 border-blue-400 bg-blue-50 dark:bg-blue-950/30'
                                : 'border-border'
                            )}
                          >
                            <div className="text-xs text-muted-foreground mb-1">月度订阅</div>
                            <div className="text-xl font-bold">
                              ¥{tier.monthly_price}<span className="text-sm font-normal text-muted-foreground">/月</span>
                            </div>
                            {!canSelectMonthly && isYearlyUser && !isLowerTier && (
                              <div className="text-xs text-orange-500 mt-1">年度用户不可选</div>
                            )}
                            {!canSelectMonthly && !isYearlyUser && isCurrentTier && (
                              <div className="text-xs text-blue-500 mt-1">当前订阅</div>
                            )}
                          </button>
                          
                          {/* 年度订阅 */}
                          <button
                            onClick={() => canSelectYearly && handleSelectUpgrade(tier.tier_name, 'yearly')}
                            disabled={!canSelectYearly}
                            className={cn(
                              'p-4 rounded-lg border-2 text-left transition-all relative',
                              !canSelectYearly
                                ? 'opacity-40 cursor-not-allowed border-border bg-muted/30'
                                : 'cursor-pointer hover:border-green-300',
                              selectedTier === tier.tier_name && selectedSubscriptionType === 'yearly'
                                ? 'ring-2 ring-green-400 border-green-400 bg-green-50 dark:bg-green-950/30'
                                : 'border-border'
                            )}
                          >
                            {discount && discount > 0 && canSelectYearly && (
                              <span className="absolute -top-2 -right-2 px-2 py-0.5 rounded-full bg-green-500 text-white text-xs">省{discount}%</span>
                            )}
                            <div className="text-xs text-muted-foreground mb-1">年度订阅</div>
                            {tier.yearly_price ? (
                              <>
                                <div className="text-xl font-bold text-green-600">
                                  ¥{tier.yearly_price}<span className="text-sm font-normal text-muted-foreground">/年</span>
                                </div>
                                <div className="text-xs text-green-600">约 ¥{Math.round(Number(tier.yearly_price) / 12)}/月</div>
                                {!canSelectYearly && isCurrentTier && isYearlyUser && (
                                  <div className="text-xs text-blue-500 mt-1">当前订阅</div>
                                )}
                              </>
                            ) : (
                              <div className="text-sm text-muted-foreground">暂不支持</div>
                            )}
                          </button>
                        </div>
                      </div>
                    )
                  })}
                </div>
              )}
            </div>
            {!result && upgradeTiers.length > 0 && (
              <div className="sticky bottom-0 bg-background border-t px-6 py-4">
                {/* 价格详情 */}
                {selectedTier && (
                  <div className="mb-4 p-4 rounded-lg bg-muted/50">
                    {isCalculating ? (
                      <div className="flex items-center justify-center gap-2 text-muted-foreground">
                        <Loader2 className="w-4 h-4 animate-spin" />
                        <span>计算中...</span>
                      </div>
                    ) : priceInfo ? (
                      priceInfo.can_upgrade ? (
                        <div className="space-y-2">
                          <div className="flex items-center justify-between text-sm">
                            <span className="text-muted-foreground">
                              升级到 {priceInfo.target_tier_display} ({priceInfo.subscription_type === 'yearly' ? '年度' : '月度'})
                            </span>
                            <span>原价 ¥{priceInfo.original_price}</span>
                          </div>
                          {priceInfo.remaining_value > 0 && (
                            <div className="flex items-center justify-between text-sm">
                              <span className="text-muted-foreground">
                                当前订阅剩余价值 ({priceInfo.remaining_days}天)
                              </span>
                              <span className="text-green-600">- ¥{priceInfo.remaining_value}</span>
                            </div>
                          )}
                          <div className="flex items-center justify-between font-semibold border-t pt-2">
                            <span>实际支付</span>
                            <span className="text-xl text-red-500">¥{priceInfo.final_price}</span>
                          </div>
                        </div>
                      ) : (
                        <div className="text-center text-orange-500">
                          {priceInfo.message}
                        </div>
                      )
                    ) : null}
                  </div>
                )}
                
                <div className="flex justify-end gap-3">
                  <Button variant="outline" onClick={() => {
                    setShowUpgradeDialog(false)
                    setPriceInfo(null)
                    setSelectedTier(null)
                  }} disabled={isProcessing}>
                    取消
                  </Button>
                  <Button 
                    onClick={handleUpgrade} 
                    disabled={!selectedTier || !priceInfo?.can_upgrade || isProcessing || isCalculating}
                    className="!bg-gradient-to-r !from-amber-500 !to-yellow-500 hover:!from-amber-600 hover:!to-yellow-600 !text-white !border-0"
                  >
                    {isProcessing ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Crown className="w-4 h-4 mr-2" />}
                    {priceInfo?.final_price !== undefined ? `支付 ¥${priceInfo.final_price}` : '确认升级'}
                  </Button>
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Credit Balance History Dialog */}
      {showHistoryDialog && (
        <div className="fixed inset-0 z-50 flex items-center justify-center">
          <div className="absolute inset-0 bg-black/50" onClick={() => setShowHistoryDialog(false)} />
          <div className="relative bg-background rounded-2xl shadow-xl w-full max-w-lg max-h-[80vh] overflow-auto m-4">
            <div className="sticky top-0 bg-background border-b px-6 py-4 flex items-center justify-between">
              <div className="flex items-center gap-3">
                <History className="w-5 h-5 text-muted-foreground" />
                <h2 className="text-lg font-semibold">历史积分记录</h2>
              </div>
              <button onClick={() => setShowHistoryDialog(false)} className="p-2 rounded-lg hover:bg-foreground/5 transition-colors">
                <X className="w-5 h-5" />
              </button>
            </div>
            <div className="p-6">
              {isLoadingHistory ? (
                <div className="flex items-center justify-center py-8">
                  <Loader2 className="w-6 h-6 animate-spin text-muted-foreground" />
                </div>
              ) : historyBalances.length === 0 ? (
                <div className="text-center py-8 text-muted-foreground">
                  暂无历史积分记录
                </div>
              ) : (
                <div className="space-y-3">
                  {historyBalances.map((balance) => (
                    <div key={balance.id} className="p-4 rounded-xl bg-muted/30 border border-border/50">
                      <div className="flex justify-between items-start mb-2">
                        <span className="text-sm font-medium">{balance.source_description || balance.source_type}</span>
                        <span className="text-xs px-2 py-0.5 rounded-full bg-muted text-muted-foreground">已过期</span>
                      </div>
                      <div className="flex items-baseline gap-2 mb-2">
                        <span className="text-lg font-semibold text-muted-foreground line-through">
                          {formatCredits(balance.remaining_credits)}/{formatCredits(balance.initial_credits)}
                        </span>
                        <span className="text-xs text-muted-foreground">积分</span>
                      </div>
                      <div className="text-xs text-muted-foreground">
                        过期时间: {new Date(balance.expires_at).toLocaleDateString('zh-CN')}
                      </div>
                    </div>
                  ))}
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Credit Package Dialog */}
      {showCreditDialog && (
        <div className="fixed inset-0 z-50 flex items-center justify-center">
          <div className="absolute inset-0 bg-black/50" onClick={() => !isProcessing && setShowCreditDialog(false)} />
          <div className="relative bg-background rounded-2xl shadow-xl w-full max-w-2xl max-h-[90vh] overflow-auto m-4">
            <div className="sticky top-0 bg-background border-b px-6 py-4 flex items-center justify-between">
              <div className="flex items-center gap-3">
                <Coins className="w-5 h-5 text-amber-600" />
                <h2 className="text-lg font-semibold">购买积分包</h2>
              </div>
              <button onClick={() => !isProcessing && setShowCreditDialog(false)} className="p-2 rounded-lg hover:bg-foreground/5 transition-colors">
                <X className="w-5 h-5" />
              </button>
            </div>
            <div className="p-6">
              {result ? (
                <div className={cn('p-8 text-center rounded-xl', result.success ? 'bg-green-50 dark:bg-green-900/20' : 'bg-red-50 dark:bg-red-900/20')}>
                  <div className={cn('w-16 h-16 rounded-full mx-auto mb-4 flex items-center justify-center', result.success ? 'bg-green-100 dark:bg-green-900' : 'bg-red-100 dark:bg-red-900')}>
                    {result.success ? <Check className="w-8 h-8 text-green-600" /> : <X className="w-8 h-8 text-red-600" />}
                  </div>
                  <h3 className={cn('text-xl font-semibold mb-2', result.success ? 'text-green-700 dark:text-green-400' : 'text-red-700 dark:text-red-400')}>
                    {result.success ? '购买成功！' : '购买失败'}
                  </h3>
                  <p className="text-muted-foreground">{result.message}</p>
                </div>
              ) : (
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                  {packages.map((pkg) => (
                    <button
                      key={pkg.id}
                      onClick={() => setSelectedPackage(pkg.id)}
                      className={cn(
                        'p-5 rounded-xl border-2 text-left transition-all relative cursor-pointer',
                        selectedPackage === pkg.id 
                          ? 'ring-2 ring-amber-400 border-amber-400 bg-gradient-to-br from-amber-50 to-yellow-50 dark:from-amber-950/30 dark:to-yellow-950/30' 
                          : 'border-border hover:border-amber-300'
                      )}
                    >
                      {pkg.bonus_credits > 0 && (
                        <div className="absolute -top-2 -right-2 px-2 py-1 rounded-full bg-red-500 text-white text-xs font-medium flex items-center gap-1">
                          <Gift className="w-3 h-3" />
                          +{formatCredits(pkg.bonus_credits)}
                        </div>
                      )}
                      <div className="flex justify-between items-center mb-2">
                        <h4 className="font-semibold">{pkg.package_name}</h4>
                        <span className="text-lg font-bold text-red-500">¥{pkg.price}</span>
                      </div>
                      <div className="flex items-baseline gap-2">
                        <span className="text-2xl font-bold">{formatCredits(pkg.credits)}</span>
                        {pkg.bonus_credits > 0 && <span className="text-green-600">+{formatCredits(pkg.bonus_credits)}</span>}
                        <span className="text-muted-foreground">积分</span>
                      </div>
                      {pkg.description && <p className="text-sm text-muted-foreground mt-2">{pkg.description}</p>}
                    </button>
                  ))}
                </div>
              )}
            </div>
            {!result && packages.length > 0 && (
              <div className="sticky bottom-0 bg-background border-t px-6 py-4 flex justify-end gap-3">
                <Button variant="outline" onClick={() => setShowCreditDialog(false)} disabled={isProcessing}>
                  取消
                </Button>
                <Button 
                  onClick={handlePurchase} 
                  disabled={!selectedPackage || isProcessing}
                  className="!bg-gradient-to-r !from-amber-500 !to-yellow-500 hover:!from-amber-600 hover:!to-yellow-600 !text-white !border-0"
                >
                  {isProcessing ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Coins className="w-4 h-4 mr-2" />}
                  确认购买
                </Button>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  )
}
