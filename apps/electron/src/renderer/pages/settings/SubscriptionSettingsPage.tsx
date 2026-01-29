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
import { Crown, Coins, TrendingUp, Calendar, Zap, RefreshCw, Check, X, Loader2, Gift, Sparkles } from 'lucide-react'
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
} from '@/api/subscription'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'subscription',
}

// ============================================
// Constants
// ============================================

const TIER_COLORS: Record<string, { bg: string; text: string }> = {
  free: { bg: 'bg-gray-100 dark:bg-gray-800', text: 'text-gray-600 dark:text-gray-400' },
  pro: { bg: 'bg-blue-100 dark:bg-blue-900', text: 'text-blue-600 dark:text-blue-400' },
  max: { bg: 'bg-purple-100 dark:bg-purple-900', text: 'text-purple-600 dark:text-purple-400' },
  ultra: { bg: 'bg-amber-100 dark:bg-amber-900', text: 'text-amber-600 dark:text-amber-400' },
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
  const [selectedTier, setSelectedTier] = useState<string | null>(null)
  const [selectedPackage, setSelectedPackage] = useState<number | null>(null)
  const [isProcessing, setIsProcessing] = useState(false)
  const [result, setResult] = useState<{ success: boolean; message: string } | null>(null)

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
  const handleUpgrade = async () => {
    if (!selectedTier) return
    setIsProcessing(true)
    setResult(null)
    try {
      const res = await subscriptionApi.upgrade(selectedTier)
      setResult({ success: res.success, message: res.message })
      if (res.success) {
        await loadData()
        setTimeout(() => {
          setShowUpgradeDialog(false)
          setResult(null)
          setSelectedTier(null)
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
                  <SettingsCard>
                    <div className={cn('p-6 rounded-xl', tierColors.bg)}>
                      <div className="flex items-center justify-between mb-4">
                        <div className="flex items-center gap-3">
                          <div className={cn('w-12 h-12 rounded-xl flex items-center justify-center', tierColors.bg)}>
                            <Crown className={cn('w-6 h-6', tierColors.text)} />
                          </div>
                          <div>
                            <h4 className={cn('text-lg font-bold', tierColors.text)}>
                              {subscription.tier_display_name || subscription.tier}
                            </h4>
                            <p className="text-sm text-muted-foreground">
                              {subscription.status === 'active' ? '已激活' : '已过期'}
                            </p>
                          </div>
                        </div>
                        {subscription.tier !== 'ultra' && (
                          <Button onClick={() => setShowUpgradeDialog(true)}>
                            <TrendingUp className="w-4 h-4 mr-2" />
                            升级
                          </Button>
                        )}
                      </div>
                      <div className="grid grid-cols-2 gap-4">
                        <div className="p-3 rounded-lg bg-background/50">
                          <div className="flex items-center gap-2 text-muted-foreground mb-1">
                            <Calendar className="w-4 h-4" />
                            <span className="text-xs">计费周期</span>
                          </div>
                          <p className="text-sm font-medium">
                            {formatDate(subscription.billing_cycle_start)} - {formatDate(subscription.billing_cycle_end)}
                          </p>
                        </div>
                        <div className="p-3 rounded-lg bg-background/50">
                          <div className="flex items-center gap-2 text-muted-foreground mb-1">
                            <Sparkles className="w-4 h-4" />
                            <span className="text-xs">每月配额</span>
                          </div>
                          <p className="text-sm font-medium">{formatCredits(subscription.monthly_credits)} 积分</p>
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
                            return (
                              <div key={balance.id} className="p-3 rounded-lg bg-muted/50">
                                {/* 第一行：类型 + 过期时间 */}
                                <div className="flex items-center justify-between mb-2">
                                  <div className="flex items-center gap-2">
                                    <span className={cn(
                                      'w-2 h-2 rounded-full',
                                      balance.credit_type === 'monthly' ? 'bg-blue-500' :
                                      balance.credit_type === 'purchased' ? 'bg-green-500' : 'bg-orange-500'
                                    )} />
                                    <span className={cn(
                                      'font-medium text-sm',
                                      balance.credit_type === 'monthly' ? 'text-blue-600' :
                                      balance.credit_type === 'purchased' ? 'text-green-600' : 'text-orange-600'
                                    )}>
                                      {balance.credit_type === 'monthly' ? '月度积分' :
                                       balance.credit_type === 'purchased' ? '购买积分' : '赠送积分'}
                                    </span>
                                    {balance.description && (
                                      <span className="text-xs text-muted-foreground">- {balance.description}</span>
                                    )}
                                  </div>
                                  {balance.expires_at ? (
                                    <span className={cn(
                                      'text-xs',
                                      new Date(balance.expires_at) < new Date() ? 'text-red-500' : 'text-muted-foreground'
                                    )}>
                                      {new Date(balance.expires_at) < new Date() ? '已过期' : `${formatDate(balance.expires_at)} 过期`}
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
          <div className="relative bg-background rounded-2xl shadow-xl w-full max-w-3xl max-h-[90vh] overflow-auto m-4">
            <div className="sticky top-0 bg-background border-b px-6 py-4 flex items-center justify-between">
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
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  {upgradeTiers.map((tier) => {
                    const isCurrentTier = tier.tier_name === subscription?.tier
                    const canUpgrade = Number(tier.monthly_price) > Number(currentTierPrice)
                    const isDisabled = !canUpgrade || isCurrentTier
                    return (
                      <button
                        key={tier.id}
                        onClick={() => !isDisabled && setSelectedTier(tier.tier_name)}
                        disabled={isDisabled}
                        className={cn(
                          'p-5 rounded-xl border-2 text-left transition-all relative',
                          isDisabled 
                            ? 'opacity-50 cursor-not-allowed border-border' 
                            : 'cursor-pointer hover:border-amber-300',
                          selectedTier === tier.tier_name 
                            ? 'ring-2 ring-amber-400 border-amber-400 bg-gradient-to-br from-amber-50 to-yellow-50 dark:from-amber-950/30 dark:to-yellow-950/30' 
                            : 'border-border'
                        )}
                      >
                        {isCurrentTier && (
                          <span className="absolute -top-2 -right-2 px-2 py-0.5 rounded-full bg-blue-500 text-white text-xs">当前</span>
                        )}
                        <h4 className="font-semibold">{tier.display_name}</h4>
                        <p className="text-2xl font-bold text-red-500">¥{tier.monthly_price}<span className="text-base">/月</span></p>
                        <p className="text-sm text-muted-foreground mt-2">{formatCredits(tier.monthly_credits)} 积分/月</p>
                      </button>
                    )
                  })}
                </div>
              )}
            </div>
            {!result && upgradeTiers.length > 0 && (
              <div className="sticky bottom-0 bg-background border-t px-6 py-4 flex justify-end gap-3">
                <Button variant="outline" onClick={() => setShowUpgradeDialog(false)} disabled={isProcessing}>
                  取消
                </Button>
                <Button 
                  onClick={handleUpgrade} 
                  disabled={!selectedTier || isProcessing}
                  className="!bg-gradient-to-r !from-amber-500 !to-yellow-500 hover:!from-amber-600 hover:!to-yellow-600 !text-white !border-0"
                >
                  {isProcessing ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Crown className="w-4 h-4 mr-2" />}
                  确认升级
                </Button>
              </div>
            )}
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
