import { useState, useEffect, useMemo } from "react"
import { ArrowLeft, Check, Cloud, Package, RefreshCw } from "lucide-react"
import { cn } from "@/lib/utils"
import { AddWorkspacePrimaryButton } from "./primitives"
import { useT } from "@/context/LocaleContext"
import { CrossfadeAvatar } from "@/components/ui/avatar"
import { ScrollArea } from "@/components/ui/scroll-area"
import type { MarketplaceApp } from "@creator-flow/shared/marketplace"

interface AppOption {
  id: string
  name: string
  description: string
  icon: string
  iconUrl?: string
  source: 'bundled' | 'marketplace'
}

interface AddWorkspaceStep_ChooseAppProps {
  onBack?: () => void
  onNext: (appId: string, appName: string, source: 'bundled' | 'marketplace') => void
  isLoading?: boolean
  /** Whether this is the first step (hides back button) */
  isFirstStep?: boolean
}

interface AppCardProps {
  app: AppOption
  selected: boolean
  onClick: () => void
  disabled?: boolean
}

function AppCard({ app, selected, onClick, disabled }: AppCardProps) {
  return (
    <button
      onClick={onClick}
      disabled={disabled}
      className={cn(
        "relative flex items-start gap-3 p-3 rounded-xl text-left w-full",
        "bg-background border border-border/60 shadow-sm",
        "transition-all duration-150",
        "focus:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2",
        selected
          ? "ring-2 ring-accent border-accent/50"
          : "hover:border-border hover:shadow-md",
        disabled && "opacity-50 cursor-not-allowed"
      )}
    >
      {/* Icon */}
      <div className={cn(
        "flex h-10 w-10 shrink-0 items-center justify-center rounded-xl text-xl",
        selected
          ? "bg-accent/10"
          : "bg-foreground/5"
      )}>
        {app.iconUrl ? (
          <CrossfadeAvatar
            src={app.iconUrl}
            alt={app.name}
            className="h-8 w-8 rounded-lg"
            fallbackClassName="bg-foreground/5 text-foreground/60 text-base rounded-lg"
            fallback={app.icon}
          />
        ) : (
          app.icon
        )}
      </div>

      {/* Content */}
      <div className="flex-1 min-w-0">
        <div className="font-medium text-sm text-foreground leading-tight">
          {app.name}
        </div>
        <div className="text-xs text-muted-foreground mt-0.5 line-clamp-2">
          {app.description}
        </div>
      </div>

      {/* Source badge */}
      {app.source === 'marketplace' && (
        <div className="shrink-0">
          <Cloud className="h-3.5 w-3.5 text-muted-foreground" />
        </div>
      )}

      {/* Selected indicator */}
      {selected && (
        <div className="absolute top-2 right-2">
          <div className="flex h-5 w-5 items-center justify-center rounded-full bg-accent text-accent-foreground">
            <Check className="h-3 w-3" />
          </div>
        </div>
      )}
    </button>
  )
}

/**
 * AddWorkspaceStep_ChooseApp - Choose an application template for the workspace
 *
 * Displays available applications (bundled and marketplace) and lets user choose one.
 * Uses a compact 4-column grid layout.
 */
export function AddWorkspaceStep_ChooseApp({
  onBack,
  onNext,
  isLoading,
  isFirstStep = true
}: AddWorkspaceStep_ChooseAppProps) {
  const t = useT()
  const [selectedAppId, setSelectedAppId] = useState<string>('app.general')
  const [bundledApps, setBundledApps] = useState<AppOption[]>([])
  const [marketplaceApps, setMarketplaceApps] = useState<AppOption[]>([])
  const [loading, setLoading] = useState(true)
  const [loadingMarketplace, setLoadingMarketplace] = useState(true)

  // Load bundled apps from main process
  useEffect(() => {
    const loadBundledApps = async () => {
      try {
        const apps = await window.electronAPI.listBundledApps()
        const appOptions: AppOption[] = apps.map(app => ({
          id: app.id,
          name: app.name,
          description: app.description,
          icon: getAppIcon(app.id),
          source: 'bundled' as const
        }))
        setBundledApps(appOptions)
      } catch (error) {
        console.error('Failed to load bundled apps:', error)
        setBundledApps([{
          id: 'app.general',
          name: t('通用工作区'),
          description: t('通用 AI 助手工作区，适合各种任务场景'),
          icon: '🤖',
          source: 'bundled'
        }])
      } finally {
        setLoading(false)
      }
    }
    loadBundledApps()
  }, [t])

  // Load marketplace apps
  useEffect(() => {
    const loadMarketplaceApps = async () => {
      try {
        const result = await window.electronAPI.marketplaceListApps()
        const appOptions: AppOption[] = result.items.map((app: MarketplaceApp) => ({
          id: `marketplace:${app.app_id}`,
          name: app.name,
          description: app.description || '',
          icon: '📦',
          iconUrl: app.icon_url || undefined,
          source: 'marketplace' as const
        }))
        setMarketplaceApps(appOptions)
      } catch (error) {
        console.error('Failed to load marketplace apps:', error)
      } finally {
        setLoadingMarketplace(false)
      }
    }
    loadMarketplaceApps()
  }, [])

  // Get default icon for bundled app based on ID
  function getAppIcon(appId: string): string {
    const iconMap: Record<string, string> = {
      'app.general': '🤖',
      'app.creator-media': '✨',
      'app.software-dev': '💻',
      'app.personal-assistant': '🗓️',
      'app.smart-office': '📊',
      'app.writing-assistant': '✍️',
    }
    return iconMap[appId] || '📦'
  }

  // Combine all apps
  const allApps = useMemo(() => {
    return [...bundledApps, ...marketplaceApps]
  }, [bundledApps, marketplaceApps])

  // Find selected app
  const selectedApp = useMemo(() => {
    return allApps.find(app => app.id === selectedAppId)
  }, [allApps, selectedAppId])

  const handleNext = () => {
    if (selectedApp) {
      onNext(selectedApp.id, selectedApp.name, selectedApp.source)
    }
  }

  const canProceed = selectedAppId && !loading && !isLoading

  return (
    <div className="flex flex-col w-full h-full max-w-4xl mx-auto">
      {/* Header */}
      <div className="text-center mb-6">
        <h1 className="text-xl font-semibold tracking-tight">
          {t('添加工作区')}
        </h1>
        <p className="mt-1 text-sm text-muted-foreground">
          {t('选择适合你的工作流的应用模板。')}
        </p>
      </div>

      {/* Apps Grid - Scrollable */}
      <ScrollArea className="flex-1 -mx-4 px-4">
        <div className="pb-4">
          {loading ? (
            <div className="flex items-center justify-center py-12">
              <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
            </div>
          ) : (
            <>
              {/* Bundled Apps Section */}
              {bundledApps.length > 0 && (
                <div className="mb-6">
                  <div className="flex items-center gap-2 mb-3">
                    <Package className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium text-muted-foreground">{t('内置应用')}</span>
                  </div>
                  <div className="grid grid-cols-2 md:grid-cols-3 gap-2">
                    {bundledApps.map((app) => (
                      <AppCard
                        key={app.id}
                        app={app}
                        selected={selectedAppId === app.id}
                        onClick={() => setSelectedAppId(app.id)}
                        disabled={isLoading}
                      />
                    ))}
                  </div>
                </div>
              )}

              {/* Marketplace Apps Section */}
              {loadingMarketplace ? (
                <div className="flex items-center justify-center py-8">
                  <RefreshCw className="h-4 w-4 animate-spin text-muted-foreground mr-2" />
                  <span className="text-sm text-muted-foreground">{t('加载市场应用...')}</span>
                </div>
              ) : marketplaceApps.length > 0 && (
                <div>
                  <div className="flex items-center gap-2 mb-3">
                    <Cloud className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium text-muted-foreground">{t('市场应用')}</span>
                  </div>
                  <div className="grid grid-cols-2 md:grid-cols-3 gap-2">
                    {marketplaceApps.map((app) => (
                      <AppCard
                        key={app.id}
                        app={app}
                        selected={selectedAppId === app.id}
                        onClick={() => setSelectedAppId(app.id)}
                        disabled={isLoading}
                      />
                    ))}
                  </div>
                </div>
              )}
            </>
          )}
        </div>
      </ScrollArea>

      {/* Footer with Next button */}
      <div className="pt-4 mt-auto">
        <AddWorkspacePrimaryButton
          onClick={handleNext}
          disabled={!canProceed}
          loading={isLoading}
          loadingText={t('加载中...')}
          className="w-full max-w-md mx-auto block"
        >
          {t('下一步')}
        </AddWorkspacePrimaryButton>
      </div>
    </div>
  )
}
