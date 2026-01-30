import { useState, useEffect } from "react"
import { ArrowLeft, Sparkles, Briefcase, Check } from "lucide-react"
import { cn } from "@/lib/utils"
import { AddWorkspaceContainer, AddWorkspaceStepHeader, AddWorkspacePrimaryButton } from "./primitives"
import { useT } from "@/context/LocaleContext"

interface AppOption {
  id: string
  name: string
  description: string
  icon: string
}

interface AddWorkspaceStep_ChooseAppProps {
  onBack?: () => void
  onNext: (appId: string, appName: string) => void
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
        "relative flex items-start gap-4 w-full p-4 rounded-lg text-left",
        "bg-background shadow-minimal",
        "transition-all duration-150",
        "focus:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2",
        selected
          ? "ring-2 ring-accent bg-accent/5"
          : "hover:bg-foreground/5",
        disabled && "opacity-50 cursor-not-allowed"
      )}
    >
      {/* Icon */}
      <div className={cn(
        "flex h-12 w-12 shrink-0 items-center justify-center rounded-lg text-2xl",
        selected
          ? "bg-accent/10"
          : "bg-foreground/5"
      )}>
        {app.icon}
      </div>

      {/* Content */}
      <div className="flex-1 min-w-0">
        <div className="font-medium text-[15px] text-foreground mb-1">
          {app.name}
        </div>
        <div className="text-[13px] text-muted-foreground leading-relaxed">
          {app.description}
        </div>
      </div>

      {/* Selected indicator */}
      {selected && (
        <div className="absolute top-3 right-3">
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
 * Displays available applications (bundled and installed) and lets user choose one.
 */
export function AddWorkspaceStep_ChooseApp({
  onBack,
  onNext,
  isLoading,
  isFirstStep = true
}: AddWorkspaceStep_ChooseAppProps) {
  const t = useT()
  const [selectedAppId, setSelectedAppId] = useState<string>('app.general')
  const [apps, setApps] = useState<AppOption[]>([])
  const [loading, setLoading] = useState(true)

  // Load available apps from main process
  useEffect(() => {
    const loadApps = async () => {
      try {
        const bundledApps = await window.electronAPI.listBundledApps()
        // Map to AppOption format with default icons for apps without iconPath
        const appOptions: AppOption[] = bundledApps.map(app => ({
          id: app.id,
          name: app.name,
          description: app.description,
          icon: getAppIcon(app.id) // Use emoji icon based on app ID
        }))
        setApps(appOptions)
      } catch (error) {
        console.error('Failed to load apps:', error)
        // Fallback to general app if loading fails
        setApps([{
          id: 'app.general',
          name: t('通用工作区'),
          description: t('通用 AI 助手工作区，适合各种任务场景'),
          icon: '🤖'
        }])
      } finally {
        setLoading(false)
      }
    }

    loadApps()
  }, [t])

  // Get default icon for app based on ID
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

  const handleNext = () => {
    if (selectedAppId) {
      const selectedApp = apps.find(app => app.id === selectedAppId)
      onNext(selectedAppId, selectedApp?.name || '')
    }
  }

  const canProceed = selectedAppId && !loading && !isLoading

  return (
    <AddWorkspaceContainer>
      {/* Back button - only show if not first step */}
      {!isFirstStep && onBack && (
        <button
          onClick={onBack}
          disabled={isLoading}
          className={cn(
            "self-start flex items-center gap-1 text-sm text-muted-foreground",
            "hover:text-foreground transition-colors mb-4",
            isLoading && "opacity-50 cursor-not-allowed"
          )}
        >
          <ArrowLeft className="h-4 w-4" />
          {t('返回')}
        </button>
      )}

      {isFirstStep && <div className="mt-2" />}

      <AddWorkspaceStepHeader
        title={t('添加工作区')}
        description={t('选择适合你的工作流的应用模板。')}
      />

      <div className="mt-6 w-full space-y-3">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <div className="text-sm text-muted-foreground">{t('加载中...')}</div>
          </div>
        ) : (
          <>
            {apps.map((app) => (
              <AppCard
                key={app.id}
                app={app}
                selected={selectedAppId === app.id}
                onClick={() => setSelectedAppId(app.id)}
                disabled={isLoading}
              />
            ))}
          </>
        )}
      </div>

      {/* Next button */}
      <div className="mt-6 w-full">
        <AddWorkspacePrimaryButton
          onClick={handleNext}
          disabled={!canProceed}
          loading={isLoading}
          loadingText={t('加载中...')}
        >
          {t('下一步')}
        </AddWorkspacePrimaryButton>
      </div>
    </AddWorkspaceContainer>
  )
}
