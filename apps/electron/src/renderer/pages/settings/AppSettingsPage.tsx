/**
 * AppSettingsPage
 *
 * Global app-level settings that apply across all workspaces.
 *
 * Settings:
 * - Appearance (Theme, Font)
 * - Notifications
 * - API Connection (opens OnboardingWizard for editing)
 */

import { useState, useEffect, useCallback } from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { useTheme } from '@/context/ThemeContext'
import { useLocale, useT } from '@/context/LocaleContext'
import { routes } from '@/lib/navigate'
import { Monitor, Sun, Moon, X, Globe } from 'lucide-react'
import { Spinner, FullscreenOverlayBase } from '@creator-flow/ui'
import { useSetAtom } from 'jotai'
import { fullscreenOverlayOpenAtom } from '@/atoms/overlay'
import type { AuthType } from '../../../shared/types'
import type { DetailsPageMeta } from '@/lib/navigation-registry'

import {
  SettingsSection,
  SettingsCard,
  SettingsRow,
  SettingsToggle,
  SettingsSegmentedControl,
  SettingsMenuSelect,
} from '@/components/settings'
import { useUpdateChecker } from '@/hooks/useUpdateChecker'
import { useOnboarding } from '@/hooks/useOnboarding'
import { OnboardingWizard } from '@/components/onboarding'
import { useAppShellContext } from '@/context/AppShellContext'
import type { PresetTheme } from '@config/theme'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'app',
}

// ============================================
// Main Component
// ============================================

export default function AppSettingsPage() {
  const { mode, setMode, colorTheme, setColorTheme, font, setFont } = useTheme()
  const { locale, setLocale, languages } = useLocale()
  const t = useT()
  const { refreshCustomModel } = useAppShellContext()

  // Preset themes state
  const [presetThemes, setPresetThemes] = useState<PresetTheme[]>([])

  // API Connection state (read-only display — editing is done via OnboardingWizard overlay)
  const [authType, setAuthType] = useState<AuthType>('api_key')
  const [hasCredential, setHasCredential] = useState(false)
  const [showApiSetup, setShowApiSetup] = useState(false)
  const setFullscreenOverlayOpen = useSetAtom(fullscreenOverlayOpenAtom)

  // Notifications state
  const [notificationsEnabled, setNotificationsEnabled] = useState(true)

  // Auto-update state
  const updateChecker = useUpdateChecker()
  const [isCheckingForUpdates, setIsCheckingForUpdates] = useState(false)

  const handleCheckForUpdates = useCallback(async () => {
    setIsCheckingForUpdates(true)
    try {
      await updateChecker.checkForUpdates()
    } finally {
      setIsCheckingForUpdates(false)
    }
  }, [updateChecker])

  // Load current API connection info and notifications on mount
  const loadConnectionInfo = useCallback(async () => {
    if (!window.electronAPI) return
    try {
      const [billing, notificationsOn] = await Promise.all([
        window.electronAPI.getApiSetup(),
        window.electronAPI.getNotificationsEnabled(),
      ])
      setAuthType(billing.authType)
      setHasCredential(billing.hasCredential)
      setNotificationsEnabled(notificationsOn)
    } catch (error) {
      console.error('Failed to load settings:', error)
    }
  }, [])

  useEffect(() => {
    loadConnectionInfo()
  }, [])

  // Load preset themes when workspace changes (themes are workspace-scoped)
  // Load preset themes (app-level, no workspace dependency)
  useEffect(() => {
    const loadThemes = async () => {
      if (!window.electronAPI) {
        setPresetThemes([])
        return
      }
      try {
        const themes = await window.electronAPI.loadPresetThemes()
        setPresetThemes(themes)
      } catch (error) {
        console.error('Failed to load preset themes:', error)
        setPresetThemes([])
      }
    }
    loadThemes()
  }, [])

  // Helpers to open/close the fullscreen API setup overlay
  const openApiSetup = useCallback(() => {
    setShowApiSetup(true)
    setFullscreenOverlayOpen(true)
  }, [setFullscreenOverlayOpen])

  const closeApiSetup = useCallback(() => {
    setShowApiSetup(false)
    setFullscreenOverlayOpen(false)
  }, [setFullscreenOverlayOpen])

  // OnboardingWizard hook for editing API connection (starts at api-setup step).
  // onConfigSaved fires immediately when billing is persisted, updating the model UI instantly.
  const apiSetupOnboarding = useOnboarding({
    initialStep: 'api-setup',
    onConfigSaved: refreshCustomModel,
    onComplete: () => {
      closeApiSetup()
      loadConnectionInfo()
      apiSetupOnboarding.reset()
    },
    onDismiss: () => {
      closeApiSetup()
      apiSetupOnboarding.reset()
    },
  })

  // Called when user completes the wizard (clicks Finish on completion step)
  const handleApiSetupFinish = useCallback(() => {
    closeApiSetup()
    loadConnectionInfo()
    apiSetupOnboarding.reset()
  }, [closeApiSetup, loadConnectionInfo, apiSetupOnboarding])

  const handleNotificationsEnabledChange = useCallback(async (enabled: boolean) => {
    setNotificationsEnabled(enabled)
    await window.electronAPI.setNotificationsEnabled(enabled)
  }, [])

  return (
    <div className="h-full flex flex-col">
      <PanelHeader title={t('应用设置')} actions={<HeaderMenu route={routes.view.settings('app')} helpFeature="app-settings" />} />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-3xl mx-auto">
          <div className="space-y-8">
            {/* Appearance */}
            <SettingsSection title={t('外观')}>
              <SettingsCard>
                <SettingsRow label={t('主题模式')}>
                  <SettingsSegmentedControl
                    value={mode}
                    onValueChange={setMode}
                    options={[
                      { value: 'system', label: t('跟随系统'), icon: <Monitor className="w-4 h-4" /> },
                      { value: 'light', label: t('浅色'), icon: <Sun className="w-4 h-4" /> },
                      { value: 'dark', label: t('深色'), icon: <Moon className="w-4 h-4" /> },
                    ]}
                  />
                </SettingsRow>
                <SettingsRow label={t('颜色主题')}>
                  <SettingsMenuSelect
                    value={colorTheme}
                    onValueChange={setColorTheme}
                    options={[
                      { value: 'default', label: t('默认') },
                      ...presetThemes
                        .filter(theme => theme.id !== 'default')
                        .map(theme => ({
                          value: theme.id,
                          label: theme.theme.name || theme.id,
                        })),
                    ]}
                  />
                </SettingsRow>
                <SettingsRow label={t('字体')}>
                  <SettingsSegmentedControl
                    value={font}
                    onValueChange={setFont}
                    options={[
                      { value: 'inter', label: 'Inter' },
                      { value: 'system', label: t('系统字体') },
                    ]}
                  />
                </SettingsRow>
                <SettingsRow label={t('界面语言')}>
                  <SettingsSegmentedControl
                    value={locale}
                    onValueChange={setLocale}
                    options={languages.map(lang => ({
                      value: lang.code,
                      label: lang.nativeName,
                      icon: <Globe className="w-4 h-4" />,
                    }))}
                  />
                </SettingsRow>
              </SettingsCard>
            </SettingsSection>

            {/* Notifications */}
            <SettingsSection title={t('通知')}>
              <SettingsCard>
                <SettingsToggle
                  label={t('桌面通知')}
                  description={t('当 AI 完成聊天任务时发送通知')}
                  checked={notificationsEnabled}
                  onCheckedChange={handleNotificationsEnabledChange}
                />
              </SettingsCard>
            </SettingsSection>

            {/* API Connection */}
            <SettingsSection title={t('API 连接')} description={t('AI 代理连接语言模型的方式')}>
              <SettingsCard>
                <SettingsRow
                  label={t('连接类型')}
                  description={
                    authType === 'oauth_token' && hasCredential
                      ? t('Claude Pro/Max — 使用您的 Claude 订阅')
                      : authType === 'api_key' && hasCredential
                        ? t('API Key — Anthropic、OpenRouter 或兼容 API')
                        : t('未配置')
                  }
                >
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={openApiSetup}
                  >
                    {t('编辑')}
                  </Button>
                </SettingsRow>
              </SettingsCard>
            </SettingsSection>

            {/* API Setup Fullscreen Overlay — reuses the OnboardingWizard starting at the api-setup step */}
            <FullscreenOverlayBase
              isOpen={showApiSetup}
              onClose={closeApiSetup}
              className="z-splash flex flex-col bg-foreground-2"
            >
              <OnboardingWizard
                state={apiSetupOnboarding.state}
                onContinue={apiSetupOnboarding.handleContinue}
                onBack={apiSetupOnboarding.handleBack}
                onSelectApiSetupMethod={apiSetupOnboarding.handleSelectApiSetupMethod}
                onSubmitCredential={apiSetupOnboarding.handleSubmitCredential}
                onStartOAuth={apiSetupOnboarding.handleStartOAuth}
                onFinish={handleApiSetupFinish}
                existingClaudeToken={apiSetupOnboarding.existingClaudeToken}
                isClaudeCliInstalled={apiSetupOnboarding.isClaudeCliInstalled}
                onUseExistingClaudeToken={apiSetupOnboarding.handleUseExistingClaudeToken}
                isWaitingForCode={apiSetupOnboarding.isWaitingForCode}
                onSubmitAuthCode={apiSetupOnboarding.handleSubmitAuthCode}
                onCancelOAuth={apiSetupOnboarding.handleCancelOAuth}
                className="h-full"
              />
              {/* Close button — rendered AFTER the wizard so it paints above its titlebar-drag-region */}
              <div
                className="fixed top-0 right-0 h-[50px] flex items-center pr-5 [-webkit-app-region:no-drag]"
                style={{ zIndex: 'var(--z-fullscreen, 350)' }}
              >
                <button
                  onClick={closeApiSetup}
                  className="p-1.5 rounded-[6px] transition-all bg-background shadow-minimal text-muted-foreground/50 hover:text-foreground focus:outline-none focus-visible:ring-1 focus-visible:ring-ring"
                  title="Close (Esc)"
                >
                  <X className="w-3.5 h-3.5" />
                </button>
              </div>
            </FullscreenOverlayBase>

            {/* About */}
            <SettingsSection title={t('关于')}>
              <SettingsCard>
                <SettingsRow label={t('版本')}>
                  <div className="flex items-center gap-2">
                    <span className="text-muted-foreground">
                      {updateChecker.updateInfo?.currentVersion ?? t('加载中...')}
                    </span>
                    {updateChecker.updateAvailable && updateChecker.updateInfo?.latestVersion && (
                      <Button
                        variant="outline"
                        size="sm"
                        onClick={updateChecker.installUpdate}
                      >
                        {t('更新到')} {updateChecker.updateInfo.latestVersion}
                      </Button>
                    )}
                  </div>
                </SettingsRow>
                <SettingsRow label={t('检查更新')}>
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={handleCheckForUpdates}
                    disabled={isCheckingForUpdates}
                  >
                    {isCheckingForUpdates ? (
                      <>
                        <Spinner className="mr-1.5" />
                        {t('检查中...')}
                      </>
                    ) : (
                      t('立即检查')
                    )}
                  </Button>
                </SettingsRow>
                {updateChecker.isReadyToInstall && (
                  <SettingsRow label={t('安装更新')}>
                    <Button
                      size="sm"
                      onClick={updateChecker.installUpdate}
                    >
                      {t('重启并更新')}
                    </Button>
                  </SettingsRow>
                )}
              </SettingsCard>
            </SettingsSection>
          </div>
        </div>
        </ScrollArea>
      </div>
    </div>
  )
}
