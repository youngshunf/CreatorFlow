/**
 * AppSettingsPage
 *
 * Global app-level settings that apply across all workspaces.
 *
 * Settings:
 * - Appearance (Theme, Font)
 * - Notifications
 * - About (version, updates)
 *
 * Note: AI settings (connections, model, thinking) have been moved to AiSettingsPage.
 * Note: Appearance settings (theme, font) have been moved to AppearanceSettingsPage.
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
import { Spinner, FullscreenOverlayBase } from '@sprouty-ai/ui'
import { useSetAtom } from 'jotai'
import { fullscreenOverlayOpenAtom } from '@/atoms/overlay'
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

  // Notifications state
  const [notificationsEnabled, setNotificationsEnabled] = useState(true)

  // Power state
  const [keepAwakeEnabled, setKeepAwakeEnabled] = useState(false)

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

  // Load settings on mount
  const loadSettings = useCallback(async () => {
    if (!window.electronAPI) return
    try {
      const [notificationsOn, keepAwakeOn] = await Promise.all([
        window.electronAPI.getNotificationsEnabled(),
        window.electronAPI.getKeepAwakeWhileRunning(),
      ])
      setNotificationsEnabled(notificationsOn)
      setKeepAwakeEnabled(keepAwakeOn)
    } catch (error) {
      console.error('Failed to load settings:', error)
    }
  }, [])

  useEffect(() => {
    loadSettings()
  }, [])

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

  const handleNotificationsEnabledChange = useCallback(async (enabled: boolean) => {
    setNotificationsEnabled(enabled)
    await window.electronAPI.setNotificationsEnabled(enabled)
  }, [])

  const handleKeepAwakeEnabledChange = useCallback(async (enabled: boolean) => {
    setKeepAwakeEnabled(enabled)
    await window.electronAPI.setKeepAwakeWhileRunning(enabled)
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

            {/* Power */}
            <SettingsSection title={t('电源')}>
              <SettingsCard>
                <SettingsToggle
                  label={t('保持屏幕常亮')}
                  description={t('会话运行时防止屏幕关闭')}
                  checked={keepAwakeEnabled}
                  onCheckedChange={handleKeepAwakeEnabledChange}
                />
              </SettingsCard>
            </SettingsSection>

            {/* About */}
            <SettingsSection title={t('关于')}>
              <SettingsCard>
                <SettingsRow label={t('版本')}>
                  <div className="flex items-center gap-2">
                    <span className="text-muted-foreground">
                      {updateChecker.updateInfo?.currentVersion ?? t('加载中...')}
                    </span>
                    {/* 下载进度指示器 */}
                    {updateChecker.isDownloading && updateChecker.updateInfo?.latestVersion && (
                      <div className="flex items-center gap-2 text-muted-foreground text-sm">
                        <Spinner className="w-3 h-3" />
                        {updateChecker.isIndeterminate ? (
                          <span>{t('正在下载')} v{updateChecker.updateInfo.latestVersion}...</span>
                        ) : (
                          <span>{t('正在下载')} v{updateChecker.updateInfo.latestVersion} ({updateChecker.downloadProgress}%)</span>
                        )}
                      </div>
                    )}
                    {updateChecker.updateAvailable && !updateChecker.isDownloading && updateChecker.updateInfo?.latestVersion && (
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
                {updateChecker.isReadyToInstall && updateChecker.updateInfo?.latestVersion && (
                  <SettingsRow label={t('安装更新')}>
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={handleCheckForUpdates}
                      disabled={isCheckingForUpdates}
                    >
                      {t('重启并更新到')} v{updateChecker.updateInfo.latestVersion}
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
