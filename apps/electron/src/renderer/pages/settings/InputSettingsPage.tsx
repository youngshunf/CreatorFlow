/**
 * InputSettingsPage
 *
 * Input behavior settings that control how the chat input works.
 *
 * Settings:
 * - Auto Capitalisation (on/off)
 * - Spell Check (on/off)
 * - Send Message Key (Enter or ⌘+Enter)
 */

import { useState, useEffect, useCallback } from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { routes } from '@/lib/navigate'
import type { DetailsPageMeta } from '@/lib/navigation-registry'
import { useT } from '@/context/LocaleContext'

import {
  SettingsSection,
  SettingsCard,
  SettingsToggle,
  SettingsMenuSelectRow,
} from '@/components/settings'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'input',
}

// ============================================
// Main Component
// ============================================

export default function InputSettingsPage() {
  const t = useT()
  
  // Auto-capitalisation state
  const [autoCapitalisation, setAutoCapitalisation] = useState(true)

  // Spell check state (default off)
  const [spellCheck, setSpellCheck] = useState(false)

  // Send message key state
  const [sendMessageKey, setSendMessageKey] = useState<'enter' | 'cmd-enter'>('enter')

  // Load settings on mount
  useEffect(() => {
    const loadSettings = async () => {
      if (!window.electronAPI) return
      try {
        const [autoCapEnabled, spellCheckEnabled, sendKey] = await Promise.all([
          window.electronAPI.getAutoCapitalisation(),
          window.electronAPI.getSpellCheck(),
          window.electronAPI.getSendMessageKey(),
        ])
        setAutoCapitalisation(autoCapEnabled)
        setSpellCheck(spellCheckEnabled)
        setSendMessageKey(sendKey)
      } catch (error) {
        console.error('Failed to load input settings:', error)
      }
    }
    loadSettings()
  }, [])

  const handleAutoCapitalisationChange = useCallback(async (enabled: boolean) => {
    setAutoCapitalisation(enabled)
    await window.electronAPI.setAutoCapitalisation(enabled)
  }, [])

  const handleSpellCheckChange = useCallback(async (enabled: boolean) => {
    setSpellCheck(enabled)
    await window.electronAPI.setSpellCheck(enabled)
  }, [])

  const handleSendMessageKeyChange = useCallback(async (key: 'enter' | 'cmd-enter') => {
    setSendMessageKey(key)
    await window.electronAPI.setSendMessageKey(key)
  }, [])

  return (
    <div className="h-full flex flex-col">
      <PanelHeader title={t('输入')} actions={<HeaderMenu route={routes.view.settings('input')} />} />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-3xl mx-auto">
            <div className="space-y-8">
              {/* Typing Behavior */}
              <SettingsSection title={t('输入行为')} description={t('控制聊天输入框中文本的输入方式。')}>
                <SettingsCard>
                  <SettingsToggle
                    label={t('自动大写')}
                    description={t('输入消息时自动将首字母大写。')}
                    checked={autoCapitalisation}
                    onCheckedChange={handleAutoCapitalisationChange}
                  />
                  <SettingsToggle
                    label={t('拼写检查')}
                    description={t('输入时为拼写错误的单词添加下划线。')}
                    checked={spellCheck}
                    onCheckedChange={handleSpellCheckChange}
                  />
                </SettingsCard>
              </SettingsSection>

              {/* Send Behavior */}
              <SettingsSection title={t('消息发送')} description={t('选择发送消息的方式。')}>
                <SettingsCard>
                  <SettingsMenuSelectRow
                    label={t('发送消息快捷键')}
                    description={t('用于发送消息的键盘快捷键')}
                    value={sendMessageKey}
                    onValueChange={handleSendMessageKeyChange}
                    options={[
                      { value: 'enter', label: 'Enter', description: t('使用 Shift+Enter 换行') },
                      { value: 'cmd-enter', label: '⌘ Enter', description: t('使用 Enter 换行') },
                    ]}
                  />
                </SettingsCard>
              </SettingsSection>
            </div>
          </div>
        </ScrollArea>
      </div>
    </div>
  )
}
