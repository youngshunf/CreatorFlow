/**
 * ShortcutsPage
 *
 * Displays keyboard shortcuts reference from the centralized action registry.
 */

import * as React from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { SettingsSection, SettingsCard, SettingsRow } from '@/components/settings'
import { useT } from '@/context/LocaleContext'
import type { DetailsPageMeta } from '@/lib/navigation-registry'
import { isMac } from '@/lib/platform'
import { actionsByCategory, useActionLabel, type ActionId } from '@/actions'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'shortcuts',
}

interface ShortcutItem {
  keys: string[]
  descriptionKey: string
}

interface ShortcutSection {
  titleKey: string
  shortcuts: ShortcutItem[]
}

// Platform-aware modifier key symbol
const cmdKey = isMac ? '⌘' : 'Ctrl'

// Component-specific shortcuts that aren't in the centralized registry
const componentSpecificSections: ShortcutSection[] = [
  {
    titleKey: '全局',
    shortcuts: [
      { keys: [cmdKey, '1'], descriptionKey: '聚焦侧边栏' },
      { keys: [cmdKey, '2'], descriptionKey: '聚焦会话列表' },
      { keys: [cmdKey, '3'], descriptionKey: '聚焦聊天输入框' },
      { keys: [cmdKey, 'N'], descriptionKey: '新建聊天' },
      { keys: [cmdKey, 'B'], descriptionKey: '切换侧边栏' },
      { keys: [cmdKey, ','], descriptionKey: '打开设置' },
    ],
  },
  {
    titleKey: '列表导航',
    shortcuts: [
      { keys: ['↑', '↓'], descriptionKey: '在列表中导航' },
      { keys: ['Home'], descriptionKey: '跳到第一项' },
      { keys: ['End'], descriptionKey: '跳到最后一项' },
    ],
  },
  {
    titleKey: '导航',
    shortcuts: [
      { keys: ['Tab'], descriptionKey: '移动到下一个区域' },
      { keys: ['Shift', 'Tab'], descriptionKey: '切换权限模式' },
      { keys: ['←', '→'], descriptionKey: '在区域之间移动（在列表中）' },
      { keys: ['Esc'], descriptionKey: '关闭对话框 / 取消输入焦点' },
    ],
  },
  {
    titleKey: '会话列表',
    shortcuts: [
      { keys: ['Enter'], descriptionKey: '聚焦聊天输入框' },
      { keys: ['Right-click'], descriptionKey: '打开上下文菜单' },
      { keys: ['Delete'], descriptionKey: '删除会话' },
    ],
  },
  {
    titleKey: '聊天',
    shortcuts: [
      { keys: ['Enter'], descriptionKey: '发送消息' },
      { keys: ['Shift', 'Enter'], descriptionKey: '换行' },
      { keys: [cmdKey, 'Enter'], descriptionKey: '发送消息' },
      { keys: ['Esc'], descriptionKey: '关闭对话框 / 取消输入焦点' },
    ],
  },
]

function Kbd({ children }: { children: React.ReactNode }) {
  return (
    <kbd className="inline-flex items-center justify-center min-w-[20px] h-5 px-1.5 text-[11px] font-medium font-sans bg-muted border border-border rounded shadow-sm">
      {children}
    </kbd>
  )
}

/**
 * Renders a shortcut row for an action from the registry
 */
function ActionShortcutRow({ actionId }: { actionId: ActionId }) {
  const { label, hotkey } = useActionLabel(actionId)

  if (!hotkey) return null

  // Split hotkey into individual keys for display
  // Mac: symbols are concatenated (⌘⇧N) - need smart splitting
  // Windows: separated by + (Ctrl+Shift+N) - split on +
  const keys = isMac
    ? hotkey.match(/[⌘⇧⌥←→]|Tab|Esc|./g) || []
    : hotkey.split('+')

  return (
    <SettingsRow label={label}>
      <div className="flex items-center gap-1">
        {keys.map((key, keyIndex) => (
          <Kbd key={keyIndex}>{key}</Kbd>
        ))}
      </div>
    </SettingsRow>
  )
}

export default function ShortcutsPage() {
  const t = useT()
  
  return (
    <div className="h-full flex flex-col">
      <PanelHeader title={t('快捷键')} />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-3xl mx-auto space-y-8">
            {/* Registry-driven sections */}
            {Object.entries(actionsByCategory).map(([category, actions]) => (
              <SettingsSection key={category} title={category}>
                <SettingsCard>
                  {actions.map(action => (
                    <ActionShortcutRow key={action.id} actionId={action.id as ActionId} />
                  ))}
                </SettingsCard>
              </SettingsSection>
            ))}

            {/* Component-specific sections */}
            {componentSpecificSections.map((section) => (
              <SettingsSection key={section.titleKey} title={t(section.titleKey)}>
                <SettingsCard>
                  {section.shortcuts.map((shortcut, index) => (
                    <SettingsRow key={index} label={t(shortcut.descriptionKey)}>
                      <div className="flex items-center gap-1">
                        {shortcut.keys.map((key, keyIndex) => (
                          <Kbd key={keyIndex}>{key}</Kbd>
                        ))}
                      </div>
                    </SettingsRow>
                  ))}
                </SettingsCard>
              </SettingsSection>
            ))}
          </div>
        </ScrollArea>
      </div>
    </div>
  )
}
