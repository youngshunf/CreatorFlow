import { useMemo } from 'react'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from "@/components/ui/dialog"
import { useRegisterModal } from "@/context/ModalContext"
import { useT } from "@/context/LocaleContext"
import { isMac } from "@/lib/platform"
import { actionsByCategory, useActionLabel, type ActionId } from "@/actions"

interface KeyboardShortcutsDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
}

interface ShortcutItem {
  keys: string[]
  description: string
}

interface ShortcutSection {
  title: string
  shortcuts: ShortcutItem[]
}

// Component-specific shortcuts that aren't in the centralized registry
// These are context-sensitive behaviors, not global actions
function useComponentSpecificSections(t: (text: string) => string): ShortcutSection[] {
  return useMemo(() => [
    {
      title: t('列表导航'),
      shortcuts: [
        { keys: ['↑', '↓'], description: t('在列表中导航') },
        { keys: ['Home'], description: t('跳到第一项') },
        { keys: ['End'], description: t('跳到最后一项') },
      ],
    },
    {
      title: t('会话列表'),
      shortcuts: [
        { keys: ['Enter'], description: t('聚焦聊天输入框') },
        { keys: ['Delete'], description: t('删除会话') },
        { keys: ['R'], description: t('重命名会话') },
        { keys: [t('右键点击')], description: t('打开上下文菜单') },
      ],
    },
    {
      title: t('智能体树'),
      shortcuts: [
        { keys: ['←'], description: t('折叠文件夹') },
        { keys: ['→'], description: t('展开文件夹') },
      ],
    },
    {
      title: t('聊天'),
      shortcuts: [
        { keys: ['Enter'], description: t('发送消息') },
        { keys: ['Shift', 'Enter'], description: t('换行') },
        { keys: ['Esc'], description: t('关闭对话框 / 取消输入焦点') },
      ],
    },
  ], [t])
}

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
    <div className="flex items-center justify-between py-1">
      <span className="text-sm">{label}</span>
      <div className="flex items-center gap-1">
        {keys.map((key, keyIndex) => (
          <Kbd key={keyIndex}>{key}</Kbd>
        ))}
      </div>
    </div>
  )
}

/**
 * Renders a section of shortcuts from the registry
 */
function RegistrySection({ category, actionIds }: { category: string; actionIds: ActionId[] }) {
  return (
    <div>
      <h3 className="text-xs font-semibold text-muted-foreground uppercase tracking-wide mb-2">
        {category}
      </h3>
      <div className="space-y-1.5">
        {actionIds.map(actionId => (
          <ActionShortcutRow key={actionId} actionId={actionId} />
        ))}
      </div>
    </div>
  )
}

/**
 * Renders a section of static shortcuts (component-specific)
 */
function StaticSection({ section }: { section: ShortcutSection }) {
  return (
    <div>
      <h3 className="text-xs font-semibold text-muted-foreground uppercase tracking-wide mb-2">
        {section.title}
      </h3>
      <div className="space-y-1.5">
        {section.shortcuts.map((shortcut, index) => (
          <div key={index} className="flex items-center justify-between py-1">
            <span className="text-sm">{shortcut.description}</span>
            <div className="flex items-center gap-1">
              {shortcut.keys.map((key, keyIndex) => (
                <Kbd key={keyIndex}>{key}</Kbd>
              ))}
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

export function KeyboardShortcutsDialog({ open, onOpenChange }: KeyboardShortcutsDialogProps) {
  const t = useT()
  const componentSpecificSections = useComponentSpecificSections(t)

  // Register with modal context so X button / Cmd+W closes this dialog first
  useRegisterModal(open, () => onOpenChange(false))

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-[500px] max-h-[80vh] overflow-y-auto" aria-describedby={undefined}>
        <DialogHeader>
          <DialogTitle>{t('键盘快捷键')}</DialogTitle>
        </DialogHeader>
        <div className="space-y-6 py-2">
          {/* Registry-driven sections */}
          {Object.entries(actionsByCategory).map(([category, actions]) => (
            <RegistrySection
              key={category}
              category={category}
              actionIds={actions.map(a => a.id as ActionId)}
            />
          ))}

          {/* Component-specific sections */}
          {componentSpecificSections.map((section) => (
            <StaticSection key={section.title} section={section} />
          ))}
        </div>
      </DialogContent>
    </Dialog>
  )
}
