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

const cmdKey = isMac ? '⌘' : 'Ctrl'

function useSections(t: (text: string) => string): ShortcutSection[] {
  return useMemo(() => [
    {
      title: t('全局'),
      shortcuts: [
        { keys: [cmdKey, '1'], description: t('聚焦侧边栏') },
        { keys: [cmdKey, '2'], description: t('聚焦会话列表') },
        { keys: [cmdKey, '3'], description: t('聚焦聊天输入框') },
        { keys: [cmdKey, 'N'], description: t('新建聊天') },
        { keys: [cmdKey, 'Shift', 'N'], description: t('新建窗口') },
        { keys: [cmdKey, '\\'], description: t('切换侧边栏') },
        { keys: [cmdKey, ','], description: t('打开设置') },
        { keys: [cmdKey, '/'], description: t('显示此对话框') },
      ],
    },
    {
      title: t('导航'),
      shortcuts: [
        { keys: ['Tab'], description: t('移动到下一个区域') },
        { keys: ['Shift', 'Tab'], description: t('移动到上一个区域') },
        { keys: ['←', '→'], description: t('在区域之间移动（在列表中）') },
        { keys: ['↑', '↓'], description: t('在列表中导航') },
        { keys: ['Home'], description: t('跳到第一项') },
        { keys: ['End'], description: t('跳到最后一项') },
        { keys: ['Esc'], description: t('关闭对话框 / 取消输入焦点') },
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
        { keys: [cmdKey, 'Enter'], description: t('发送消息') },
      ],
    },
  ], [t])
}

function Kbd({ children }: { children: React.ReactNode }) {
  return (
    <kbd className="inline-flex items-center justify-center min-w-[20px] h-5 px-1.5 text-[11px] font-medium bg-muted border border-border rounded shadow-sm">
      {children}
    </kbd>
  )
}

export function KeyboardShortcutsDialog({ open, onOpenChange }: KeyboardShortcutsDialogProps) {
  const t = useT()
  const sections = useSections(t)
  
  // Register with modal context so X button / Cmd+W closes this dialog first
  useRegisterModal(open, () => onOpenChange(false))

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-[500px] max-h-[80vh] overflow-y-auto">
        <DialogHeader>
          <DialogTitle>{t('键盘快捷键')}</DialogTitle>
        </DialogHeader>
        <div className="space-y-6 py-2">
          {sections.map((section) => (
            <div key={section.title}>
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
          ))}
        </div>
      </DialogContent>
    </Dialog>
  )
}
