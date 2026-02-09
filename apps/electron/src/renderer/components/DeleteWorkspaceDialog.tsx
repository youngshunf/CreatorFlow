import { useState, useEffect, useCallback } from 'react'
import { createPortal } from 'react-dom'
import { Button } from "@/components/ui/button"
import { useT } from "@/context/LocaleContext"
import logoIcon from "@/assets/logo-icon.svg"

interface DeleteWorkspaceDialogProps {
  open: boolean
  workspaceName: string
  onConfirm: (mode: 'delete' | 'backup') => void
  onCancel: () => void
}

/**
 * 删除工作区确认弹窗 - 提供备份和彻底删除两个选项
 *
 * Uses a plain portal instead of Radix Dialog to avoid focus-trap conflicts
 * when rendered inside a FullscreenOverlayBase (which is itself a Radix Dialog).
 */
export function DeleteWorkspaceDialog({
  open,
  workspaceName,
  onConfirm,
  onCancel,
}: DeleteWorkspaceDialogProps) {
  const t = useT()
  const [mode, setMode] = useState<'delete' | 'backup'>('backup')

  // ESC key to close
  useEffect(() => {
    if (!open) return
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        e.stopPropagation()
        onCancel()
      }
    }
    window.addEventListener('keydown', handleKeyDown, true)
    return () => window.removeEventListener('keydown', handleKeyDown, true)
  }, [open, onCancel])

  // Reset mode when dialog opens
  useEffect(() => {
    if (open) setMode('backup')
  }, [open])

  const handleBackdropClick = useCallback((e: React.MouseEvent) => {
    if (e.target === e.currentTarget) onCancel()
  }, [onCancel])

  if (!open) return null

  return createPortal(
    <div
      className="fixed inset-0 z-[400] flex items-center justify-center"
      onClick={handleBackdropClick}
    >
      {/* Backdrop */}
      <div className="absolute inset-0 bg-black/50" />

      {/* Content */}
      <div
        role="alertdialog"
        aria-modal="true"
        aria-labelledby="delete-ws-title"
        aria-describedby="delete-ws-desc"
        className="popover-styled relative z-10 w-full max-w-md mx-4 p-6 grid gap-4 animate-in fade-in-0 zoom-in-95"
      >
        {/* Header */}
        <div className="flex flex-col space-y-1.5 text-center">
          <div className="flex justify-center mb-2">
            <img src={logoIcon} alt="" className="h-12 w-12" />
          </div>
          <h2 id="delete-ws-title" className="text-lg font-semibold leading-none tracking-tight text-center">
            {t('确认删除工作区')}
          </h2>
          <p id="delete-ws-desc" className="text-sm text-muted-foreground text-center">
            {t('确定要删除工作区')}{' "'}
            <span className="font-medium text-foreground">{workspaceName}</span>
            {'" ？'}
          </p>
        </div>

        {/* 选项 */}
        <div className="space-y-2 py-2">
          <label className="flex items-start gap-3 p-3 rounded-lg border cursor-pointer hover:bg-foreground/[0.02] transition-colors has-[:checked]:border-accent has-[:checked]:bg-accent/5">
            <input
              type="radio"
              name="deleteMode"
              value="backup"
              checked={mode === 'backup'}
              onChange={() => setMode('backup')}
              className="mt-0.5 accent-accent"
            />
            <div>
              <div className="text-sm font-medium">{t('备份后删除')}</div>
              <div className="text-xs text-muted-foreground">{t('将 .sprouty-ai 重命名为 .sprouty-ai-backup')}</div>
            </div>
          </label>
          <label className="flex items-start gap-3 p-3 rounded-lg border cursor-pointer hover:bg-foreground/[0.02] transition-colors has-[:checked]:border-destructive has-[:checked]:bg-destructive/5">
            <input
              type="radio"
              name="deleteMode"
              value="delete"
              checked={mode === 'delete'}
              onChange={() => setMode('delete')}
              className="mt-0.5 accent-destructive"
            />
            <div>
              <div className="text-sm font-medium">{t('彻底删除')}</div>
              <div className="text-xs text-muted-foreground">{t('永久删除所有工作区数据')}</div>
            </div>
          </label>
        </div>

        <p className="text-xs text-muted-foreground text-center">
          {t('此操作将删除所有会话、配置和数据。不会影响您的项目文件。')}
        </p>

        {/* Footer */}
        <div className="flex flex-col-reverse sm:flex-row sm:justify-end gap-2 sm:gap-2">
          <Button variant="outline" onClick={onCancel}>
            {t('取消')}
          </Button>
          <Button variant="destructive" onClick={() => onConfirm(mode)}>
            {mode === 'backup' ? t('备份后删除') : t('彻底删除')}
          </Button>
        </div>
      </div>
    </div>,
    document.body
  )
}
