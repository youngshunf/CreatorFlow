import { useState } from 'react'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
} from "@/components/ui/dialog"
import { Button } from "@/components/ui/button"
import { useRegisterModal } from "@/context/ModalContext"
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
 */
export function DeleteWorkspaceDialog({
  open,
  workspaceName,
  onConfirm,
  onCancel,
}: DeleteWorkspaceDialogProps) {
  const t = useT()
  const [mode, setMode] = useState<'delete' | 'backup'>('backup')

  // 注册到 modal context，使 X 按钮 / Cmd+W 先关闭此弹窗
  useRegisterModal(open, onCancel)

  return (
    <Dialog open={open} onOpenChange={(isOpen) => !isOpen && onCancel()}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <div className="flex justify-center mb-2">
            <img src={logoIcon} alt="" className="h-12 w-12" />
          </div>
          <DialogTitle className="text-center">
            {t('确认删除工作区')}
          </DialogTitle>
          <DialogDescription className="text-center">
            {t('确定要删除工作区')}{' "'}
            <span className="font-medium text-foreground">{workspaceName}</span>
            {'" ？'}
          </DialogDescription>
        </DialogHeader>

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

        <DialogFooter className="gap-2 sm:gap-0">
          <Button variant="outline" onClick={onCancel}>
            {t('取消')}
          </Button>
          <Button variant="destructive" onClick={() => onConfirm(mode)}>
            {mode === 'backup' ? t('备份后删除') : t('彻底删除')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
