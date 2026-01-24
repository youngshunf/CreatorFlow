import { useState, useMemo } from "react"
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
} from "@/components/ui/dialog"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { AlertTriangle } from "lucide-react"
import { useRegisterModal } from "@/context/ModalContext"
import { useT } from "@/context/LocaleContext"

interface ResetConfirmationDialogProps {
  open: boolean
  onConfirm: () => void
  onCancel: () => void
}

/**
 * ResetConfirmationDialog - Destructive action confirmation with math problem
 *
 * Shows a warning about data loss and requires the user to solve a random
 * math problem to confirm the reset action.
 */
export function ResetConfirmationDialog({
  open,
  onConfirm,
  onCancel,
}: ResetConfirmationDialogProps) {
  const t = useT()
  const [answer, setAnswer] = useState("")

  // Register with modal context so X button / Cmd+W closes this dialog first
  useRegisterModal(open, onCancel)

  // Generate a random math problem when dialog opens
  const problem = useMemo(() => {
    const a = Math.floor(Math.random() * 50) + 10
    const b = Math.floor(Math.random() * 50) + 10
    return { a, b, sum: a + b }
  }, [open]) // Regenerate when dialog opens

  const isCorrect = parseInt(answer) === problem.sum

  const handleConfirm = () => {
    if (isCorrect) {
      setAnswer("")
      onConfirm()
    }
  }

  const handleCancel = () => {
    setAnswer("")
    onCancel()
  }

  return (
    <Dialog open={open} onOpenChange={(isOpen) => !isOpen && handleCancel()}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2 text-destructive">
            <AlertTriangle className="h-5 w-5" />
            {t('重置应用')}
          </DialogTitle>
          <DialogDescription className="text-left pt-2">
            {t('这将')}<strong>{t('永久删除')}</strong>{t('以下内容：')}
          </DialogDescription>
        </DialogHeader>

        <ul className="list-disc list-inside text-sm text-muted-foreground space-y-1 pl-2">
          <li>{t('所有工作区及其设置')}</li>
          <li>{t('所有凭据和 API 密钥')}</li>
          <li>{t('所有偏好设置和会话数据')}</li>
        </ul>

        <div className="bg-amber-500/10 border border-amber-500/30 rounded-md p-3 text-sm">
          <strong className="text-amber-600 dark:text-amber-400">{t('请先备份重要数据！')}</strong>
          <p className="text-muted-foreground mt-1">
            {t('此操作无法撤销。')}
          </p>
        </div>

        <div className="space-y-2 pt-2">
          <label className="text-sm font-medium">
            {t('请计算以确认：')} {problem.a} + {problem.b} =
          </label>
          <Input
            type="text"
            inputMode="numeric"
            pattern="[0-9]*"
            placeholder={t('输入答案')}
            value={answer}
            onChange={(e) => setAnswer(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter" && isCorrect) {
                handleConfirm()
              }
            }}
            className="max-w-32"
          />
        </div>

        <DialogFooter className="gap-2 sm:gap-0">
          <Button variant="outline" onClick={handleCancel}>
            {t('取消')}
          </Button>
          <Button
            variant="destructive"
            disabled={!isCorrect}
            onClick={handleConfirm}
          >
            {t('重置应用')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
