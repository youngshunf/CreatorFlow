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

interface DeleteSessionDialogProps {
  open: boolean
  sessionName: string
  onConfirm: () => void
  onCancel: () => void
}

/**
 * 删除会话确认弹窗 - 使用品牌 logo 和国际化文本
 */
export function DeleteSessionDialog({
  open,
  sessionName,
  onConfirm,
  onCancel,
}: DeleteSessionDialogProps) {
  const t = useT()

  // 注册到 modal context，使 X 按钮 / Cmd+W 先关闭此弹窗
  useRegisterModal(open, onCancel)

  return (
    <Dialog open={open} onOpenChange={(isOpen) => !isOpen && onCancel()}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <div className="flex justify-center mb-2">
            <img src={logoIcon} alt="智小芽" className="h-12 w-12" />
          </div>
          <DialogTitle className="text-center">
            {t('确认删除会话')}
          </DialogTitle>
          <DialogDescription className="text-center">
            {t('确定要删除')}{' "'}
            <span className="font-medium text-foreground">{sessionName}</span>
            {'" ？'}
          </DialogDescription>
        </DialogHeader>

        <p className="text-sm text-muted-foreground text-center">
          {t('此操作无法撤销。')}
        </p>

        <DialogFooter className="gap-2 sm:gap-0">
          <Button variant="outline" onClick={onCancel}>
            {t('取消')}
          </Button>
          <Button variant="destructive" onClick={onConfirm}>
            {t('删除')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
