import { useT } from '@/context/LocaleContext'
import { Dialog, DialogContent, DialogHeader, DialogTitle, DialogDescription, DialogFooter } from '@/components/ui/dialog'
import { Button } from '@/components/ui/button'
import { CheckCircle, XCircle, FolderOpen, Video } from 'lucide-react'

interface VideoPreviewDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  videoOutputPath: string
  contentTitle: string
  onApprove: () => void
  onReject: () => void
}

export function VideoPreviewDialog({
  open,
  onOpenChange,
  videoOutputPath,
  contentTitle,
  onApprove,
  onReject,
}: VideoPreviewDialogProps) {
  const t = useT()

  /** 在系统文件管理器中打开视频所在目录 */
  const handleOpenInExplorer = () => {
    window.electron?.shell?.showItemInFolder?.(videoOutputPath)
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2">
            <Video className="h-5 w-5 text-muted-foreground" />
            {t('视频审核')}
          </DialogTitle>
          <DialogDescription>{contentTitle}</DialogDescription>
        </DialogHeader>

        <div className="space-y-3 py-4">
          <p className="text-sm text-muted-foreground">{t('视频文件路径')}:</p>
          <div className="flex items-center gap-2 rounded-md border border-border/60 bg-muted/30 px-3 py-2">
            <code className="flex-1 truncate text-xs text-foreground">{videoOutputPath}</code>
            <Button variant="ghost" size="sm" onClick={handleOpenInExplorer} className="shrink-0">
              <FolderOpen className="mr-1 h-3.5 w-3.5" />
              {t('打开文件夹')}
            </Button>
          </div>
        </div>

        <DialogFooter className="gap-2 sm:gap-0">
          <Button variant="outline" onClick={onReject}>
            <XCircle className="mr-1.5 h-4 w-4" />
            {t('退回修改')}
          </Button>
          <Button onClick={onApprove}>
            <CheckCircle className="mr-1.5 h-4 w-4" />
            {t('审核通过')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
