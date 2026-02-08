import { Button } from "@/components/ui/button"
import { Spinner } from "@sprouty-ai/ui"
import { CreatorFlowSymbol } from "@/components/icons/CreatorFlowSymbol"
import { StepFormLayout } from "./primitives"
import { useT } from "@/context/LocaleContext"

interface CompletionStepProps {
  status: 'saving' | 'complete'
  spaceName?: string
  onFinish: () => void
}

/**
 * CompletionStep - Success screen after onboarding
 *
 * Shows:
 * - saving: Spinner while saving configuration
 * - complete: Success message with option to start
 */
export function CompletionStep({
  status,
  spaceName,
  onFinish
}: CompletionStepProps) {
  const t = useT()
  const isSaving = status === 'saving'

  return (
    <StepFormLayout
      iconElement={isSaving ? (
        <div className="flex size-16 items-center justify-center">
          <Spinner className="text-2xl text-foreground" />
        </div>
      ) : (
        <div className="flex size-16 items-center justify-center">
          <CreatorFlowSymbol className="size-10 text-accent" />
        </div>
      )}
      title={isSaving ? t('正在设置...') : t('一切就绪！')}
      description={
        isSaving ? (
          t('正在保存您的配置...')
        ) : (
          t('开始一个聊天，开始工作吧。')
        )
      }
      actions={
        status === 'complete' ? (
          <Button onClick={onFinish} className="w-full max-w-[320px] bg-background shadow-minimal text-foreground hover:bg-foreground/5 rounded-lg" size="lg">
            {t('开始使用')}
          </Button>
        ) : undefined
      }
    />
  )
}
