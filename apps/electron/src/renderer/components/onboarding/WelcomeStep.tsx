import { CreatorFlowSymbol } from "@/components/icons/CreatorFlowSymbol"
import { StepFormLayout, ContinueButton } from "./primitives"
import { useT } from "@/context/LocaleContext"

interface WelcomeStepProps {
  onContinue: () => void
  /** Whether this is an existing user updating settings */
  isExistingUser?: boolean
}

/**
 * WelcomeStep - Initial welcome screen for onboarding
 *
 * Shows different messaging for new vs existing users:
 * - New users: Welcome to CreatorFlow
 * - Existing users: Update your API connection settings
 */
export function WelcomeStep({
  onContinue,
  isExistingUser = false
}: WelcomeStepProps) {
  const t = useT()
  return (
    <StepFormLayout
      iconElement={
        <div className="flex size-16 items-center justify-center">
          <CreatorFlowSymbol className="size-10 text-accent" />
        </div>
      }
      title={isExistingUser ? t('更新设置') : t('欢迎使用 CreatorFlow')}
      description={
        isExistingUser
          ? t('更新您的 API 连接或更改设置。')
          : t('为智能体提供应有的用户体验。连接任何服务。组织您的会话。一切所需，尽在此处！')
      }
      actions={
        <ContinueButton onClick={onContinue} className="w-full">
          {isExistingUser ? t('继续') : t('开始使用')}
        </ContinueButton>
      }
    />
  )
}
