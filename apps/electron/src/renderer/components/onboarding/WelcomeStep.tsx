import { SproutySymbol } from "@/components/icons/SproutySymbol"
import { StepFormLayout, ContinueButton } from "./primitives"
import { useT } from "@/context/LocaleContext"

interface WelcomeStepProps {
  onContinue: () => void
  /** Whether this is an existing user updating settings */
  isExistingUser?: boolean
  /** Whether the app is loading (e.g., checking Git Bash on Windows) */
  isLoading?: boolean
}

/**
 * WelcomeStep - Initial welcome screen for onboarding
 *
 * Shows different messaging for new vs existing users:
 * - New users: Welcome to Sprouty AI
 * - Existing users: Update your API connection settings
 */
export function WelcomeStep({
  onContinue,
  isExistingUser = false,
  isLoading = false
}: WelcomeStepProps) {
  const t = useT()
  return (
    <StepFormLayout
      iconElement={
        <div className="flex size-16 items-center justify-center">
          <SproutySymbol className="size-10 text-accent" />
        </div>
      }
      title={isExistingUser ? t('更新设置') : t('欢迎使用智小芽')}
      description={
        isExistingUser
          ? t('更新您的 API 连接或更改设置。')
          : t('为智能体提供应有的用户体验。连接任何服务。组织您的会话。一切所需，尽在此处！')
      }
      actions={
        <ContinueButton onClick={onContinue} className="w-full" loading={isLoading} loadingText={t('检查中...')}>
          {isExistingUser ? t('继续') : t('开始使用')}
        </ContinueButton>
      }
    />
  )
}
