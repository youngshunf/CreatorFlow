/**
 * CredentialsStep - Onboarding step wrapper for API key or OAuth flow
 *
 * Thin wrapper that composes ApiKeyInput or OAuthConnect controls
 * with StepFormLayout for the onboarding wizard context.
 */

import { ExternalLink, CheckCircle2 } from "lucide-react"
import type { ApiSetupMethod } from "./APISetupStep"
import { StepFormLayout, BackButton, ContinueButton } from "./primitives"
import {
  ApiKeyInput,
  type ApiKeyStatus,
  type ApiKeySubmitData,
  OAuthConnect,
  type OAuthStatus,
} from "../apisetup"
import { useT } from "@/context/LocaleContext"

export type CredentialStatus = ApiKeyStatus | OAuthStatus

interface CredentialsStepProps {
  apiSetupMethod: ApiSetupMethod
  status: CredentialStatus
  errorMessage?: string
  onSubmit: (data: ApiKeySubmitData) => void
  onStartOAuth?: () => void
  onBack: () => void
  // Claude OAuth specific
  existingClaudeToken?: string | null
  isClaudeCliInstalled?: boolean
  onUseExistingClaudeToken?: () => void
  // Two-step OAuth flow
  isWaitingForCode?: boolean
  onSubmitAuthCode?: (code: string) => void
  onCancelOAuth?: () => void
}

export function CredentialsStep({
  apiSetupMethod,
  status,
  errorMessage,
  onSubmit,
  onStartOAuth,
  onBack,
  existingClaudeToken,
  onUseExistingClaudeToken,
  isWaitingForCode,
  onSubmitAuthCode,
  onCancelOAuth,
}: CredentialsStepProps) {
  const t = useT()
  const isOAuth = apiSetupMethod === 'claude_oauth'

  // --- OAuth flow ---
  if (isOAuth) {
    const hasExistingToken = !!existingClaudeToken

    // Waiting for authorization code entry
    if (isWaitingForCode) {
      return (
        <StepFormLayout
          title={t('输入授权码')}
          description={t('从浏览器页面复制授权码并粘贴到下方。')}
          actions={
            <>
              <BackButton onClick={onCancelOAuth} disabled={status === 'validating'}>{t('取消')}</BackButton>
              <ContinueButton
                type="submit"
                form="auth-code-form"
                disabled={false}
                loading={status === 'validating'}
                loadingText={t('连接中...')}
              />
            </>
          }
        >
          <OAuthConnect
            status={status as OAuthStatus}
            errorMessage={errorMessage}
            existingClaudeToken={existingClaudeToken}
            isWaitingForCode={true}
            onStartOAuth={onStartOAuth!}
            onUseExistingClaudeToken={onUseExistingClaudeToken}
            onSubmitAuthCode={onSubmitAuthCode}
            onCancelOAuth={onCancelOAuth}
          />
        </StepFormLayout>
      )
    }

    // Static layout matching the API key step pattern:
    // Fixed title/description, button shows loading, error below content
    const description = hasExistingToken
      ? t('找到现有的 Claude 令牌。使用它或使用其他账户登录。')
      : t('使用您的 Claude 订阅来支持多智能体工作流。')

    return (
      <StepFormLayout
        title={t('连接 Claude 账户')}
        description={description}
        actions={
          <>
            <BackButton onClick={onBack} disabled={status === 'validating'} />
            {hasExistingToken ? (
              <ContinueButton
                onClick={onUseExistingClaudeToken}
                className="gap-2"
                loading={status === 'validating'}
                loadingText={t('连接中...')}
              >
                <CheckCircle2 className="size-4" />
                {t('使用现有令牌')}
              </ContinueButton>
            ) : (
              <ContinueButton
                onClick={onStartOAuth}
                className="gap-2"
                loading={status === 'validating'}
                loadingText={t('连接中...')}
              >
                <ExternalLink className="size-4" />
                {t('使用 Claude 登录')}
              </ContinueButton>
            )}
          </>
        }
      >
        <OAuthConnect
          status={status as OAuthStatus}
          errorMessage={errorMessage}
          existingClaudeToken={existingClaudeToken}
          isWaitingForCode={false}
          onStartOAuth={onStartOAuth!}
          onUseExistingClaudeToken={onUseExistingClaudeToken}
          onSubmitAuthCode={onSubmitAuthCode}
          onCancelOAuth={onCancelOAuth}
        />
      </StepFormLayout>
    )
  }

  // --- API Key flow ---
  return (
    <StepFormLayout
      title={t('API 配置')}
      description={t('输入您的 API 密钥。可选择性地配置 OpenRouter、Ollama 或兼容 API 的自定义端点。')}
      actions={
        <>
          <BackButton onClick={onBack} disabled={status === 'validating'} />
          <ContinueButton
            type="submit"
            form="api-key-form"
            disabled={false}
            loading={status === 'validating'}
            loadingText={t('验证中...')}
          />
        </>
      }
    >
      <ApiKeyInput
        status={status as ApiKeyStatus}
        errorMessage={errorMessage}
        onSubmit={onSubmit}
      />
    </StepFormLayout>
  )
}
