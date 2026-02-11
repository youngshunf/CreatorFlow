import { useState, useEffect, useCallback, useRef } from 'react'
import { useT } from '@/context/LocaleContext'
import type { Platform } from '@sprouty-ai/shared/db/types'
import { PLATFORM_LIST } from '@sprouty-ai/shared/db/types'

// ============================================================
// 常量
// ============================================================

/** 平台配置 — 从共享定义导入 */
const PLATFORMS = PLATFORM_LIST

/** 登录方式 */
type AuthMethod = 'browser_profile' | 'cookie' | 'api_key'

const AUTH_METHODS: Array<{ id: AuthMethod; label: string; desc: string; recommended?: boolean }> = [
  { id: 'browser_profile', label: '浏览器登录', desc: '使用反检测浏览器自动登录，推荐方式', recommended: true },
  { id: 'cookie', label: 'Cookie 导入', desc: '从浏览器复制 Cookie 粘贴导入' },
  { id: 'api_key', label: 'API Key', desc: '使用平台开放 API 密钥' },
]

/** 步骤 */
type Step = 'platform' | 'method' | 'login' | 'done'

const LOGIN_TIMEOUT_MS = 5 * 60 * 1000 // 5 分钟

// ============================================================
// Props
// ============================================================

interface AddPlatformAccountDialogProps {
  open: boolean
  onClose: () => void
  onSuccess: () => void
  /** 预选平台（重新登录模式） */
  preselectedPlatform?: Platform
  /** 预选账号 ID（重新登录模式） */
  existingAccountId?: string
  /** Hook 方法 */
  createPlatformAccount: (data: {
    project_id: string
    platform: Platform
    auth_method: string
    auth_status: string
    auth_data: string | null
    platform_uid: string | null
    nickname: string | null
    avatar_url: string | null
    bio: string | null
    home_url: string | null
    followers: number
    following: number
    total_likes: number
    total_favorites: number
    total_comments: number
    total_posts: number
    metrics_json: string | null
    metrics_updated_at: string | null
    auth_expires_at: string | null
    last_login_at: string | null
    last_login_check: string | null
    is_primary: number
    notes: string | null
    profile_path: string | null
  }) => Promise<unknown>
  launchBrowserLogin: (platformAccountId: string, platform: Platform) => Promise<{ success: boolean; error?: string }>
  generateBrowserFingerprint: (platformAccountId: string) => Promise<unknown>
  updatePlatformAccount?: (id: string, data: Record<string, unknown>) => Promise<unknown>
  projectId: string
}

// ============================================================
// Component
// ============================================================

export function AddPlatformAccountDialog({
  open,
  onClose,
  onSuccess,
  preselectedPlatform,
  existingAccountId,
  createPlatformAccount,
  launchBrowserLogin,
  generateBrowserFingerprint,
  updatePlatformAccount,
  projectId,
}: AddPlatformAccountDialogProps) {
  const t = useT()

  // 步骤状态
  const [step, setStep] = useState<Step>('platform')
  const [selectedPlatform, setSelectedPlatform] = useState<Platform | null>(null)
  const [authMethod, setAuthMethod] = useState<AuthMethod>('browser_profile')

  // 登录状态
  const [loginStatus, setLoginStatus] = useState<'idle' | 'launching' | 'waiting' | 'success' | 'failed'>('idle')
  const [loginError, setLoginError] = useState<string | null>(null)
  const [countdown, setCountdown] = useState(300) // 5 分钟倒计时（秒）
  const countdownRef = useRef<ReturnType<typeof setInterval> | null>(null)

  // Cookie / API Key 输入
  const [cookieText, setCookieText] = useState('')
  const [apiKey, setApiKey] = useState('')

  // 账号 ID（新建或已有）
  const [accountId, setAccountId] = useState<string | null>(null)

  // 指纹信息
  const [fingerprintGenerated, setFingerprintGenerated] = useState(false)

  /** 重置状态 */
  const resetState = useCallback(() => {
    setStep(preselectedPlatform ? 'method' : 'platform')
    setSelectedPlatform(preselectedPlatform ?? null)
    setAuthMethod('browser_profile')
    setLoginStatus('idle')
    setLoginError(null)
    setCountdown(300)
    setCookieText('')
    setApiKey('')
    setAccountId(existingAccountId ?? null)
    setFingerprintGenerated(false)
    if (countdownRef.current) {
      clearInterval(countdownRef.current)
      countdownRef.current = null
    }
  }, [preselectedPlatform, existingAccountId])

  useEffect(() => {
    if (open) resetState()
    return () => {
      if (countdownRef.current) clearInterval(countdownRef.current)
    }
  }, [open, resetState])

  /** 开始倒计时 */
  const startCountdown = useCallback(() => {
    setCountdown(300)
    if (countdownRef.current) clearInterval(countdownRef.current)
    countdownRef.current = setInterval(() => {
      setCountdown((prev) => {
        if (prev <= 1) {
          if (countdownRef.current) clearInterval(countdownRef.current)
          return 0
        }
        return prev - 1
      })
    }, 1000)
  }, [])

  /** 格式化倒计时 */
  const formatCountdown = (seconds: number) => {
    const m = Math.floor(seconds / 60)
    const s = seconds % 60
    return `${m}:${s.toString().padStart(2, '0')}`
  }

  /** 创建账号记录（如果是新建模式） */
  const ensureAccount = useCallback(async (platform: Platform): Promise<string> => {
    if (existingAccountId) return existingAccountId

    const id = crypto.randomUUID()
    await createPlatformAccount({
      project_id: projectId,
      platform,
      auth_method: authMethod,
      auth_status: 'not_logged_in',
      auth_data: null,
      platform_uid: null,
      nickname: null,
      avatar_url: null,
      bio: null,
      home_url: null,
      followers: 0,
      following: 0,
      total_likes: 0,
      total_favorites: 0,
      total_comments: 0,
      total_posts: 0,
      metrics_json: null,
      metrics_updated_at: null,
      auth_expires_at: null,
      last_login_at: null,
      last_login_check: null,
      is_primary: 0,
      notes: null,
      profile_path: null,
    })
    setAccountId(id)
    return id
  }, [existingAccountId, createPlatformAccount, projectId, authMethod])

  /** 执行浏览器登录 */
  const handleBrowserLogin = useCallback(async () => {
    if (!selectedPlatform) return

    setLoginStatus('launching')
    setLoginError(null)

    try {
      const id = await ensureAccount(selectedPlatform)

      setLoginStatus('waiting')
      startCountdown()

      const result = await launchBrowserLogin(id, selectedPlatform)

      if (countdownRef.current) clearInterval(countdownRef.current)

      if (result.success) {
        setLoginStatus('success')
        // 自动生成指纹
        try {
          await generateBrowserFingerprint(id)
          setFingerprintGenerated(true)
        } catch {
          // 指纹生成失败不影响登录成功
        }
        setStep('done')
      } else {
        setLoginStatus('failed')
        setLoginError(result.error || t('登录失败，请重试'))
      }
    } catch (err) {
      if (countdownRef.current) clearInterval(countdownRef.current)
      setLoginStatus('failed')
      setLoginError(err instanceof Error ? err.message : t('启动浏览器失败'))
    }
  }, [selectedPlatform, ensureAccount, launchBrowserLogin, generateBrowserFingerprint, startCountdown, t])

  /** 提交 Cookie */
  const handleCookieSubmit = useCallback(async () => {
    if (!selectedPlatform || !cookieText.trim()) return

    setLoginStatus('launching')
    setLoginError(null)

    try {
      const id = await ensureAccount(selectedPlatform)

      if (updatePlatformAccount) {
        await updatePlatformAccount(id, {
          auth_status: 'logged_in',
          auth_method: 'cookie',
          auth_data: cookieText.trim(),
          last_login_at: new Date().toISOString(),
        })
      }

      setLoginStatus('success')
      setStep('done')
    } catch (err) {
      setLoginStatus('failed')
      setLoginError(err instanceof Error ? err.message : t('保存 Cookie 失败'))
    }
  }, [selectedPlatform, cookieText, ensureAccount, updatePlatformAccount, t])

  /** 提交 API Key */
  const handleApiKeySubmit = useCallback(async () => {
    if (!selectedPlatform || !apiKey.trim()) return

    setLoginStatus('launching')
    setLoginError(null)

    try {
      const id = await ensureAccount(selectedPlatform)

      if (updatePlatformAccount) {
        await updatePlatformAccount(id, {
          auth_status: 'logged_in',
          auth_method: 'api_key',
          auth_data: apiKey.trim(),
          last_login_at: new Date().toISOString(),
        })
      }

      setLoginStatus('success')
      setStep('done')
    } catch (err) {
      setLoginStatus('failed')
      setLoginError(err instanceof Error ? err.message : t('保存 API Key 失败'))
    }
  }, [selectedPlatform, apiKey, ensureAccount, updatePlatformAccount, t])

  if (!open) return null

  const isReloginMode = !!existingAccountId

  return (
    <div className="fixed inset-0 z-modal flex items-center justify-center">
      {/* 遮罩 */}
      <div className="absolute inset-0 bg-black/50" onClick={onClose} />

      {/* 对话框 */}
      <div className="relative w-full max-w-lg rounded-lg border border-border bg-background shadow-lg mx-4">
        {/* 头部 */}
        <div className="px-6 py-4 border-b border-border/40">
          <h2 className="text-base font-semibold text-foreground">
            {isReloginMode ? t('重新登录') : t('添加平台账号')}
          </h2>
          <p className="text-xs text-muted-foreground mt-1">
            {step === 'platform' && t('选择要添加的社交媒体平台')}
            {step === 'method' && t('选择登录方式')}
            {step === 'login' && t('完成账号登录')}
            {step === 'done' && t('账号添加成功')}
          </p>
          {/* 步骤指示器 */}
          <div className="flex items-center gap-1.5 mt-3">
            {(['platform', 'method', 'login', 'done'] as Step[]).map((s, i) => (
              <div
                key={s}
                className={`h-1 flex-1 rounded-full transition-colors ${
                  i <= ['platform', 'method', 'login', 'done'].indexOf(step)
                    ? 'bg-foreground'
                    : 'bg-border/60'
                }`}
              />
            ))}
          </div>
        </div>

        {/* 内容 */}
        <div className="px-6 py-5 max-h-[60vh] overflow-auto">
          {/* 步骤 1：选择平台 */}
          {step === 'platform' && (
            <div className="grid grid-cols-2 gap-3">
              {PLATFORMS.map((p) => (
                <button
                  key={p.id}
                  type="button"
                  onClick={() => {
                    setSelectedPlatform(p.id)
                    setStep('method')
                  }}
                  className={`flex flex-col items-start gap-1 rounded-lg border p-3 text-left transition-colors hover:bg-muted/30 ${
                    selectedPlatform === p.id
                      ? `${p.color} bg-muted/20`
                      : 'border-border/60'
                  }`}
                >
                  <span className={`text-sm font-semibold ${p.color.split(' ')[0]}`}>{t(p.label)}</span>
                  <span className="text-[11px] text-muted-foreground">{t(p.desc)}</span>
                </button>
              ))}
            </div>
          )}

          {/* 步骤 2：选择登录方式 */}
          {step === 'method' && (
            <div className="space-y-3">
              {AUTH_METHODS.map((m) => (
                <button
                  key={m.id}
                  type="button"
                  onClick={() => setAuthMethod(m.id)}
                  className={`w-full flex items-start gap-3 rounded-lg border p-3 text-left transition-colors hover:bg-muted/30 ${
                    authMethod === m.id
                      ? 'border-foreground/40 bg-muted/20'
                      : 'border-border/60'
                  }`}
                >
                  {/* 单选圆点 */}
                  <div className={`mt-0.5 w-4 h-4 rounded-full border-2 flex items-center justify-center shrink-0 ${
                    authMethod === m.id ? 'border-foreground' : 'border-border'
                  }`}>
                    {authMethod === m.id && <div className="w-2 h-2 rounded-full bg-foreground" />}
                  </div>
                  <div>
                    <div className="flex items-center gap-2">
                      <span className="text-sm font-medium text-foreground">{t(m.label)}</span>
                      {m.recommended && (
                        <span className="inline-flex rounded-full bg-green-500/15 px-1.5 py-0.5 text-[10px] font-medium text-green-600 dark:text-green-400">
                          {t('推荐')}
                        </span>
                      )}
                    </div>
                    <span className="text-xs text-muted-foreground">{t(m.desc)}</span>
                  </div>
                </button>
              ))}
            </div>
          )}

          {/* 步骤 3：执行登录 */}
          {step === 'login' && authMethod === 'browser_profile' && (
            <div className="space-y-4">
              {loginStatus === 'idle' && (
                <div className="text-center space-y-3 py-4">
                  <p className="text-sm text-foreground">
                    {t('点击下方按钮启动反检测浏览器')}
                  </p>
                  <p className="text-xs text-muted-foreground">
                    {t('浏览器窗口将自动打开平台登录页，请在 5 分钟内完成登录')}
                  </p>
                  <button
                    type="button"
                    onClick={handleBrowserLogin}
                    className="inline-flex items-center gap-2 rounded-md bg-foreground px-4 py-2 text-sm font-medium text-background hover:bg-foreground/90 transition-colors"
                  >
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="M12 21a9.004 9.004 0 0 0 8.716-6.747M12 21a9.004 9.004 0 0 1-8.716-6.747M12 21c2.485 0 4.5-4.03 4.5-9S14.485 3 12 3m0 18c-2.485 0-4.5-4.03-4.5-9S9.515 3 12 3m0 0a8.997 8.997 0 0 1 7.843 4.582M12 3a8.997 8.997 0 0 0-7.843 4.582m15.686 0A11.953 11.953 0 0 1 12 10.5c-2.998 0-5.74-1.1-7.843-2.918m15.686 0A8.959 8.959 0 0 1 21 12c0 .778-.099 1.533-.284 2.253m0 0A17.919 17.919 0 0 1 12 16.5c-3.162 0-6.133-.815-8.716-2.247m0 0A9.015 9.015 0 0 1 3 12c0-1.605.42-3.113 1.157-4.418" />
                    </svg>
                    {t('启动浏览器登录')}
                  </button>
                </div>
              )}

              {loginStatus === 'launching' && (
                <div className="text-center space-y-3 py-8">
                  <div className="mx-auto w-8 h-8 border-2 border-foreground/20 border-t-foreground rounded-full animate-spin" />
                  <p className="text-sm text-foreground">{t('正在启动浏览器...')}</p>
                </div>
              )}

              {loginStatus === 'waiting' && (
                <div className="text-center space-y-3 py-4">
                  <div className="mx-auto w-8 h-8 border-2 border-foreground/20 border-t-foreground rounded-full animate-spin" />
                  <p className="text-sm text-foreground">{t('请在弹出的浏览器窗口中完成登录')}</p>
                  <p className="text-xs text-muted-foreground">
                    {t('登录成功后将自动检测，剩余时间：')}{formatCountdown(countdown)}
                  </p>
                  {/* 进度条 */}
                  <div className="w-full h-1.5 bg-border/40 rounded-full overflow-hidden">
                    <div
                      className="h-full bg-foreground rounded-full transition-all duration-1000"
                      style={{ width: `${(countdown / 300) * 100}%` }}
                    />
                  </div>
                </div>
              )}

              {loginStatus === 'failed' && (
                <div className="text-center space-y-3 py-4">
                  <div className="mx-auto w-10 h-10 rounded-full bg-red-500/10 flex items-center justify-center">
                    <svg className="w-5 h-5 text-red-500" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="M6 18 18 6M6 6l12 12" />
                    </svg>
                  </div>
                  <p className="text-sm text-foreground">{t('登录失败')}</p>
                  {loginError && (
                    <p className="text-xs text-red-500">{loginError}</p>
                  )}
                  <button
                    type="button"
                    onClick={handleBrowserLogin}
                    className="inline-flex items-center gap-1.5 rounded-md border border-border/60 px-3 py-1.5 text-xs text-foreground hover:bg-muted/40 transition-colors"
                  >
                    {t('重试')}
                  </button>
                </div>
              )}
            </div>
          )}

          {step === 'login' && authMethod === 'cookie' && (
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-foreground mb-1.5">
                  {t('粘贴 Cookie')}
                </label>
                <textarea
                  value={cookieText}
                  onChange={(e) => setCookieText(e.target.value)}
                  placeholder={t('从浏览器开发者工具复制 Cookie 字符串...')}
                  rows={6}
                  className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-xs text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-foreground/30 font-mono resize-none"
                />
                <p className="text-[11px] text-muted-foreground mt-1">
                  {t('打开浏览器开发者工具 (F12) > Application > Cookies，复制所有 Cookie')}
                </p>
              </div>

              {loginError && (
                <p className="text-xs text-red-500">{loginError}</p>
              )}

              <button
                type="button"
                onClick={handleCookieSubmit}
                disabled={!cookieText.trim() || loginStatus === 'launching'}
                className="w-full rounded-md bg-foreground px-4 py-2 text-sm font-medium text-background hover:bg-foreground/90 transition-colors disabled:opacity-50"
              >
                {loginStatus === 'launching' ? t('保存中...') : t('保存 Cookie')}
              </button>
            </div>
          )}

          {step === 'login' && authMethod === 'api_key' && (
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-foreground mb-1.5">
                  {t('API Key')}
                </label>
                <input
                  type="password"
                  value={apiKey}
                  onChange={(e) => setApiKey(e.target.value)}
                  placeholder={t('输入平台 API Key...')}
                  className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-xs text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-foreground/30 font-mono"
                />
                <p className="text-[11px] text-muted-foreground mt-1">
                  {t('在平台开放平台获取 API Key，密钥将加密存储')}
                </p>
              </div>

              {loginError && (
                <p className="text-xs text-red-500">{loginError}</p>
              )}

              <button
                type="button"
                onClick={handleApiKeySubmit}
                disabled={!apiKey.trim() || loginStatus === 'launching'}
                className="w-full rounded-md bg-foreground px-4 py-2 text-sm font-medium text-background hover:bg-foreground/90 transition-colors disabled:opacity-50"
              >
                {loginStatus === 'launching' ? t('保存中...') : t('保存 API Key')}
              </button>
            </div>
          )}

          {/* 步骤 4：完成 */}
          {step === 'done' && (
            <div className="text-center space-y-4 py-4">
              <div className="mx-auto w-12 h-12 rounded-full bg-green-500/10 flex items-center justify-center">
                <svg className="w-6 h-6 text-green-500" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="m4.5 12.75 6 6 9-13.5" />
                </svg>
              </div>
              <div>
                <p className="text-sm font-medium text-foreground">
                  {isReloginMode ? t('重新登录成功') : t('账号添加成功')}
                </p>
                <p className="text-xs text-muted-foreground mt-1">
                  {selectedPlatform && t(PLATFORMS.find(p => p.id === selectedPlatform)?.label || '')}
                  {' '}{t('账号已就绪')}
                </p>
              </div>
              {fingerprintGenerated && (
                <div className="rounded-md border border-border/40 bg-muted/20 p-3 text-left">
                  <p className="text-[11px] font-medium text-foreground mb-1">{t('浏览器指纹已生成')}</p>
                  <p className="text-[11px] text-muted-foreground">
                    {t('指纹配置已保存，后续自动发布和数据采集将使用此指纹')}
                  </p>
                </div>
              )}
            </div>
          )}
        </div>

        {/* 底部按钮 */}
        <div className="flex items-center justify-between px-6 py-4 border-t border-border/40">
          <div>
            {step === 'method' && (
              <button
                type="button"
                onClick={() => {
                  if (preselectedPlatform) {
                    onClose()
                  } else {
                    setStep('platform')
                  }
                }}
                className="rounded-md border border-border/60 bg-background px-3 py-1.5 text-xs text-foreground hover:bg-muted/40 transition-colors"
              >
                {t('上一步')}
              </button>
            )}
            {step === 'login' && loginStatus === 'idle' && (
              <button
                type="button"
                onClick={() => setStep('method')}
                className="rounded-md border border-border/60 bg-background px-3 py-1.5 text-xs text-foreground hover:bg-muted/40 transition-colors"
              >
                {t('上一步')}
              </button>
            )}
          </div>
          <div className="flex items-center gap-2">
            {step !== 'done' && (
              <button
                type="button"
                onClick={onClose}
                className="rounded-md border border-border/60 bg-background px-3 py-1.5 text-xs text-foreground hover:bg-muted/40 transition-colors"
              >
                {t('取消')}
              </button>
            )}
            {step === 'method' && (
              <button
                type="button"
                onClick={() => setStep('login')}
                className="rounded-md bg-foreground px-3 py-1.5 text-xs font-medium text-background hover:bg-foreground/90 transition-colors"
              >
                {t('下一步')}
              </button>
            )}
            {step === 'done' && (
              <button
                type="button"
                onClick={() => {
                  onSuccess()
                  onClose()
                }}
                className="rounded-md bg-foreground px-3 py-1.5 text-xs font-medium text-background hover:bg-foreground/90 transition-colors"
              >
                {t('完成')}
              </button>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}
