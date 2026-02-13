import React, { useState, useEffect } from 'react'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { authApi, type LoginResult } from '@/api/auth'
import { getUserLlmConfig } from '@/api/models'
import { cn } from '@/lib/utils'
import { enableCloudMode, type CloudConfig, getCloudApiUrl, getCurrentEnv, isDebugMode } from '@sprouty-ai/shared/cloud'
import { toast } from 'sonner'

interface LoginScreenProps {
  onLoginSuccess: () => void
}

export function LoginScreen({ onLoginSuccess }: LoginScreenProps) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  
  // Phone Login State
  const [phone, setPhone] = useState('')
  const [code, setCode] = useState('')
  const [countdown, setCountdown] = useState(0)
  const [sendingCode, setSendingCode] = useState(false)
  
  // Password Login State
  const [username, setUsername] = useState('')
  const [password, setPassword] = useState('')

  useEffect(() => {
    let timer: NodeJS.Timeout
    if (countdown > 0) {
      timer = setInterval(() => {
        setCountdown((prev) => prev - 1)
      }, 1000)
    }
    return () => clearInterval(timer)
  }, [countdown])

  const handleSendCode = async () => {
    if (!phone) {
      setError('请输入手机号')
      return
    }
    setSendingCode(true)
    setError(null)
    try {
      await authApi.sendSmsCode(phone)
      setCountdown(60)
      toast.success('验证码已发送')
    } catch (err: any) {
      const msg = err.message || '发送验证码失败'
      setError(msg)
      toast.error(msg)
    } finally {
      setSendingCode(false)
    }
  }

  const handlePhoneLogin = async (e: React.FormEvent) => {
    e.preventDefault()
    if (!phone || !code) {
      setError('请输入手机号和验证码')
      return
    }
    
    setLoading(true)
    setError(null)
    try {
      const res = await authApi.loginByPhone(phone, code)
      handleSuccess(res)
    } catch (err: any) {
      const msg = err.message || '登录失败'
      setError(msg)
      toast.error(msg)
    } finally {
      setLoading(false)
    }
  }

  const handlePasswordLogin = async (e: React.FormEvent) => {
    e.preventDefault()
    if (!username || !password) {
      setError('请输入手机号和密码')
      return
    }

    setLoading(true)
    setError(null)
    try {
      const res = await authApi.loginByPassword(username, password)
      handleSuccess(res)
    } catch (err: any) {
      const msg = err.message || '登录失败'
      setError(msg)
      toast.error(msg)
    } finally {
      setLoading(false)
    }
  }

  const handleSuccess = async (res: LoginResult) => {
    // Save tokens to localStorage
    localStorage.setItem('access_token', res.access_token)
    if (res.llm_token) {
      localStorage.setItem('llm_token', res.llm_token)
    }

    // Enable cloud mode for LLM calls
    if (res.llm_token) {
      // Get API URL from environment configuration
      const apiBaseUrl = getCloudApiUrl()

      if (isDebugMode()) {
        console.log(`[Login] Environment: ${getCurrentEnv()}, API URL: ${apiBaseUrl}`)
      }

      const cloudConfig: CloudConfig = {
        apiBaseUrl,
        accessToken: res.access_token,
        llmToken: res.llm_token,
        expiresAt: res.access_token_expire_time ? new Date(res.access_token_expire_time).getTime() : undefined,
      }
      // Save to localStorage for renderer access
      enableCloudMode(cloudConfig)

      // Notify main process to use cloud gateway for SDK
      // Build full gateway URL: apiBaseUrl + /llm/proxy
      const gatewayUrl = `${apiBaseUrl.replace(/\/$/, '')}/llm/proxy`
      try {
        // 使用新版 setCloudAuth 传递完整令牌（含 refreshToken）
        await window.electronAPI.setCloudAuth({
          accessToken: res.access_token,
          refreshToken: res.refresh_token,
          llmToken: res.llm_token,
          gatewayUrl,
          expiresAt: res.access_token_expire_time ? new Date(res.access_token_expire_time).getTime() : undefined,
          refreshExpiresAt: res.refresh_token_expire_time ? new Date(res.refresh_token_expire_time).getTime() : undefined,
        })
        if (isDebugMode()) {
          console.log(`[Login] Cloud auth set in main process: ${gatewayUrl}`)
        }

        // 获取用户的 LLM 配置（包括默认模型）
        try {
          const userConfig = await getUserLlmConfig()
          if (userConfig.default_anthropic_model || userConfig.default_openai_model) {
            await window.electronAPI.updateCloudConnectionModels(
              userConfig.default_anthropic_model,
              userConfig.default_openai_model
            )
            if (isDebugMode()) {
              console.log('[Login] Updated cloud connection models:', userConfig)
            }
          }
        } catch (err) {
          // 如果获取用户配置失败，不影响登录流程
          console.warn('[Login] Failed to fetch user LLM config:', err)
        }
      } catch (err) {
        console.error('[Login] Failed to set cloud auth in main process:', err)
      }
    }

    onLoginSuccess()
  }

  return (
    <div className="flex min-h-screen flex-col items-center justify-center bg-background p-4">
      {/* Draggable title bar region */}
      <div className="titlebar-drag-region fixed top-0 left-0 right-0 h-[50px] z-50" />
      
      <div className="w-full max-w-md space-y-8 z-10">
        <div className="text-center">
          <h2 className="text-2xl font-bold tracking-tight text-foreground">
            登录智小芽
          </h2>
          <p className="mt-2 text-sm text-muted-foreground">
            使用云端服务账号登录
          </p>
        </div>

        <Tabs defaultValue="phone" className="w-full">
          <TabsList className="grid w-full grid-cols-2">
            <TabsTrigger value="phone">手机号登录</TabsTrigger>
            <TabsTrigger value="password">账号密码</TabsTrigger>
          </TabsList>
          
          <TabsContent value="phone">
            <form onSubmit={handlePhoneLogin} className="space-y-4 pt-4">
              <div className="space-y-2">
                <Label htmlFor="phone">手机号</Label>
                <Input
                  id="phone"
                  placeholder="请输入手机号"
                  value={phone}
                  onChange={(e) => setPhone(e.target.value)}
                />
              </div>
              <div className="space-y-2">
                <Label htmlFor="code">验证码</Label>
                <div className="flex gap-2">
                  <Input
                    id="code"
                    placeholder="6位验证码"
                    value={code}
                    onChange={(e) => setCode(e.target.value)}
                  />
                  <Button 
                    type="button" 
                    variant="outline" 
                    disabled={countdown > 0 || !phone || sendingCode}
                    onClick={handleSendCode}
                    className="w-32"
                  >
                    {sendingCode ? '发送中...' : countdown > 0 ? `${countdown}s` : '获取验证码'}
                  </Button>
                </div>
              </div>
              {error && <p className="text-sm text-red-500">{error}</p>}
              <Button type="submit" className="w-full" disabled={loading}>
                {loading ? '登录中...' : '登录 / 注册'}
              </Button>
            </form>
          </TabsContent>
          
          <TabsContent value="password">
            <form onSubmit={handlePasswordLogin} className="space-y-4 pt-4">
              <div className="space-y-2">
                <Label htmlFor="username">手机号</Label>
                <Input
                  id="username"
                  placeholder="请输入手机号"
                  value={username}
                  onChange={(e) => setUsername(e.target.value)}
                />
              </div>
              <div className="space-y-2">
                <Label htmlFor="password">密码</Label>
                <Input
                  id="password"
                  type="password"
                  placeholder="请输入密码"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
              </div>
              {error && <p className="text-sm text-red-500">{error}</p>}
              <Button type="submit" className="w-full" disabled={loading}>
                {loading ? '登录中...' : '登录'}
              </Button>
            </form>
          </TabsContent>
        </Tabs>
      </div>
    </div>
  )
}
