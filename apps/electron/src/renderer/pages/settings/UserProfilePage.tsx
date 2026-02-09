/**
 * UserProfilePage
 *
 * User profile view page - displays user profile card.
 * Click edit button to navigate to edit page.
 * Includes logout functionality.
 */

import { useState, useEffect, useCallback } from 'react'
import { format } from 'date-fns'
import { zhCN } from 'date-fns/locale'
import { Edit2, LogOut, User, Mail, Phone, MapPin, Briefcase, Calendar, Key, Eye, EyeOff } from 'lucide-react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { Button } from '@/components/ui/button'
import { Avatar, AvatarImage, AvatarFallback } from '@/components/ui/avatar'
import { Dialog, DialogContent, DialogHeader, DialogTitle } from '@/components/ui/dialog'
import { Input } from '@/components/ui/input'
import { useT } from '@/context/LocaleContext'
import { navigate, routes } from '@/lib/navigate'
import { Spinner } from '@sprouty-ai/ui'
import type { DetailsPageMeta } from '@/lib/navigation-registry'
import { useAppShellContext } from '@/context/AppShellContext'

import {
  SettingsSection,
  SettingsCard,
  SettingsCardContent,
} from '@/components/settings'

import { userApi, type UserProfile } from '@/api/user'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'user-profile',
}

export default function UserProfilePage() {
  const t = useT()
  const { onLogout } = useAppShellContext()
  const [isLoading, setIsLoading] = useState(true)
  const [profile, setProfile] = useState<UserProfile | null>(null)
  const [error, setError] = useState<string | null>(null)
  
  // 修改密码弹窗状态
  const [isPasswordDialogOpen, setIsPasswordDialogOpen] = useState(false)
  const [smsCode, setSmsCode] = useState('')
  const [newPassword, setNewPassword] = useState('')
  const [confirmPassword, setConfirmPassword] = useState('')
  const [showNewPassword, setShowNewPassword] = useState(false)
  const [showConfirmPassword, setShowConfirmPassword] = useState(false)
  const [isSendingCode, setIsSendingCode] = useState(false)
  const [isChangingPassword, setIsChangingPassword] = useState(false)
  const [countdown, setCountdown] = useState(0)
  const [passwordError, setPasswordError] = useState<string | null>(null)
  const [passwordSuccess, setPasswordSuccess] = useState<string | null>(null)

  // Load user profile on mount
  useEffect(() => {
    const loadProfile = async () => {
      try {
        setIsLoading(true)
        const data = await userApi.getCurrentUser()
        setProfile(data)
      } catch (err: any) {
        console.error('Failed to load user profile:', err)
        setError(err.message || t('加载用户信息失败'))
      } finally {
        setIsLoading(false)
      }
    }
    loadProfile()
  }, [t])

  // Navigate to edit page
  const handleEdit = useCallback(() => {
    navigate(routes.view.settings('user-profile-edit'))
  }, [])

  // Handle logout - uses onLogout from AppShellContext to clear auth and return to login
  const handleLogout = useCallback(() => {
    onLogout()
  }, [onLogout])

  // 倒计时效果
  useEffect(() => {
    if (countdown > 0) {
      const timer = setTimeout(() => setCountdown(countdown - 1), 1000)
      return () => clearTimeout(timer)
    }
  }, [countdown])

  // 发送短信验证码
  const handleSendCode = useCallback(async () => {
    if (!profile?.phone) {
      setPasswordError(t('您的账号未绑定手机号'))
      return
    }
    try {
      setIsSendingCode(true)
      setPasswordError(null)
      await userApi.sendSmsCode(profile.phone)
      setCountdown(60)
    } catch (err: any) {
      setPasswordError(err.message || t('发送验证码失败'))
    } finally {
      setIsSendingCode(false)
    }
  }, [profile?.phone, t])

  // 修改密码
  const handleChangePassword = useCallback(async () => {
    setPasswordError(null)
    setPasswordSuccess(null)

    if (!smsCode) {
      setPasswordError(t('请输入验证码'))
      return
    }
    if (!newPassword) {
      setPasswordError(t('请输入新密码'))
      return
    }
    if (newPassword.length < 6) {
      setPasswordError(t('新密码长度不能少于 6 位'))
      return
    }
    if (newPassword !== confirmPassword) {
      setPasswordError(t('两次输入的密码不一致'))
      return
    }

    try {
      setIsChangingPassword(true)
      await userApi.changePasswordBySms({
        code: smsCode,
        new_password: newPassword,
        confirm_password: confirmPassword,
      })
      setPasswordSuccess(t('密码修改成功'))
      // 清空输入
      setSmsCode('')
      setNewPassword('')
      setConfirmPassword('')
      // 2秒后关闭弹窗
      setTimeout(() => {
        setIsPasswordDialogOpen(false)
        setPasswordSuccess(null)
      }, 2000)
    } catch (err: any) {
      setPasswordError(err.message || t('密码修改失败'))
    } finally {
      setIsChangingPassword(false)
    }
  }, [smsCode, newPassword, confirmPassword, t])

  // 打开密码弹窗
  const openPasswordDialog = useCallback(() => {
    setSmsCode('')
    setNewPassword('')
    setConfirmPassword('')
    setPasswordError(null)
    setPasswordSuccess(null)
    setIsPasswordDialogOpen(true)
  }, [])

  if (isLoading) {
    return (
      <div className="h-full flex items-center justify-center">
        <Spinner className="text-lg text-muted-foreground" />
      </div>
    )
  }

  if (error && !profile) {
    return (
      <div className="h-full flex items-center justify-center">
        <div className="text-center">
          <p className="text-muted-foreground">{error}</p>
        </div>
      </div>
    )
  }

  if (!profile) {
    return null
  }

  return (
    <div className="h-full flex flex-col">
      <PanelHeader 
        title={t('用户资料')} 
        actions={<HeaderMenu route={routes.view.settings('user-profile')} helpFeature="user-profile" />} 
      />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-2xl mx-auto space-y-6">
            
            {/* Profile Card */}
            <SettingsCard divided={false}>
              <div className="p-6">
                {/* Header with Avatar and Basic Info */}
                <div className="flex items-start gap-5">
                  <Avatar className="size-20 ring-2 ring-background shadow-lg">
                    <AvatarImage src={profile.avatar || undefined} alt={profile.nickname || profile.username} />
                    <AvatarFallback className="text-2xl bg-primary/10 text-primary">
                      {(profile.nickname || profile.username || 'U').charAt(0).toUpperCase()}
                    </AvatarFallback>
                  </Avatar>
                  <div className="flex-1 min-w-0">
                    <h2 className="text-xl font-semibold truncate">
                      {profile.nickname || profile.username}
                    </h2>
                    <p className="text-sm text-muted-foreground">@{profile.username}</p>
                    {profile.bio && (
                      <p className="mt-2 text-sm text-muted-foreground line-clamp-2">
                        {profile.bio}
                      </p>
                    )}
                  </div>
                  <Button variant="outline" size="sm" onClick={handleEdit}>
                    <Edit2 className="size-4 mr-1.5" />
                    {t('编辑')}
                  </Button>
                </div>

                {/* Info Grid */}
                <div className="mt-6 pt-6 border-t grid grid-cols-2 gap-4">
                  <ProfileInfoItem
                    icon={<User className="size-4" />}
                    label={t('性别')}
                    value={profile.gender ? getGenderLabel(profile.gender, t) : t('未设置')}
                  />
                  <ProfileInfoItem
                    icon={<Calendar className="size-4" />}
                    label={t('生日')}
                    value={profile.birthday ? formatBirthday(profile.birthday) : t('未设置')}
                  />
                  <ProfileInfoItem
                    icon={<Mail className="size-4" />}
                    label={t('邮箱')}
                    value={profile.email || t('未设置')}
                  />
                  <ProfileInfoItem
                    icon={<Phone className="size-4" />}
                    label={t('手机号')}
                    value={profile.phone || t('未设置')}
                  />
                  <ProfileInfoItem
                    icon={<MapPin className="size-4" />}
                    label={t('地区')}
                    value={formatLocation(profile) || t('未设置')}
                  />
                  <ProfileInfoItem
                    icon={<Briefcase className="size-4" />}
                    label={t('行业')}
                    value={profile.industry || t('未设置')}
                  />
                </div>
              </div>
            </SettingsCard>

            {/* Security Section */}
            <SettingsSection
              title={t('账号安全')}
              description={t('管理您的账号安全设置')}
            >
              <SettingsCard divided={false}>
                <SettingsCardContent className="flex items-center justify-between">
                  <div>
                    <p className="font-medium">{t('修改密码')}</p>
                    <p className="text-sm text-muted-foreground">{t('通过手机验证码修改登录密码')}</p>
                  </div>
                  <Button variant="outline" size="sm" onClick={openPasswordDialog}>
                    <Key className="size-4 mr-1.5" />
                    {t('修改密码')}
                  </Button>
                </SettingsCardContent>
              </SettingsCard>
            </SettingsSection>

            {/* Logout Section */}
            <SettingsSection
              title={t('账号操作')}
              variant="danger"
            >
              <SettingsCard divided={false}>
                <SettingsCardContent className="flex items-center justify-between">
                  <div>
                    <p className="font-medium">{t('退出登录')}</p>
                    <p className="text-sm text-muted-foreground">{t('退出当前账号')}</p>
                  </div>
                  <Button variant="destructive" size="sm" onClick={handleLogout}>
                    <LogOut className="size-4 mr-1.5" />
                    {t('退出登录')}
                  </Button>
                </SettingsCardContent>
              </SettingsCard>
            </SettingsSection>

          </div>
        </ScrollArea>
      </div>

      {/* 修改密码弹窗 */}
      <Dialog open={isPasswordDialogOpen} onOpenChange={setIsPasswordDialogOpen}>
        <DialogContent className="max-w-md" aria-describedby={undefined}>
          <DialogHeader>
            <DialogTitle>{t('修改密码')}</DialogTitle>
          </DialogHeader>

          <div className="space-y-4">
            {passwordError && (
              <div className="text-sm text-destructive bg-destructive/10 px-3 py-2 rounded">
                {passwordError}
              </div>
            )}
            {passwordSuccess && (
              <div className="text-sm text-green-600 bg-green-50 dark:bg-green-900/20 px-3 py-2 rounded">
                {passwordSuccess}
              </div>
            )}

            {/* 手机号显示 */}
            <div className="space-y-2">
              <label className="text-sm font-medium">{t('手机号')}</label>
              <div className="text-sm text-muted-foreground bg-muted px-3 py-2 rounded">
                {profile?.phone ? profile.phone.replace(/(\d{3})\d{4}(\d{4})/, '$1****$2') : t('未绑定')}
              </div>
            </div>

            {/* 验证码 */}
            <div className="space-y-2">
              <label className="text-sm font-medium">{t('验证码')}</label>
              <div className="flex gap-2">
                <Input
                  value={smsCode}
                  onChange={(e) => setSmsCode(e.target.value)}
                  placeholder={t('请输入验证码')}
                  maxLength={6}
                />
                <Button
                  variant="outline"
                  onClick={handleSendCode}
                  disabled={isSendingCode || countdown > 0 || !profile?.phone}
                  className="whitespace-nowrap"
                >
                  {countdown > 0 ? `${countdown}s` : t('获取验证码')}
                </Button>
              </div>
            </div>

            {/* 新密码 */}
            <div className="space-y-2">
              <label className="text-sm font-medium">{t('新密码')}</label>
              <div className="relative">
                <Input
                  type={showNewPassword ? 'text' : 'password'}
                  value={newPassword}
                  onChange={(e) => setNewPassword(e.target.value)}
                  placeholder={t('请输入新密码（至少 6 位）')}
                />
                <button
                  type="button"
                  className="absolute right-3 top-1/2 -translate-y-1/2 text-muted-foreground hover:text-foreground"
                  onClick={() => setShowNewPassword(!showNewPassword)}
                >
                  {showNewPassword ? <EyeOff className="size-4" /> : <Eye className="size-4" />}
                </button>
              </div>
            </div>

            {/* 确认密码 */}
            <div className="space-y-2">
              <label className="text-sm font-medium">{t('确认新密码')}</label>
              <div className="relative">
                <Input
                  type={showConfirmPassword ? 'text' : 'password'}
                  value={confirmPassword}
                  onChange={(e) => setConfirmPassword(e.target.value)}
                  placeholder={t('请再次输入新密码')}
                />
                <button
                  type="button"
                  className="absolute right-3 top-1/2 -translate-y-1/2 text-muted-foreground hover:text-foreground"
                  onClick={() => setShowConfirmPassword(!showConfirmPassword)}
                >
                  {showConfirmPassword ? <EyeOff className="size-4" /> : <Eye className="size-4" />}
                </button>
              </div>
            </div>

            <Button
              className="w-full"
              onClick={handleChangePassword}
              disabled={isChangingPassword}
            >
              {isChangingPassword && <Spinner className="size-4 mr-2" />}
              {t('确认修改')}
            </Button>
          </div>
        </DialogContent>
      </Dialog>
    </div>
  )
}

// Helper Components

interface ProfileInfoItemProps {
  icon: React.ReactNode
  label: string
  value: string
}

function ProfileInfoItem({ icon, label, value }: ProfileInfoItemProps) {
  return (
    <div className="flex items-start gap-3">
      <div className="text-muted-foreground mt-0.5">{icon}</div>
      <div className="min-w-0">
        <p className="text-xs text-muted-foreground">{label}</p>
        <p className="text-sm truncate">{value}</p>
      </div>
    </div>
  )
}

// Helper Functions

function getGenderLabel(gender: string, t: (key: string) => string): string {
  switch (gender) {
    case 'male': return t('男')
    case 'female': return t('女')
    case 'other': return t('其他')
    default: return t('未设置')
  }
}

function formatBirthday(birthday: string): string {
  try {
    const date = new Date(birthday)
    return format(date, 'yyyy年MM月dd日', { locale: zhCN })
  } catch {
    return birthday
  }
}

function formatLocation(profile: UserProfile): string | null {
  const parts = [profile.province, profile.city, profile.district].filter(Boolean)
  return parts.length > 0 ? parts.join(' ') : null
}
