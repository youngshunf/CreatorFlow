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
import { Edit2, LogOut, User, Mail, Phone, MapPin, Briefcase, Calendar } from 'lucide-react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { Button } from '@/components/ui/button'
import { Avatar, AvatarImage, AvatarFallback } from '@/components/ui/avatar'
import { useT } from '@/context/LocaleContext'
import { navigate, routes } from '@/lib/navigate'
import { Spinner } from '@creator-flow/ui'
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
