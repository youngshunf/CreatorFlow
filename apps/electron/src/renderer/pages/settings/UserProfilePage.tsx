/**
 * UserProfilePage
 *
 * User profile settings page - displays and edits user information
 * including avatar, nickname, gender, birthday, location, industry, and bio.
 */

import { useState, useEffect, useCallback } from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { useT } from '@/context/LocaleContext'
import { routes } from '@/lib/navigate'
import { Spinner } from '@creator-flow/ui'
import type { DetailsPageMeta } from '@/lib/navigation-registry'

import {
  SettingsSection,
  SettingsCard,
  SettingsInput,
  SettingsSelect,
  SettingsTextarea,
} from '@/components/settings'

import { userApi, type UserProfile, type UpdateUserProfileParams } from '@/api/user'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'user-profile',
}

export default function UserProfilePage() {
  const t = useT()
  const [isLoading, setIsLoading] = useState(true)
  const [isSaving, setIsSaving] = useState(false)
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
        setError(err.message || '加载用户信息失败')
      } finally {
        setIsLoading(false)
      }
    }
    loadProfile()
  }, [])

  // Update field and save to server
  const updateField = useCallback(async <K extends keyof UpdateUserProfileParams>(
    field: K,
    value: UpdateUserProfileParams[K]
  ) => {
    if (!profile) return

    // Optimistic update
    setProfile(prev => prev ? ({ ...prev, [field]: value }) : null)

    try {
      setIsSaving(true)
      await userApi.updateProfile({ [field]: value })
      setError(null)
    } catch (err: any) {
      console.error(`Failed to update ${field}:`, err)
      setError(err.message || '更新失败')
      // Reload profile on error
      const data = await userApi.getCurrentUser()
      setProfile(data)
    } finally {
      setIsSaving(false)
    }
  }, [profile])

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
          <div className="px-5 py-7 max-w-3xl mx-auto space-y-8">
            
            {error && (
              <div className="bg-destructive/10 text-destructive px-4 py-2 rounded-md text-sm">
                {error}
              </div>
            )}

            {/* Basic Info */}
            <SettingsSection
              title={t('基本信息')}
              description={t('您的基本个人信息')}
            >
              <SettingsCard divided>
                <SettingsInput
                  label={t('用户名')}
                  description={t('您的登录用户名（不可更改）')}
                  value={profile.username}
                  disabled
                  inCard
                />
                <SettingsInput
                  label={t('昵称')}
                  description={t('其他人看到的名称')}
                  value={profile.nickname || ''}
                  onChange={(v) => updateField('nickname', v)}
                  placeholder={t('请输入昵称')}
                  inCard
                />
                <SettingsInput
                  label={t('头像')}
                  description={t('头像 URL 地址')}
                  value={profile.avatar || ''}
                  onChange={(v) => updateField('avatar', v)}
                  placeholder={t('https://example.com/avatar.png')}
                  inCard
                />
                <SettingsSelect
                  label={t('性别')}
                  description={t('您的性别')}
                  value={profile.gender || ''}
                  onValueChange={(v) => updateField('gender', v)}
                  options={[
                    { value: '', label: t('未设置') },
                    { value: 'male', label: t('男') },
                    { value: 'female', label: t('女') },
                    { value: 'other', label: t('其他') },
                  ]}
                  inCard
                />
                <SettingsInput
                  label={t('生日')}
                  description={t('格式：YYYY-MM-DD')}
                  value={profile.birthday || ''}
                  onChange={(v) => updateField('birthday', v)}
                  placeholder="1990-01-01"
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Contact Info */}
            <SettingsSection
              title={t('联系信息')}
              description={t('您的联系方式')}
            >
              <SettingsCard divided>
                <SettingsInput
                  label={t('邮箱')}
                  description={t('您的邮箱地址（需验证后修改）')}
                  value={profile.email || ''}
                  disabled
                  inCard
                />
                <SettingsInput
                  label={t('手机号')}
                  description={t('您的手机号（需验证后修改）')}
                  value={profile.phone || ''}
                  disabled
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Location */}
            <SettingsSection
              title={t('地区信息')}
              description={t('您所在的地理位置')}
            >
              <SettingsCard divided>
                <SettingsInput
                  label={t('省份')}
                  description={t('您所在的省份')}
                  value={profile.province || ''}
                  onChange={(v) => updateField('province', v)}
                  placeholder={t('例如：广东省')}
                  inCard
                />
                <SettingsInput
                  label={t('城市')}
                  description={t('您所在的城市')}
                  value={profile.city || ''}
                  onChange={(v) => updateField('city', v)}
                  placeholder={t('例如：深圳市')}
                  inCard
                />
                <SettingsInput
                  label={t('区')}
                  description={t('您所在的区')}
                  value={profile.district || ''}
                  onChange={(v) => updateField('district', v)}
                  placeholder={t('例如：南山区')}
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Professional Info */}
            <SettingsSection
              title={t('职业信息')}
              description={t('您的工作和专业背景')}
            >
              <SettingsCard divided>
                <SettingsInput
                  label={t('行业')}
                  description={t('您所从事的行业')}
                  value={profile.industry || ''}
                  onChange={(v) => updateField('industry', v)}
                  placeholder={t('例如：互联网/科技')}
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Bio */}
            <SettingsSection
              title={t('个人简介')}
              description={t('介绍一下您自己')}
            >
              <SettingsCard divided={false}>
                <SettingsTextarea
                  value={profile.bio || ''}
                  onChange={(v) => updateField('bio', v)}
                  placeholder={t('在这里写下您的个人简介...')}
                  rows={5}
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

          </div>
        </ScrollArea>
      </div>
    </div>
  )
}
