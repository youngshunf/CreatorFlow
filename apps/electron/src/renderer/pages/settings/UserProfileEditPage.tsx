/**
 * UserProfileEditPage
 *
 * Separate page for editing user profile with region picker and date picker.
 * Save button is at the bottom of the page.
 */

import { useState, useEffect, useCallback } from 'react'
import { ArrowLeft } from 'lucide-react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { useT } from '@/context/LocaleContext'
import { navigate, routes } from '@/lib/navigate'
import { Spinner } from '@creator-flow/ui'
import type { DetailsPageMeta } from '@/lib/navigation-registry'

import {
  SettingsSection,
  SettingsCard,
  SettingsCardContent,
  SettingsInput,
  SettingsSelect,
  SettingsTextarea,
} from '@/components/settings'

import { DatePicker } from '@/components/ui/date-picker'
import { RegionPicker, type RegionValue } from '@/components/ui/region-picker'
import { AvatarUploader } from '@/components/ui/avatar-uploader'
import { userApi, type UserProfile, type UpdateUserProfileParams } from '@/api/user'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'user-profile-edit',
}

export default function UserProfileEditPage() {
  const t = useT()
  const [isLoading, setIsLoading] = useState(true)
  const [isSaving, setIsSaving] = useState(false)
  const [profile, setProfile] = useState<UserProfile | null>(null)
  const [editedProfile, setEditedProfile] = useState<Partial<UpdateUserProfileParams>>({})
  const [error, setError] = useState<string | null>(null)

  // Load user profile on mount
  useEffect(() => {
    const loadProfile = async () => {
      try {
        setIsLoading(true)
        const data = await userApi.getCurrentUser()
        setProfile(data)
        // Initialize edited profile with current values
        setEditedProfile({
          nickname: data.nickname || '',
          avatar: data.avatar || '',
          gender: data.gender || '',
          birthday: data.birthday || '',
          province: data.province || '',
          city: data.city || '',
          district: data.district || '',
          industry: data.industry || '',
          bio: data.bio || '',
        })
      } catch (err: any) {
        console.error('Failed to load user profile:', err)
        setError(err.message || t('加载用户信息失败'))
      } finally {
        setIsLoading(false)
      }
    }
    loadProfile()
  }, [t])

  // Navigate back to profile view
  const handleBack = useCallback(() => {
    navigate(routes.view.settings('user-profile'))
  }, [])

  // Save changes
  const handleSave = useCallback(async () => {
    if (!profile) return

    try {
      setIsSaving(true)
      setError(null)
      await userApi.updateProfile(editedProfile)
      // Navigate back to profile view
      navigate(routes.view.settings('user-profile'))
    } catch (err: any) {
      console.error('Failed to save profile:', err)
      setError(err.message || t('保存失败'))
    } finally {
      setIsSaving(false)
    }
  }, [profile, editedProfile, t])

  // Update edit field
  const updateEditField = useCallback(<K extends keyof UpdateUserProfileParams>(
    field: K,
    value: UpdateUserProfileParams[K]
  ) => {
    setEditedProfile(prev => ({ ...prev, [field]: value }))
  }, [])

  // Handle region change
  const handleRegionChange = useCallback((value: RegionValue) => {
    setEditedProfile(prev => ({
      ...prev,
      province: value.province || '',
      city: value.city || '',
      district: value.district || '',
    }))
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
          <Button variant="outline" className="mt-4" onClick={handleBack}>
            {t('返回')}
          </Button>
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
        title={t('编辑资料')}
        actions={
          <Button variant="ghost" size="sm" onClick={handleBack}>
            <ArrowLeft className="size-4 mr-1.5" />
            {t('返回')}
          </Button>
        }
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
                  onChange={() => {}}
                  disabled
                  inCard
                />
                <SettingsInput
                  label={t('昵称')}
                  description={t('其他人看到的名称')}
                  value={editedProfile.nickname || ''}
                  onChange={(v) => updateEditField('nickname', v)}
                  placeholder={t('请输入昵称')}
                  inCard
                />
                {/* 头像上传 */}
                <SettingsCardContent>
                  <div className="space-y-2">
                    <div className="space-y-0.5">
                      <label className="text-sm font-medium">{t('头像')}</label>
                      <p className="text-sm text-muted-foreground">{t('点击上传或更换您的头像')}</p>
                    </div>
                    <AvatarUploader
                      value={editedProfile.avatar || ''}
                      onChange={(url) => updateEditField('avatar', url)}
                      onUpload={userApi.uploadAvatar}
                      size={80}
                    />
                  </div>
                </SettingsCardContent>
                <SettingsSelect
                  label={t('性别')}
                  description={t('您的性别')}
                  value={editedProfile.gender || 'unset'}
                  onValueChange={(v) => updateEditField('gender', v === 'unset' ? '' : v)}
                  options={[
                    { value: 'unset', label: t('未设置') },
                    { value: 'male', label: t('男') },
                    { value: 'female', label: t('女') },
                    { value: 'other', label: t('其他') },
                  ]}
                  inCard
                />
                {/* Birthday with DatePicker */}
                <SettingsCardContent>
                  <div className="space-y-2">
                    <div className="space-y-0.5">
                      <label className="text-sm font-medium">{t('生日')}</label>
                      <p className="text-sm text-muted-foreground">{t('选择您的出生日期')}</p>
                    </div>
                    <DatePicker
                      value={editedProfile.birthday || ''}
                      onChange={(v) => updateEditField('birthday', v)}
                      placeholder={t('选择日期')}
                      maxDate={new Date()}
                    />
                  </div>
                </SettingsCardContent>
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
                  onChange={() => {}}
                  disabled
                  inCard
                />
                <SettingsInput
                  label={t('手机号')}
                  description={t('您的手机号（需验证后修改）')}
                  value={profile.phone || ''}
                  onChange={() => {}}
                  disabled
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Location with RegionPicker */}
            <SettingsSection
              title={t('地区信息')}
              description={t('您所在的地理位置')}
            >
              <SettingsCard divided={false}>
                <SettingsCardContent>
                  <div className="space-y-2">
                    <div className="space-y-0.5">
                      <label className="text-sm font-medium">{t('所在地区')}</label>
                      <p className="text-sm text-muted-foreground">{t('选择省市区')}</p>
                    </div>
                    <RegionPicker
                      value={{
                        province: editedProfile.province,
                        city: editedProfile.city,
                        district: editedProfile.district,
                      }}
                      onChange={handleRegionChange}
                      placeholder={t('请选择地区')}
                    />
                  </div>
                </SettingsCardContent>
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
                  value={editedProfile.industry || ''}
                  onChange={(v) => updateEditField('industry', v)}
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
                  value={editedProfile.bio || ''}
                  onChange={(v) => updateEditField('bio', v)}
                  placeholder={t('在这里写下您的个人简介...')}
                  rows={5}
                  inCard
                />
              </SettingsCard>
            </SettingsSection>

            {/* Save Button at Bottom */}
            <div className="pt-4 pb-8">
              <Button
                className="w-full"
                size="lg"
                onClick={handleSave}
                disabled={isSaving}
              >
                {isSaving && <Spinner className="size-4 mr-2" />}
                {t('保存修改')}
              </Button>
            </div>

          </div>
        </ScrollArea>
      </div>
    </div>
  )
}
