import { useState, useEffect } from 'react'
import { useT } from '@/context/LocaleContext'
import {
  Dialog, DialogContent, DialogHeader, DialogTitle, DialogDescription, DialogFooter,
} from '@/components/ui/dialog'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Textarea } from '@/components/ui/textarea'
import {
  Select, SelectTrigger, SelectValue, SelectContent, SelectItem,
} from '@/components/ui/select'
import type { AccountProfile } from '@sprouty-ai/shared/db/types'

const FREQUENCIES = [
  { value: 'daily', label: '每天' },
  { value: '3_per_week', label: '每周3次' },
  { value: 'weekly', label: '每周1次' },
  { value: 'biweekly', label: '每两周1次' },
  { value: 'monthly', label: '每月1次' },
]

interface ProfileEditDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  projectId: string
  profile: AccountProfile | null
  onSave: (data: {
    project_id: string
    niche: string
    sub_niche: string | null
    persona: string | null
    target_audience: string | null
    tone: string | null
    keywords: string | null
    bio: string | null
    content_pillars: string | null
    posting_frequency: string | null
    best_posting_time: string | null
    style_references: string | null
    taboo_topics: string | null
  }) => Promise<unknown>
}

export function ProfileEditDialog({ open, onOpenChange, projectId, profile, onSave }: ProfileEditDialogProps) {
  const t = useT()
  const [saving, setSaving] = useState(false)

  const [niche, setNiche] = useState('')
  const [subNiche, setSubNiche] = useState('')
  const [persona, setPersona] = useState('')
  const [targetAudience, setTargetAudience] = useState('')
  const [tone, setTone] = useState('')
  const [keywords, setKeywords] = useState('')
  const [bio, setBio] = useState('')
  const [contentPillars, setContentPillars] = useState('')
  const [postingFrequency, setPostingFrequency] = useState('')
  const [tabooTopics, setTabooTopics] = useState('')

  useEffect(() => {
    if (profile) {
      setNiche(profile.niche || '')
      setSubNiche(profile.sub_niche || '')
      setPersona(profile.persona || '')
      setTargetAudience(profile.target_audience || '')
      setTone(profile.tone || '')
      setKeywords(profile.keywords || '')
      setBio(profile.bio || '')
      setContentPillars(profile.content_pillars || '')
      setPostingFrequency(profile.posting_frequency || '')
      setTabooTopics(profile.taboo_topics || '')
    } else {
      setNiche('')
      setSubNiche('')
      setPersona('')
      setTargetAudience('')
      setTone('')
      setKeywords('')
      setBio('')
      setContentPillars('')
      setPostingFrequency('')
      setTabooTopics('')
    }
  }, [profile, open])

  const handleSave = async () => {
    if (!niche.trim()) return
    setSaving(true)
    try {
      await onSave({
        project_id: projectId,
        niche: niche.trim(),
        sub_niche: subNiche.trim() || null,
        persona: persona.trim() || null,
        target_audience: targetAudience.trim() || null,
        tone: tone.trim() || null,
        keywords: keywords.trim() || null,
        bio: bio.trim() || null,
        content_pillars: contentPillars.trim() || null,
        posting_frequency: postingFrequency || null,
        best_posting_time: null,
        style_references: null,
        taboo_topics: tabooTopics.trim() || null,
      })
      onOpenChange(false)
    } finally {
      setSaving(false)
    }
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-lg max-h-[85vh] overflow-y-auto">
        <DialogHeader>
          <DialogTitle>{t('编辑账号画像')}</DialogTitle>
          <DialogDescription>{t('设置账号定位和内容策略')}</DialogDescription>
        </DialogHeader>

        <div className="space-y-4">
          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label>{t('领域')} *</Label>
              <Input
                value={niche}
                onChange={(e) => setNiche(e.target.value)}
                placeholder={t('例如：美妆')}
                autoFocus
              />
            </div>
            <div className="space-y-2">
              <Label>{t('细分领域')}</Label>
              <Input
                value={subNiche}
                onChange={(e) => setSubNiche(e.target.value)}
                placeholder={t('例如：护肤成分')}
              />
            </div>
          </div>

          <div className="space-y-2">
            <Label>{t('人设')}</Label>
            <Input
              value={persona}
              onChange={(e) => setPersona(e.target.value)}
              placeholder={t('例如：专业测评博主')}
            />
          </div>

          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label>{t('目标受众')}</Label>
              <Input
                value={targetAudience}
                onChange={(e) => setTargetAudience(e.target.value)}
                placeholder={t('例如：18-30岁女性')}
              />
            </div>
            <div className="space-y-2">
              <Label>{t('调性')}</Label>
              <Input
                value={tone}
                onChange={(e) => setTone(e.target.value)}
                placeholder={t('例如：专业、亲切')}
              />
            </div>
          </div>

          <div className="space-y-2">
            <Label>{t('简介')}</Label>
            <Textarea
              value={bio}
              onChange={(e) => setBio(e.target.value)}
              placeholder={t('账号简介...')}
              rows={2}
            />
          </div>

          <div className="space-y-2">
            <Label>{t('关键词')}</Label>
            <Input
              value={keywords}
              onChange={(e) => setKeywords(e.target.value)}
              placeholder={t('逗号分隔，例如：护肤,成分,测评')}
            />
            <p className="text-xs text-muted-foreground">{t('多个关键词用逗号分隔')}</p>
          </div>

          <div className="space-y-2">
            <Label>{t('内容支柱')}</Label>
            <Input
              value={contentPillars}
              onChange={(e) => setContentPillars(e.target.value)}
              placeholder={t('逗号分隔，例如：产品测评,成分科普,护肤教程')}
            />
            <p className="text-xs text-muted-foreground">{t('多个内容支柱用逗号分隔')}</p>
          </div>

          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label>{t('发布频率')}</Label>
              <Select value={postingFrequency} onValueChange={setPostingFrequency}>
                <SelectTrigger>
                  <SelectValue placeholder={t('选择频率')} />
                </SelectTrigger>
                <SelectContent className="z-modal">
                  {FREQUENCIES.map((f) => (
                    <SelectItem key={f.value} value={f.value}>{t(f.label)}</SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>
            <div className="space-y-2">
              <Label>{t('禁忌话题')}</Label>
              <Input
                value={tabooTopics}
                onChange={(e) => setTabooTopics(e.target.value)}
                placeholder={t('逗号分隔')}
              />
            </div>
          </div>
        </div>

        <DialogFooter>
          <Button variant="outline" onClick={() => onOpenChange(false)} disabled={saving}>
            {t('取消')}
          </Button>
          <Button onClick={handleSave} disabled={!niche.trim() || saving}>
            {saving ? t('保存中...') : t('保存')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
