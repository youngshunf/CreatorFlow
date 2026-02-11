import { useState } from 'react'
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
import type { Platform } from '@sprouty-ai/shared/db/types'
import { PLATFORM_LIST } from '@sprouty-ai/shared/db/types'

const PLATFORMS = PLATFORM_LIST.map(p => ({ value: p.id, label: p.label }))

/** 读取文件为 data URL */
function readFileAsDataURL(file: File): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader()
    reader.onload = () => resolve(reader.result as string)
    reader.onerror = reject
    reader.readAsDataURL(file)
  })
}

interface CreateProjectDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  onCreateProject: (data: {
    name: string
    platform: Platform
    description: string | null
    platforms: string | null
    avatar_path: string | null
  }) => Promise<unknown>
  onUpsertProfile?: (projectId: string, data: {
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

export function CreateProjectDialog({ open, onOpenChange, onCreateProject, onUpsertProfile }: CreateProjectDialogProps) {
  const t = useT()
  const [step, setStep] = useState(1)
  const [saving, setSaving] = useState(false)

  // 步骤 1：基本信息
  const [name, setName] = useState('')
  const [platform, setPlatform] = useState<Platform>('xiaohongshu')
  const [selectedPlatforms, setSelectedPlatforms] = useState<Platform[]>([])
  const [description, setDescription] = useState('')
  const [avatarPreview, setAvatarPreview] = useState('')

  // 步骤 2：账号画像
  const [niche, setNiche] = useState('')
  const [persona, setPersona] = useState('')
  const [targetAudience, setTargetAudience] = useState('')
  const [tone, setTone] = useState('')
  const [keywords, setKeywords] = useState('')

  const reset = () => {
    setStep(1)
    setName('')
    setPlatform('xiaohongshu')
    setSelectedPlatforms([])
    setDescription('')
    setAvatarPreview('')
    setNiche('')
    setPersona('')
    setTargetAudience('')
    setTone('')
    setKeywords('')
    setSaving(false)
  }

  const handleClose = (v: boolean) => {
    if (!v) reset()
    onOpenChange(v)
  }

  const handleNext = () => {
    if (!name.trim()) return
    setStep(2)
  }

  const handleAvatarChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    if (!file.type.startsWith('image/')) return
    if (file.size > 5 * 1024 * 1024) return
    const dataUrl = await readFileAsDataURL(file)
    setAvatarPreview(dataUrl)
    e.target.value = ''
  }

  const togglePlatform = (p: Platform) => {
    setSelectedPlatforms((prev) =>
      prev.includes(p) ? prev.filter((x) => x !== p) : [...prev, p]
    )
  }

  const handleCreate = async () => {
    if (!name.trim()) return
    setSaving(true)
    try {
      const project = await onCreateProject({
        name: name.trim(),
        platform,
        description: description.trim() || null,
        platforms: selectedPlatforms.length > 0 ? JSON.stringify(selectedPlatforms) : null,
        avatar_path: avatarPreview || null,
      }) as { id: string } | null

      // 如果有画像数据且创建成功，保存画像
      if (project && onUpsertProfile && niche.trim()) {
        await onUpsertProfile(project.id, {
          project_id: project.id,
          niche: niche.trim(),
          sub_niche: null,
          persona: persona.trim() || null,
          target_audience: targetAudience.trim() || null,
          tone: tone.trim() || null,
          keywords: keywords.trim() || null,
          bio: null,
          content_pillars: null,
          posting_frequency: null,
          best_posting_time: null,
          style_references: null,
          taboo_topics: null,
        })
      }
      handleClose(false)
    } finally {
      setSaving(false)
    }
  }

  return (
    <Dialog open={open} onOpenChange={handleClose}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle>{t('新建项目')}</DialogTitle>
          <DialogDescription>
            {step === 1 ? t('填写项目基本信息') : t('设置账号画像（可选）')}
          </DialogDescription>
        </DialogHeader>

        {step === 1 ? (
          <div className="space-y-4">
            {/* 头像 */}
            <div className="space-y-2">
              <Label>{t('项目头像')}</Label>
              <div className="flex items-center gap-3">
                <label className="relative w-14 h-14 rounded-full bg-muted/60 overflow-hidden cursor-pointer group flex-shrink-0">
                  {avatarPreview ? (
                    <img src={avatarPreview} alt="" className="w-full h-full object-cover" />
                  ) : (
                    <div className="w-full h-full flex items-center justify-center text-muted-foreground">
                      <svg className="w-6 h-6" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                        <path strokeLinecap="round" strokeLinejoin="round" d="M6.827 6.175A2.31 2.31 0 0 1 5.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 0 0 2.25 2.25h15A2.25 2.25 0 0 0 21.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 0 0-1.134-.175 2.31 2.31 0 0 1-1.64-1.055l-.822-1.316a2.192 2.192 0 0 0-1.736-1.039 48.774 48.774 0 0 0-5.232 0 2.192 2.192 0 0 0-1.736 1.039l-.821 1.316Z" />
                        <path strokeLinecap="round" strokeLinejoin="round" d="M16.5 12.75a4.5 4.5 0 1 1-9 0 4.5 4.5 0 0 1 9 0Z" />
                      </svg>
                    </div>
                  )}
                  <div className="absolute inset-0 bg-black/40 opacity-0 group-hover:opacity-100 transition-opacity flex items-center justify-center">
                    <svg className="w-4 h-4 text-white" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="M6.827 6.175A2.31 2.31 0 0 1 5.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 0 0 2.25 2.25h15A2.25 2.25 0 0 0 21.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 0 0-1.134-.175 2.31 2.31 0 0 1-1.64-1.055l-.822-1.316a2.192 2.192 0 0 0-1.736-1.039 48.774 48.774 0 0 0-5.232 0 2.192 2.192 0 0 0-1.736 1.039l-.821 1.316Z" />
                      <path strokeLinecap="round" strokeLinejoin="round" d="M16.5 12.75a4.5 4.5 0 1 1-9 0 4.5 4.5 0 0 1 9 0Z" />
                    </svg>
                  </div>
                  <input
                    type="file"
                    accept="image/jpeg,image/png,image/gif,image/webp"
                    className="hidden"
                    onChange={handleAvatarChange}
                  />
                </label>
                <span className="text-xs text-muted-foreground">{t('点击上传头像')}</span>
              </div>
            </div>

            <div className="space-y-2">
              <Label>{t('项目名称')} *</Label>
              <Input
                value={name}
                onChange={(e) => setName(e.target.value)}
                placeholder={t('例如：我的小红书账号')}
                autoFocus
              />
            </div>

            <div className="space-y-2">
              <Label>{t('主平台')}</Label>
              <Select value={platform} onValueChange={(v) => setPlatform(v as Platform)}>
                <SelectTrigger>
                  <SelectValue />
                </SelectTrigger>
                <SelectContent className="z-modal">
                  {PLATFORMS.map((p) => (
                    <SelectItem key={p.value} value={p.value}>{t(p.label)}</SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label>{t('多平台分发')}</Label>
              <div className="flex flex-wrap gap-2">
                {PLATFORMS.filter((p) => p.value !== platform).map((p) => (
                  <button
                    key={p.value}
                    type="button"
                    onClick={() => togglePlatform(p.value)}
                    className={`inline-flex items-center rounded-full px-2.5 py-1 text-xs font-medium transition-colors ${
                      selectedPlatforms.includes(p.value)
                        ? 'bg-foreground text-background'
                        : 'bg-muted/60 text-muted-foreground hover:bg-muted'
                    }`}
                  >
                    {t(p.label)}
                  </button>
                ))}
              </div>
              <p className="text-xs text-muted-foreground">{t('选择需要同步分发的其他平台')}</p>
            </div>

            <div className="space-y-2">
              <Label>{t('描述')}</Label>
              <Textarea
                value={description}
                onChange={(e) => setDescription(e.target.value)}
                placeholder={t('简要描述项目定位...')}
                rows={2}
              />
            </div>
          </div>
        ) : (
          <div className="space-y-4">
            <div className="space-y-2">
              <Label>{t('领域')}</Label>
              <Input
                value={niche}
                onChange={(e) => setNiche(e.target.value)}
                placeholder={t('例如：美妆、科技、美食')}
                autoFocus
              />
            </div>
            <div className="space-y-2">
              <Label>{t('人设')}</Label>
              <Input
                value={persona}
                onChange={(e) => setPersona(e.target.value)}
                placeholder={t('例如：专业测评博主')}
              />
            </div>
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
                placeholder={t('例如：专业、亲切、幽默')}
              />
            </div>
            <div className="space-y-2">
              <Label>{t('关键词')}</Label>
              <Input
                value={keywords}
                onChange={(e) => setKeywords(e.target.value)}
                placeholder={t('逗号分隔，例如：护肤,成分,测评')}
              />
            </div>
          </div>
        )}

        <DialogFooter>
          {step === 2 && (
            <Button variant="outline" onClick={() => setStep(1)} disabled={saving}>
              {t('上一步')}
            </Button>
          )}
          {step === 1 ? (
            <Button onClick={handleNext} disabled={!name.trim()}>
              {t('下一步')}
            </Button>
          ) : (
            <Button onClick={handleCreate} disabled={saving}>
              {saving ? t('创建中...') : t('创建项目')}
            </Button>
          )}
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
