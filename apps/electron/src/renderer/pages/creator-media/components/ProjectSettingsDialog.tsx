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
import type { Project, Platform, UpdateProject, AccountProfile } from '@sprouty-ai/shared/db/types'

/** 构建 AI 编辑项目的 prompt，返回完整 prompt 和用于隐藏元数据的 badge */
function buildAIEditPrompt(project: Project, profile: AccountProfile | null, userInput: string) {
  const currentInfo = JSON.stringify({
    project: { name: project.name, platform: project.platform, description: project.description, platforms: project.platforms },
    profile: profile ? {
      niche: profile.niche, sub_niche: profile.sub_niche, persona: profile.persona,
      target_audience: profile.target_audience, tone: profile.tone, keywords: profile.keywords,
      content_pillars: profile.content_pillars, posting_frequency: profile.posting_frequency,
      taboo_topics: profile.taboo_topics,
    } : null,
  }, null, 2)

  const metadataSection = `<creator_media_task>
<action>update_project</action>
<project_id>${project.id}</project_id>
<current_data>${currentInfo}</current_data>
<context>用户想要通过 AI 更新自媒体项目信息和账号画像。请根据用户提供的信息（可能是账号链接、账号描述、或修改指令），分析并更新项目信息和画像。

你需要：
1. 如果用户提供了账号链接，尝试访问并分析该账号的公开信息
2. 根据分析结果，更新项目和画像数据
3. 项目字段：name、platform（xiaohongshu/douyin/bilibili/wechat/zhihu/weibo/x）、description、platforms（JSON 数组）
4. 画像字段：niche、sub_niche、persona、target_audience、tone、keywords（逗号分隔）、content_pillars（逗号分隔）、posting_frequency（daily/3_per_week/weekly/biweekly/monthly）、taboo_topics

请使用 creatorMedia 技能完成更新。如果无法直接更新，请输出结构化的 JSON 数据供用户参考。</context>
</creator_media_task>

`
  const prompt = metadataSection + userInput
  const badges = [{
    type: 'context' as const,
    label: '编辑项目',
    rawText: metadataSection,
    start: 0,
    end: metadataSection.length,
    collapsedLabel: '编辑项目',
  }]
  return { prompt, badges }
}

const PLATFORMS: { value: Platform; label: string }[] = [
  { value: 'xiaohongshu', label: '小红书' },
  { value: 'douyin', label: '抖音' },
  { value: 'bilibili', label: 'B站' },
  { value: 'wechat', label: '微信公众号' },
  { value: 'zhihu', label: '知乎' },
  { value: 'weibo', label: '微博' },
  { value: 'x', label: 'X (Twitter)' },
]

/** 读取文件为 data URL */
function readFileAsDataURL(file: File): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader()
    reader.onload = () => resolve(reader.result as string)
    reader.onerror = reject
    reader.readAsDataURL(file)
  })
}

interface ProjectSettingsDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  project: Project
  profile?: AccountProfile | null
  onUpdate: (projectId: string, data: UpdateProject) => Promise<unknown>
  onDelete: (projectId: string) => Promise<unknown>
}

export function ProjectSettingsDialog({ open, onOpenChange, project, profile, onUpdate, onDelete }: ProjectSettingsDialogProps) {
  const t = useT()
  const [mode, setMode] = useState<'manual' | 'ai'>('manual')
  const [aiInput, setAiInput] = useState('')
  const [name, setName] = useState(project.name)
  const [platform, setPlatform] = useState<Platform>(project.platform)
  const [selectedPlatforms, setSelectedPlatforms] = useState<Platform[]>([])
  const [description, setDescription] = useState(project.description || '')
  const [avatarPreview, setAvatarPreview] = useState(project.avatar_path || '')
  const [saving, setSaving] = useState(false)
  const [confirmDelete, setConfirmDelete] = useState(false)

  useEffect(() => {
    setName(project.name)
    setPlatform(project.platform)
    setDescription(project.description || '')
    setAvatarPreview(project.avatar_path || '')
    setConfirmDelete(false)
    setMode('manual')
    setAiInput('')
    // 解析已有的多平台数据
    try {
      const parsed = project.platforms ? JSON.parse(project.platforms) : []
      setSelectedPlatforms(Array.isArray(parsed) ? parsed : [])
    } catch {
      setSelectedPlatforms([])
    }
  }, [project, open])

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

  const handleSave = async () => {
    if (!name.trim()) return
    setSaving(true)
    try {
      await onUpdate(project.id, {
        name: name.trim(),
        platform,
        description: description.trim() || null,
        platforms: selectedPlatforms.length > 0 ? JSON.stringify(selectedPlatforms) : null,
        avatar_path: avatarPreview || null,
      })
      onOpenChange(false)
    } finally {
      setSaving(false)
    }
  }

  const handleDelete = async () => {
    if (!confirmDelete) {
      setConfirmDelete(true)
      return
    }
    setSaving(true)
    try {
      await onDelete(project.id)
      onOpenChange(false)
    } finally {
      setSaving(false)
    }
  }

  const handleAIEdit = () => {
    if (!aiInput.trim()) return
    const { prompt, badges } = buildAIEditPrompt(project, profile || null, aiInput.trim())
    const encodedInput = encodeURIComponent(prompt)
    const encodedBadges = encodeURIComponent(JSON.stringify(badges))
    const url = `sproutyai://action/new-chat?input=${encodedInput}&send=true&mode=allow-all&badges=${encodedBadges}`
    window.electronAPI.openUrl(url)
    onOpenChange(false)
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle>{t('项目设置')}</DialogTitle>
          <DialogDescription>
            {mode === 'ai'
              ? t('描述修改内容或粘贴账号链接，AI 自动更新项目信息')
              : t('编辑项目信息或删除项目')}
          </DialogDescription>
        </DialogHeader>

        {/* 模式切换 */}
        <div className="flex items-center gap-1 rounded-lg bg-muted/40 p-0.5">
          <button
            type="button"
            onClick={() => setMode('manual')}
            className={`flex-1 rounded-md px-3 py-1.5 text-xs font-medium transition-colors ${
              mode === 'manual'
                ? 'bg-background text-foreground shadow-sm'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            {t('手动编辑')}
          </button>
          <button
            type="button"
            onClick={() => setMode('ai')}
            className={`flex-1 rounded-md px-3 py-1.5 text-xs font-medium transition-colors ${
              mode === 'ai'
                ? 'bg-background text-foreground shadow-sm'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            <span className="flex items-center justify-center gap-1">
              <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456Z" />
              </svg>
              {t('AI 智能编辑')}
            </span>
          </button>
        </div>

        {mode === 'ai' ? (
          <div className="space-y-4">
            <div className="space-y-2">
              <Label>{t('描述你想要的修改')}</Label>
              <Textarea
                value={aiInput}
                onChange={(e) => setAiInput(e.target.value)}
                placeholder={t('例如：\n• 根据这个链接更新项目信息 https://...\n• 把定位改为科技数码测评\n• 增加 B站 和知乎作为分发平台')}
                rows={5}
                autoFocus
              />
              <p className="text-xs text-muted-foreground">
                {t('支持输入修改描述或账号链接，AI 将自动分析并更新项目')}
              </p>
            </div>
          </div>
        ) : (
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
            <Input value={name} onChange={(e) => setName(e.target.value)} />
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
                  {t(PLATFORMS.find((x) => x.value === p.value)?.label || p.value)}
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
              rows={2}
            />
          </div>
        </div>
        )}

        <DialogFooter className="flex items-center justify-between sm:justify-between">
          {mode === 'ai' ? (
            <>
              <Button
                variant="destructive"
                onClick={handleDelete}
                disabled={saving}
                className="mr-auto"
              >
                {confirmDelete ? t('确认删除') : t('删除项目')}
              </Button>
              <Button onClick={handleAIEdit} disabled={!aiInput.trim()}>
                <svg className="w-4 h-4 mr-1.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456Z" />
                </svg>
                {t('AI 更新')}
              </Button>
            </>
          ) : (
            <>
              <Button
                variant="destructive"
                onClick={handleDelete}
                disabled={saving}
                className="mr-auto"
              >
                {confirmDelete ? t('确认删除') : t('删除项目')}
              </Button>
              <div className="flex gap-2">
                <Button variant="outline" onClick={() => onOpenChange(false)} disabled={saving}>
                  {t('取消')}
                </Button>
                <Button onClick={handleSave} disabled={!name.trim() || saving}>
                  {saving ? t('保存中...') : t('保存')}
                </Button>
              </div>
            </>
          )}
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
