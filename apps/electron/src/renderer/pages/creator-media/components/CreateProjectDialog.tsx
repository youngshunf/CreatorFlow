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
import { PLATFORM_LIST, PLATFORM_IDS } from '@sprouty-ai/shared/db/types'

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

/** 构建 AI 创建项目的 prompt，返回完整 prompt 和用于隐藏元数据的 badge */
function buildAICreatePrompt(userInput: string) {
  const metadataSection = `<creator_media_task>
<action>create_project</action>
<context>用户想要通过 AI 智能创建一个自媒体项目。请根据用户提供的信息（可能是账号链接、账号描述、或领域描述），分析并生成完整的项目信息和账号画像。

你需要：
1. 如果用户提供了账号链接，尝试访问并分析该账号的公开信息
2. 根据分析结果，使用 sqlite3 命令直接将项目和画像写入数据库
3. 所有字段都必须填写，不能为 NULL —— 如果用户未提供某些信息，请根据上下文智能推断合理的值

数据库路径：.sprouty-ai/db/creator.db（相对于工作目录）

操作步骤（必须严格按照以下顺序执行）：

步骤1 - 生成 UUID 和头像：
\`\`\`bash
PROJECT_ID=$(uuidgen | tr '[:upper:]' '[:lower:]')
PROFILE_ID=$(uuidgen | tr '[:upper:]' '[:lower:]')
\`\`\`

头像处理（avatar_path）：
- 如果用户提供了图片 URL，直接使用该 URL 作为 avatar_path
- 否则，根据项目名称和领域生成一个 SVG 头像文件，保存到 .sprouty-ai/avatars/<PROJECT_ID>.svg，avatar_path 填写该相对路径
- SVG 头像要求：简洁美观，使用领域相关图标，配合适合该领域的渐变背景色，尺寸 128x128

步骤2 - 创建项目并设为活跃：
\`\`\`bash
sqlite3 .sprouty-ai/db/creator.db <<'SQL'
UPDATE projects SET is_active = 0, updated_at = CURRENT_TIMESTAMP WHERE is_active = 1;
INSERT INTO projects (id, name, description, platform, platforms, avatar_path, is_active)
VALUES ('<PROJECT_ID>', '<项目名称>', '<项目描述>', '<主平台>', '<平台JSON数组>', '<头像路径>', 1);
SQL
\`\`\`

步骤3 - 创建账号画像（所有字段必填）：
\`\`\`bash
sqlite3 .sprouty-ai/db/creator.db <<'SQL'
INSERT INTO account_profiles (id, project_id, niche, sub_niche, persona, target_audience, tone, keywords, bio, content_pillars, posting_frequency, best_posting_time, style_references, taboo_topics, pillar_weights)
VALUES ('<PROFILE_ID>', '<PROJECT_ID>', '<领域>', '<细分领域>', '<人设定位>', '<目标受众>', '<内容调性>', '<关键词逗号分隔>', '<账号简介>', '<内容支柱逗号分隔>', '<发布频率>', '<最佳发布时间>', '<风格参考>', '<禁忌话题>', '<支柱权重JSON>');
SQL
\`\`\`

字段说明：
- platform 取值：${PLATFORM_IDS}
- platforms：JSON 数组字符串，如 '["xiaohongshu","douyin"]'
- posting_frequency 取值：daily / 3_per_week / weekly / biweekly / monthly
- best_posting_time：如 '20:00' 或 '12:00,20:00'
- pillar_weights：JSON 对象，键为内容支柱名称，值为权重(0-1)，所有权重之和为1，如 '{"知识分享":0.4,"产品评测":0.3,"日常vlog":0.3}'
- style_references：参考的同领域优秀账号或内容风格，逗号分隔
- taboo_topics：该领域应避免的敏感话题，逗号分隔
- 文本中的单引号需要转义为两个单引号（SQL 标准）

重要：
- 所有字段都必须有值，不允许 NULL，请根据用户信息和领域知识智能推断
- 必须使用 sqlite3 命令操作数据库，不要调用任何 MCP 工具
- 创建 SVG 头像前先确保目录存在：mkdir -p .sprouty-ai/avatars
- 创建完成后告知用户项目已创建成功</context>
</creator_media_task>

`
  const prompt = metadataSection + userInput
  const badges = [{
    type: 'context' as const,
    label: '创建项目',
    rawText: metadataSection,
    start: 0,
    end: metadataSection.length,
    collapsedLabel: '创建项目',
  }]
  return { prompt, badges }
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
  const [mode, setMode] = useState<'manual' | 'ai'>('ai')
  const [step, setStep] = useState(1)
  const [saving, setSaving] = useState(false)

  // AI 模式
  const [aiInput, setAiInput] = useState('')

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
    setMode('manual')
    setStep(1)
    setAiInput('')
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

  const handleAICreate = () => {
    if (!aiInput.trim()) return
    const { prompt, badges } = buildAICreatePrompt(aiInput.trim())
    const encodedInput = encodeURIComponent(prompt)
    const encodedBadges = encodeURIComponent(JSON.stringify(badges))
    const url = `sproutyai://action/new-chat?input=${encodedInput}&send=true&mode=allow-all&badges=${encodedBadges}`
    window.electronAPI.openUrl(url)
    handleClose(false)
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
            {mode === 'ai'
              ? t('描述你的账号或粘贴主页链接，AI 自动生成项目信息')
              : step === 1 ? t('填写项目基本信息') : t('设置账号画像（可选）')}
          </DialogDescription>
        </DialogHeader>

        {/* 模式切换 */}
        <div className="flex items-center gap-1 rounded-lg bg-muted/40 p-0.5">
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
              {t('AI 智能填写')}
            </span>
          </button>
          <button
            type="button"
            onClick={() => setMode('manual')}
            className={`flex-1 rounded-md px-3 py-1.5 text-xs font-medium transition-colors ${
              mode === 'manual'
                ? 'bg-background text-foreground shadow-sm'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            {t('手动填写')}
          </button>
        </div>

        {mode === 'ai' ? (
          <div className="space-y-4">
            <div className="space-y-2">
              <Label>{t('描述你的账号')}</Label>
              <Textarea
                value={aiInput}
                onChange={(e) => setAiInput(e.target.value)}
                placeholder={t('例如：\n• 我是一个小红书美妆博主，主要做护肤品测评\n• https://www.xiaohongshu.com/user/profile/xxx\n• 我想做一个抖音科技数码账号，面向年轻男性')}
                rows={5}
                autoFocus
              />
              <p className="text-xs text-muted-foreground">
                {t('支持输入账号描述、主页链接或定位方向，AI 将自动分析并创建项目')}
              </p>
            </div>
          </div>
        ) : step === 1 ? (
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
          {mode === 'ai' ? (
            <Button onClick={handleAICreate} disabled={!aiInput.trim()}>
              <svg className="w-4 h-4 mr-1.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456Z" />
              </svg>
              {t('AI 创建')}
            </Button>
          ) : (
            <>
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
            </>
          )}
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
