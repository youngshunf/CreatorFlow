import { useState, useEffect } from 'react'
import { useT } from '@/context/LocaleContext'
import type { ViralPattern, ViralPatternCategory, ViralPatternSource, Platform } from '@sprouty-ai/shared/db/types'
import { PLATFORM_LIST } from '@sprouty-ai/shared/db/types'

const CATEGORIES: { value: ViralPatternCategory; label: string }[] = [
  { value: 'hook', label: '开头钩子' },
  { value: 'structure', label: '内容结构' },
  { value: 'title', label: '标题模式' },
  { value: 'cta', label: '行动号召' },
  { value: 'visual', label: '视觉风格' },
  { value: 'rhythm', label: '节奏韵律' },
]

const SOURCES: { value: ViralPatternSource; label: string }[] = [
  { value: 'manual', label: '手动添加' },
  { value: 'competitor_analysis', label: '竞品分析' },
  { value: 'ai_discovered', label: 'AI 发现' },
]

const PLATFORMS: { value: Platform | ''; label: string }[] = [
  { value: '', label: '全平台' },
  ...PLATFORM_LIST.map(p => ({ value: p.id as Platform, label: p.shortLabel })),
]

interface ViralPatternDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  pattern: ViralPattern | null // null = 新增模式
  projectId: string | null
  onSave: (data: {
    project_id: string | null
    platform: Platform | null
    category: ViralPatternCategory
    name: string
    description: string | null
    template: string | null
    source: ViralPatternSource | null
  }) => void
}

/**
 * 爆款模式新增/编辑对话框
 */
export function ViralPatternDialog({ open, onOpenChange, pattern, projectId, onSave }: ViralPatternDialogProps) {
  const t = useT()

  const [name, setName] = useState('')
  const [category, setCategory] = useState<ViralPatternCategory>('hook')
  const [platform, setPlatform] = useState<Platform | ''>('')
  const [description, setDescription] = useState('')
  const [template, setTemplate] = useState('')
  const [source, setSource] = useState<ViralPatternSource>('manual')

  // 编辑模式时填充数据
  useEffect(() => {
    if (pattern) {
      setName(pattern.name)
      setCategory(pattern.category)
      setPlatform(pattern.platform || '')
      setDescription(pattern.description || '')
      setTemplate(pattern.template || '')
      setSource(pattern.source || 'manual')
    } else {
      setName('')
      setCategory('hook')
      setPlatform('')
      setDescription('')
      setTemplate('')
      setSource('manual')
    }
  }, [pattern, open])

  const handleSubmit = () => {
    if (!name.trim()) return
    onSave({
      project_id: projectId,
      platform: platform || null,
      category,
      name: name.trim(),
      description: description.trim() || null,
      template: template.trim() || null,
      source,
    })
    onOpenChange(false)
  }

  if (!open) return null

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center">
      {/* 遮罩 */}
      <div className="absolute inset-0 bg-black/50" onClick={() => onOpenChange(false)} />

      {/* 对话框 */}
      <div className="relative w-full max-w-lg rounded-lg border border-border bg-background shadow-lg mx-4">
        <div className="px-6 py-4 border-b border-border/40">
          <h2 className="text-base font-semibold text-foreground">
            {pattern ? t('编辑爆款模式') : t('新增爆款模式')}
          </h2>
        </div>

        <div className="px-6 py-4 space-y-4 max-h-[60vh] overflow-auto">
          {/* 名称 */}
          <div>
            <label className="block text-xs font-medium text-muted-foreground mb-1">{t('模式名称')} *</label>
            <input
              type="text"
              value={name}
              onChange={(e) => setName(e.target.value)}
              placeholder={t('例如：3秒悬念开头')}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-foreground/20"
            />
          </div>

          {/* 分类 + 平台 */}
          <div className="grid grid-cols-2 gap-3">
            <div>
              <label className="block text-xs font-medium text-muted-foreground mb-1">{t('分类')}</label>
              <select
                value={category}
                onChange={(e) => setCategory(e.target.value as ViralPatternCategory)}
                className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-foreground/20"
              >
                {CATEGORIES.map((c) => (
                  <option key={c.value} value={c.value}>{t(c.label)}</option>
                ))}
              </select>
            </div>
            <div>
              <label className="block text-xs font-medium text-muted-foreground mb-1">{t('适用平台')}</label>
              <select
                value={platform}
                onChange={(e) => setPlatform(e.target.value as Platform | '')}
                className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-foreground/20"
              >
                {PLATFORMS.map((p) => (
                  <option key={p.value} value={p.value}>{t(p.label)}</option>
                ))}
              </select>
            </div>
          </div>

          {/* 描述 */}
          <div>
            <label className="block text-xs font-medium text-muted-foreground mb-1">{t('描述')}</label>
            <textarea
              value={description}
              onChange={(e) => setDescription(e.target.value)}
              placeholder={t('描述这个模式的特点和使用场景')}
              rows={2}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-foreground/20 resize-none"
            />
          </div>

          {/* 模板 */}
          <div>
            <label className="block text-xs font-medium text-muted-foreground mb-1">{t('模板')}</label>
            <textarea
              value={template}
              onChange={(e) => setTemplate(e.target.value)}
              placeholder={t('模式模板，例如：[悬念问题] + [反转事实] + [引导互动]')}
              rows={3}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-foreground/20 resize-none"
            />
          </div>

          {/* 来源 */}
          <div>
            <label className="block text-xs font-medium text-muted-foreground mb-1">{t('来源')}</label>
            <select
              value={source}
              onChange={(e) => setSource(e.target.value as ViralPatternSource)}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-2 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-foreground/20"
            >
              {SOURCES.map((s) => (
                <option key={s.value} value={s.value}>{t(s.label)}</option>
              ))}
            </select>
          </div>
        </div>

        {/* 底部按钮 */}
        <div className="flex items-center justify-end gap-2 px-6 py-4 border-t border-border/40">
          <button
            type="button"
            onClick={() => onOpenChange(false)}
            className="rounded-md border border-border/60 bg-background px-4 py-2 text-sm text-foreground hover:bg-muted/40 transition-colors"
          >
            {t('取消')}
          </button>
          <button
            type="button"
            onClick={handleSubmit}
            disabled={!name.trim()}
            className="rounded-md bg-foreground px-4 py-2 text-sm font-medium text-background hover:bg-foreground/90 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {pattern ? t('保存') : t('创建')}
          </button>
        </div>
      </div>
    </div>
  )
}
