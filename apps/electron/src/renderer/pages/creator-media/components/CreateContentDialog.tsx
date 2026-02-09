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
import type { ContentType, ContentStatus, PipelineMode } from '@sprouty-ai/shared/db/types'

const CONTENT_TYPES: { value: ContentType; label: string }[] = [
  { value: 'image-text', label: '图文' },
  { value: 'video', label: '视频' },
  { value: 'short-video', label: '短视频' },
  { value: 'article', label: '文章' },
  { value: 'live', label: '直播' },
]

const PLATFORM_OPTIONS = [
  { value: 'xiaohongshu', label: '小红书' },
  { value: 'douyin', label: '抖音' },
  { value: 'bilibili', label: 'B站' },
  { value: 'wechat', label: '微信公众号' },
  { value: 'zhihu', label: '知乎' },
  { value: 'weibo', label: '微博' },
  { value: 'x', label: 'X (Twitter)' },
]

interface CreateContentDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  onCreateContent: (data: {
    title: string | null
    topic: string | null
    topic_source: null
    source_topic_id: null
    script_path: null
    status: ContentStatus
    content_type: ContentType | null
    target_platforms: string | null
    pipeline_mode: PipelineMode
    pipeline_state: null
    viral_pattern_id: null
    tags: null
    scheduled_at: string | null
    files: null
    metadata: null
  }) => Promise<unknown>
}

export function CreateContentDialog({ open, onOpenChange, onCreateContent }: CreateContentDialogProps) {
  const t = useT()
  const [saving, setSaving] = useState(false)

  const [title, setTitle] = useState('')
  const [topic, setTopic] = useState('')
  const [contentType, setContentType] = useState<ContentType>('image-text')
  const [selectedPlatforms, setSelectedPlatforms] = useState<string[]>([])
  const [scheduledAt, setScheduledAt] = useState('')

  const reset = () => {
    setTitle('')
    setTopic('')
    setContentType('image-text')
    setSelectedPlatforms([])
    setScheduledAt('')
    setSaving(false)
  }

  const handleClose = (v: boolean) => {
    if (!v) reset()
    onOpenChange(v)
  }

  const togglePlatform = (platform: string) => {
    setSelectedPlatforms((prev) =>
      prev.includes(platform) ? prev.filter((p) => p !== platform) : [...prev, platform]
    )
  }

  const handleCreate = async () => {
    if (!title.trim()) return
    setSaving(true)
    try {
      await onCreateContent({
        title: title.trim(),
        topic: topic.trim() || null,
        topic_source: null,
        source_topic_id: null,
        script_path: null,
        status: 'idea',
        content_type: contentType,
        target_platforms: selectedPlatforms.length > 0 ? selectedPlatforms.join(',') : null,
        pipeline_mode: 'manual',
        pipeline_state: null,
        viral_pattern_id: null,
        tags: null,
        scheduled_at: scheduledAt || null,
        files: null,
        metadata: null,
      })
      handleClose(false)
    } finally {
      setSaving(false)
    }
  }

  return (
    <Dialog open={open} onOpenChange={handleClose}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle>{t('新建内容')}</DialogTitle>
          <DialogDescription>{t('创建新的内容选题')}</DialogDescription>
        </DialogHeader>

        <div className="space-y-4">
          <div className="space-y-2">
            <Label>{t('标题')} *</Label>
            <Input
              value={title}
              onChange={(e) => setTitle(e.target.value)}
              placeholder={t('内容标题...')}
              autoFocus
            />
          </div>

          <div className="space-y-2">
            <Label>{t('选题描述')}</Label>
            <Textarea
              value={topic}
              onChange={(e) => setTopic(e.target.value)}
              placeholder={t('描述选题方向和要点...')}
              rows={2}
            />
          </div>

          <div className="space-y-2">
            <Label>{t('内容类型')}</Label>
            <Select value={contentType} onValueChange={(v) => setContentType(v as ContentType)}>
              <SelectTrigger>
                <SelectValue />
              </SelectTrigger>
              <SelectContent className="z-modal">
                {CONTENT_TYPES.map((ct) => (
                  <SelectItem key={ct.value} value={ct.value}>{t(ct.label)}</SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          <div className="space-y-2">
            <Label>{t('目标平台')}</Label>
            <div className="flex flex-wrap gap-2">
              {PLATFORM_OPTIONS.map((p) => (
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
          </div>

          <div className="space-y-2">
            <Label>{t('计划发布时间')}</Label>
            <Input
              type="datetime-local"
              value={scheduledAt}
              onChange={(e) => setScheduledAt(e.target.value)}
            />
          </div>
        </div>

        <DialogFooter>
          <Button variant="outline" onClick={() => handleClose(false)} disabled={saving}>
            {t('取消')}
          </Button>
          <Button onClick={handleCreate} disabled={!title.trim() || saving}>
            {saving ? t('创建中...') : t('创建')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
