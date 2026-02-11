import { useState, useCallback, useEffect, useMemo } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { PublishRecordTable } from './components/PublishRecordTable'
import type {
  Draft, Content, PublishRecord, PublishQueueItem,
} from '@sprouty-ai/shared/db/types'

type TabId = 'create' | 'drafts' | 'queue' | 'published'

/**
 * 发布工作台 — 创建发布、草稿管理、发布队列、已发布记录
 */
export default function PublishWorkbench() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const {
    projects, activeProject, contents, publishRecords, platformAccounts, loading,
    switchProject,
  } = useCreatorMedia()

  const [activeTab, setActiveTab] = useState<TabId>('create')
  const [drafts, setDrafts] = useState<Draft[]>([])
  const [queueItems, setQueueItems] = useState<PublishQueueItem[]>([])
  const [draftsLoading, setDraftsLoading] = useState(false)

  // 标签页配置
  const tabs: { id: TabId; label: string }[] = [
    { id: 'create', label: '创建发布' },
    { id: 'drafts', label: '草稿箱' },
    { id: 'queue', label: '发布队列' },
    { id: 'published', label: '已发布' },
  ]

  /** 加载草稿列表 */
  const loadDrafts = useCallback(async () => {
    if (!workspaceId || !activeProject) return
    setDraftsLoading(true)
    try {
      const list = await window.electronAPI.creatorMedia.drafts.list(workspaceId, activeProject.id)
      setDrafts(list)
    } catch {
      setDrafts([])
    } finally {
      setDraftsLoading(false)
    }
  }, [workspaceId, activeProject])

  /** 加载发布队列 */
  const loadQueue = useCallback(async () => {
    if (!workspaceId || !activeProject) return
    try {
      const allItems: PublishQueueItem[] = []
      for (const c of contents) {
        try {
          const items = await window.electronAPI.creatorMedia.publishQueue.list(workspaceId, c.id)
          allItems.push(...items)
        } catch {
          // 单条失败不影响整体
        }
      }
      setQueueItems(allItems)
    } catch {
      setQueueItems([])
    }
  }, [workspaceId, activeProject, contents])

  useEffect(() => {
    if (activeProject) {
      loadDrafts()
      loadQueue()
    } else {
      setDrafts([])
      setQueueItems([])
    }
  }, [activeProject, loadDrafts, loadQueue])

  /** 删除草稿 */
  const handleDeleteDraft = useCallback(async (id: string) => {
    if (!workspaceId) return
    try {
      await window.electronAPI.creatorMedia.drafts.delete(workspaceId, id)
      setDrafts(prev => prev.filter(d => d.id !== id))
    } catch {
      // 静默处理
    }
  }, [workspaceId])

  /** 待处理的队列项 */
  const pendingQueue = useMemo(() => {
    return queueItems.filter(q => q.status === 'pending' || q.status === 'processing')
  }, [queueItems])

  /** 已发布的内容 */
  const publishedContents = useMemo(() => {
    return contents.filter(c => c.status === 'published')
  }, [contents])

  if (!activeProject) {
    return (
      <div className="flex flex-col items-center justify-center h-full gap-4 text-muted-foreground">
        <p>{t('请先选择或创建一个项目')}</p>
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* 顶部：项目切换 + 标题 */}
      <div className="relative z-panel flex items-center justify-between px-6 py-4 border-b border-border/40">
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('发布工作台')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject.name}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <ProjectSwitcher
            projects={projects}
            activeProject={activeProject}
            onSwitch={switchProject}
          />
        </div>
      </div>

      {/* 标签页导航 */}
      <div className="flex items-center gap-0 px-6 border-b border-border/40">
        {tabs.map((tab) => (
          <button
            key={tab.id}
            type="button"
            onClick={() => setActiveTab(tab.id)}
            className={`px-4 py-2.5 text-sm transition-colors relative ${
              activeTab === tab.id
                ? 'text-foreground font-medium'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            {t(tab.label)}
            {activeTab === tab.id && (
              <span className="absolute bottom-0 left-0 right-0 h-0.5 bg-foreground rounded-full" />
            )}
            {/* 数量徽标 */}
            {tab.id === 'drafts' && drafts.length > 0 && (
              <span className="ml-1.5 inline-flex items-center justify-center min-w-[18px] h-[18px] px-1 text-[10px] font-medium rounded-full bg-muted text-muted-foreground">
                {drafts.length}
              </span>
            )}
            {tab.id === 'queue' && pendingQueue.length > 0 && (
              <span className="ml-1.5 inline-flex items-center justify-center min-w-[18px] h-[18px] px-1 text-[10px] font-medium rounded-full bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400">
                {pendingQueue.length}
              </span>
            )}
            {tab.id === 'published' && publishedContents.length > 0 && (
              <span className="ml-1.5 inline-flex items-center justify-center min-w-[18px] h-[18px] px-1 text-[10px] font-medium rounded-full bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400">
                {publishedContents.length}
              </span>
            )}
          </button>
        ))}
      </div>

      {/* 标签页内容 */}
      <div className="flex-1 overflow-auto px-6 py-4">
        {activeTab === 'create' && (
          <PublishCreateSection
            projectId={activeProject.id}
            platformAccounts={platformAccounts}
            contents={contents}
          />
        )}
        {activeTab === 'drafts' && (
          <DraftsSection
            drafts={drafts}
            loading={draftsLoading}
            onDelete={handleDeleteDraft}
            onRefresh={loadDrafts}
          />
        )}
        {activeTab === 'queue' && (
          <QueueSection items={queueItems} />
        )}
        {activeTab === 'published' && (
          <PublishRecordTable
            records={publishRecords}
            loading={loading}
          />
        )}
      </div>
    </div>
  )
}

// ============================================================
// 创建发布区域
// ============================================================

function PublishCreateSection({ projectId, platformAccounts, contents }: {
  projectId: string
  platformAccounts: any[]
  contents: Content[]
}) {
  const t = useT()
  const [title, setTitle] = useState('')
  const [content, setContent] = useState('')
  const [tags, setTags] = useState<string[]>([])
  const [tagInput, setTagInput] = useState('')
  const [selectedPlatforms, setSelectedPlatforms] = useState<string[]>([])

  const readyContents = contents.filter(c => c.status === 'reviewing' || c.status === 'scheduled')

  const platforms = [
    { id: 'xiaohongshu', name: '小红书', color: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400' },
    { id: 'douyin', name: '抖音', color: 'bg-gray-100 text-gray-700 dark:bg-gray-800/50 dark:text-gray-300' },
    { id: 'wechat', name: '公众号', color: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400' },
    { id: 'bilibili', name: 'B站', color: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400' },
    { id: 'weibo', name: '微博', color: 'bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-400' },
    { id: 'zhihu', name: '知乎', color: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400' },
    { id: 'x', name: 'X', color: 'bg-gray-100 text-gray-700 dark:bg-gray-800/50 dark:text-gray-300' },
  ]

  const togglePlatform = (id: string) => {
    setSelectedPlatforms((prev) =>
      prev.includes(id) ? prev.filter((p) => p !== id) : [...prev, id]
    )
  }

  const addTag = () => {
    const trimmed = tagInput.trim()
    if (trimmed && !tags.includes(trimmed)) {
      setTags([...tags, trimmed])
      setTagInput('')
    }
  }

  const removeTag = (tag: string) => {
    setTags(tags.filter((t) => t !== tag))
  }

  return (
    <div className="max-w-4xl mx-auto space-y-6">
      {/* 待发布内容提示 */}
      {readyContents.length > 0 && (
        <div>
          <h2 className="text-sm font-medium text-foreground mb-2">{t('待发布内容')}</h2>
          <div className="space-y-2">
            {readyContents.map(c => (
              <div key={c.id} className="flex items-center justify-between rounded-lg border border-border/60 bg-background/40 px-4 py-3">
                <div className="flex-1 min-w-0">
                  <p className="text-sm font-medium text-foreground truncate">{c.title || t('无标题')}</p>
                  <p className="text-xs text-muted-foreground mt-0.5">
                    {c.content_type} · {c.status === 'reviewing' ? t('审核中') : t('已排期')}
                    {c.scheduled_at && ` · ${c.scheduled_at.slice(0, 16)}`}
                  </p>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* 标题 */}
      <div>
        <label className="block text-sm font-medium mb-1.5">{t('标题')}</label>
        <input
          type="text"
          value={title}
          onChange={(e) => setTitle(e.target.value)}
          placeholder={t('输入内容标题...')}
          className="w-full px-3 py-2 rounded-md border border-input bg-background text-sm"
        />
      </div>

      {/* 正文 */}
      <div>
        <label className="block text-sm font-medium mb-1.5">{t('正文')}</label>
        <textarea
          value={content}
          onChange={(e) => setContent(e.target.value)}
          placeholder={t('输入正文内容...')}
          rows={10}
          className="w-full px-3 py-2 rounded-md border border-input bg-background text-sm resize-y"
        />
      </div>

      {/* 标签 */}
      <div>
        <label className="block text-sm font-medium mb-1.5">{t('标签')}</label>
        <div className="flex flex-wrap gap-2 mb-2">
          {tags.map((tag) => (
            <span
              key={tag}
              className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs"
            >
              #{tag}
              <button type="button" onClick={() => removeTag(tag)} className="hover:text-destructive">x</button>
            </span>
          ))}
        </div>
        <div className="flex gap-2">
          <input
            type="text"
            value={tagInput}
            onChange={(e) => setTagInput(e.target.value)}
            onKeyDown={(e) => e.key === 'Enter' && (e.preventDefault(), addTag())}
            placeholder={t('输入标签，回车添加')}
            className="flex-1 px-3 py-1.5 rounded-md border border-input bg-background text-sm"
          />
          <button
            type="button"
            onClick={addTag}
            className="px-3 py-1.5 rounded-md bg-secondary text-secondary-foreground text-sm"
          >
            {t('添加')}
          </button>
        </div>
      </div>

      {/* 平台选择 */}
      <div>
        <label className="block text-sm font-medium mb-1.5">{t('目标平台')}</label>
        <div className="flex flex-wrap gap-2">
          {platforms.map((p) => (
            <button
              key={p.id}
              type="button"
              onClick={() => togglePlatform(p.id)}
              className={`px-3 py-1.5 rounded-full text-xs font-medium transition-colors ${
                selectedPlatforms.includes(p.id)
                  ? p.color + ' ring-2 ring-primary/50'
                  : 'bg-muted text-muted-foreground'
              }`}
            >
              {p.name}
            </button>
          ))}
        </div>
        {/* 已登录账号提示 */}
        {platformAccounts.length > 0 && (
          <p className="mt-2 text-xs text-muted-foreground">
            {platformAccounts.filter((a: any) => a.auth_status === 'logged_in').length} / {platformAccounts.length} {t('个账号已登录')}
          </p>
        )}
      </div>

      {/* 操作按钮 */}
      <div className="flex gap-3 pt-4 border-t border-border">
        <button type="button" className="px-4 py-2 rounded-md bg-secondary text-secondary-foreground text-sm">
          {t('保存草稿')}
        </button>
        <button type="button" className="px-4 py-2 rounded-md bg-primary text-primary-foreground text-sm">
          {t('发布')}
        </button>
      </div>
    </div>
  )
}

// ============================================================
// 草稿箱区域
// ============================================================

function DraftsSection({ drafts, loading, onDelete, onRefresh }: {
  drafts: Draft[]
  loading: boolean
  onDelete: (id: string) => void
  onRefresh: () => void
}) {
  const t = useT()

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  if (drafts.length === 0) {
    return (
      <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-8 text-center">
        <p className="text-sm text-muted-foreground">{t('草稿箱为空')}</p>
        <p className="mt-1 text-xs text-muted-foreground/60">{t('通过 AI 创作的草稿会保存在这里')}</p>
      </div>
    )
  }

  return (
    <div className="space-y-2">
      {drafts.map(d => (
        <div key={d.id} className="flex items-center justify-between rounded-lg border border-border/60 bg-background/40 px-4 py-3">
          <div className="flex-1 min-w-0">
            <p className="text-sm font-medium text-foreground truncate">{d.title || t('无标题')}</p>
            <p className="text-xs text-muted-foreground mt-0.5">
              {d.content_type || t('文本')}
              {d.target_platforms && ` · ${d.target_platforms}`}
              {d.updated_at && ` · ${t('更新于')} ${d.updated_at.slice(0, 16)}`}
            </p>
          </div>
          <button
            type="button"
            onClick={() => onDelete(d.id)}
            className="ml-4 text-xs text-muted-foreground hover:text-destructive transition-colors"
          >
            {t('删除')}
          </button>
        </div>
      ))}
    </div>
  )
}

// ============================================================
// 发布队列区域
// ============================================================

function QueueSection({ items }: { items: PublishQueueItem[] }) {
  const t = useT()

  const statusLabel: Record<string, string> = {
    pending: '等待中',
    processing: '发布中',
    success: '成功',
    failed: '失败',
    cancelled: '已取消',
  }

  const statusColor: Record<string, string> = {
    pending: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
    processing: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400',
    success: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
    failed: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
    cancelled: 'bg-muted text-muted-foreground',
  }

  if (items.length === 0) {
    return (
      <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-8 text-center">
        <p className="text-sm text-muted-foreground">{t('发布队列为空')}</p>
        <p className="mt-1 text-xs text-muted-foreground/60">{t('排期发布的内容会出现在这里')}</p>
      </div>
    )
  }

  return (
    <div className="space-y-2">
      {items.map(item => (
        <div key={item.id} className="flex items-center justify-between rounded-lg border border-border/60 bg-background/40 px-4 py-3">
          <div className="flex items-center gap-3 min-w-0">
            <div className="min-w-0">
              <p className="text-sm text-foreground truncate">
                {item.content_id.slice(0, 8)}...
              </p>
              <p className="text-xs text-muted-foreground mt-0.5">
                {item.scheduled_at ? `${t('排期')} ${item.scheduled_at.slice(0, 16)}` : t('立即发布')}
              </p>
            </div>
          </div>
          <span className={`text-xs px-2 py-0.5 rounded-full ${statusColor[item.status] || 'bg-muted text-muted-foreground'}`}>
            {t(statusLabel[item.status] || item.status)}
          </span>
        </div>
      ))}
    </div>
  )
}
