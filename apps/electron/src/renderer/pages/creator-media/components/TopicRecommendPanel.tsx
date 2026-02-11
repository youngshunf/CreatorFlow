import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import type { RecommendedTopic } from '@sprouty-ai/shared/db/types'

/** 状态标签 */
const STATUS_MAP: Record<number, { label: string; color: string }> = {
  0: { label: '待选', color: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400' },
  1: { label: '已采纳', color: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400' },
  2: { label: '已忽略', color: 'bg-gray-100 text-gray-600 dark:bg-gray-800/30 dark:text-gray-400' },
}

/** 筛选 tab */
const FILTER_TABS = [
  { status: undefined, label: '全部' },
  { status: 0, label: '待选' },
  { status: 1, label: '已采纳' },
] as const

export function TopicRecommendPanel({ projectId }: { projectId: string }) {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''

  const [topics, setTopics] = useState<RecommendedTopic[]>([])
  const [loading, setLoading] = useState(false)
  const [filterStatus, setFilterStatus] = useState<number | undefined>(0)
  const [adoptingId, setAdoptingId] = useState<string | null>(null)

  /** 加载选题列表 */
  const loadTopics = useCallback(async () => {
    if (!workspaceId || !projectId) return
    setLoading(true)
    try {
      const filters: { status?: number; limit?: number } = { limit: 50 }
      if (filterStatus != null) filters.status = filterStatus
      const list = await window.electronAPI.creatorMedia.topics.list(workspaceId, projectId, filters)
      setTopics(list)
    } catch {
      setTopics([])
    } finally {
      setLoading(false)
    }
  }, [workspaceId, projectId, filterStatus])

  useEffect(() => { loadTopics() }, [loadTopics])

  /** 采纳选题 */
  const handleAdopt = useCallback(async (topicId: string) => {
    if (!workspaceId || adoptingId) return
    setAdoptingId(topicId)
    try {
      await window.electronAPI.creatorMedia.topics.adopt(workspaceId, topicId, projectId)
      await loadTopics()
    } finally {
      setAdoptingId(null)
    }
  }, [workspaceId, projectId, adoptingId, loadTopics])

  /** 忽略选题 */
  const handleIgnore = useCallback(async (topicId: string) => {
    if (!workspaceId) return
    await window.electronAPI.creatorMedia.topics.ignore(workspaceId, topicId)
    await loadTopics()
  }, [workspaceId, loadTopics])

  return (
    <div className="flex flex-col rounded-lg border border-border/60 bg-background/40">
      {/* 头部 */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-border/40">
        <h2 className="text-sm font-medium text-foreground">{t('选题推荐')}</h2>
        {/* 筛选 tab */}
        <div className="flex items-center gap-0.5 rounded-md bg-muted/40 p-0.5">
          {FILTER_TABS.map((tab) => (
            <button
              key={tab.label}
              type="button"
              onClick={() => setFilterStatus(tab.status)}
              className={`rounded px-2 py-0.5 text-[11px] transition-colors ${
                filterStatus === tab.status
                  ? 'bg-background text-foreground shadow-sm'
                  : 'text-muted-foreground hover:text-foreground'
              }`}
            >
              {t(tab.label)}
            </button>
          ))}
        </div>
      </div>

      {/* 内容 */}
      <div className="flex-1 overflow-auto max-h-[420px] px-4 py-2">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <p className="text-xs text-muted-foreground">{t('加载中...')}</p>
          </div>
        ) : topics.length === 0 ? (
          <div className="flex flex-col items-center justify-center py-8 text-center">
            <p className="text-sm text-muted-foreground">{t('暂无选题推荐')}</p>
            <p className="mt-1 text-xs text-muted-foreground/70">{t('请先拉取热榜数据，系统将自动生成选题')}</p>
          </div>
        ) : (
          <div className="space-y-2">
            {topics.map((topic) => (
              <TopicCard
                key={topic.id}
                topic={topic}
                adopting={adoptingId === topic.id}
                onAdopt={handleAdopt}
                onIgnore={handleIgnore}
              />
            ))}
          </div>
        )}
      </div>
    </div>
  )
}

/** 单条选题卡片 */
function TopicCard({ topic, adopting, onAdopt, onIgnore }: {
  topic: RecommendedTopic
  adopting: boolean
  onAdopt: (id: string) => void
  onIgnore: (id: string) => void
}) {
  const t = useT()
  const statusInfo = STATUS_MAP[topic.status] ?? STATUS_MAP[0]

  // 解析关键词
  let keywords: string[] = []
  if (topic.keywords) {
    try { keywords = JSON.parse(topic.keywords) } catch { /* ignore */ }
  }

  return (
    <div className="rounded-lg border border-border/40 bg-background/60 px-3 py-2.5 space-y-1.5">
      {/* 标题行 */}
      <div className="flex items-start justify-between gap-2">
        <p className="text-xs font-medium text-foreground leading-relaxed flex-1">{topic.title}</p>
        <span className={`shrink-0 inline-flex items-center rounded-full px-1.5 py-0.5 text-[10px] font-medium ${statusInfo.color}`}>
          {t(statusInfo.label)}
        </span>
      </div>

      {/* 分数 + 关键词 */}
      <div className="flex items-center gap-3 text-[10px] text-muted-foreground">
        <span>{t('潜力')} {topic.potential_score}</span>
        <span>{t('热度')} {topic.heat_index}</span>
        {keywords.length > 0 && (
          <div className="flex items-center gap-1 flex-1 min-w-0 overflow-hidden">
            {keywords.slice(0, 3).map((kw, i) => (
              <span key={i} className="inline-flex rounded bg-muted/60 px-1.5 py-0.5 text-[10px] text-muted-foreground truncate">
                {kw}
              </span>
            ))}
          </div>
        )}
      </div>

      {/* 推荐理由 */}
      {topic.reason && (
        <p className="text-[11px] text-muted-foreground/80 leading-relaxed line-clamp-2">{topic.reason}</p>
      )}

      {/* 操作按钮 */}
      {topic.status === 0 && (
        <div className="flex items-center gap-2 pt-0.5">
          <button
            type="button"
            onClick={() => onAdopt(topic.id)}
            disabled={adopting}
            className="inline-flex items-center gap-0.5 rounded px-2 py-0.5 text-[10px] font-medium text-green-600 hover:bg-green-50 dark:text-green-400 dark:hover:bg-green-900/20 transition-colors disabled:opacity-50"
          >
            <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="m4.5 12.75 6 6 9-13.5" />
            </svg>
            {adopting ? t('采纳中...') : t('采纳')}
          </button>
          <button
            type="button"
            onClick={() => onIgnore(topic.id)}
            className="inline-flex items-center gap-0.5 rounded px-2 py-0.5 text-[10px] font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors"
          >
            <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M6 18 18 6M6 6l12 12" />
            </svg>
            {t('忽略')}
          </button>
        </div>
      )}

      {/* 已采纳：显示关联内容 */}
      {topic.status === 1 && topic.content_id && (
        <p className="text-[10px] text-green-600 dark:text-green-400">
          {t('已创建内容')}
        </p>
      )}
    </div>
  )
}
