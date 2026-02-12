import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import type { HotTopic } from '@sprouty-ai/shared/db/types'
import { HOT_TOPIC_PLATFORMS } from '@sprouty-ai/shared/services/hot-topic-service'

/** 平台筛选列表（含"全部"选项） */
const PLATFORMS = [
  { id: 'all', name: '全部平台' },
  ...HOT_TOPIC_PLATFORMS,
]

/** 格式化热度值 */
function formatHeat(score: number | null): string {
  if (score == null) return '-'
  if (score >= 10000) return `${(score / 10000).toFixed(1)}万`
  return String(score)
}

export function HotTopicsPanel() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''

  const [topics, setTopics] = useState<HotTopic[]>([])
  const [loading, setLoading] = useState(false)
  const [fetching, setFetching] = useState(false)
  const [selectedPlatform, setSelectedPlatform] = useState('all')
  const [batchInfo, setBatchInfo] = useState<{ batchDate: string; fetchedAt: string; count: number } | null>(null)

  /** 加载热榜列表 */
  const loadTopics = useCallback(async () => {
    if (!workspaceId) return
    setLoading(true)
    try {
      const filters: { platformId?: string; limit?: number } = { limit: 100 }
      if (selectedPlatform !== 'all') filters.platformId = selectedPlatform
      const [list, batch] = await Promise.all([
        window.electronAPI.creatorMedia.hotTopics.list(workspaceId, filters),
        window.electronAPI.creatorMedia.hotTopics.getLatestBatch(workspaceId),
      ])
      setTopics(list)
      setBatchInfo(batch)
    } catch {
      setTopics([])
    } finally {
      setLoading(false)
    }
  }, [workspaceId, selectedPlatform])

  useEffect(() => { loadTopics() }, [loadTopics])

  /** 手动拉取热榜 */
  const handleFetch = useCallback(async () => {
    if (!workspaceId || fetching) return
    setFetching(true)
    try {
      await window.electronAPI.creatorMedia.hotTopics.fetch(workspaceId)
      await loadTopics()
    } finally {
      setFetching(false)
    }
  }, [workspaceId, fetching, loadTopics])

  /** 按平台分组 */
  const grouped = topics.reduce<Record<string, HotTopic[]>>((acc, item) => {
    const key = item.platform_id
    if (!acc[key]) acc[key] = []
    acc[key].push(item)
    return acc
  }, {})

  return (
    <div className="flex flex-col rounded-lg border border-border/60 bg-background/40">
      {/* 头部 */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-border/40">
        <div className="flex items-center gap-2">
          <h2 className="text-sm font-medium text-foreground">{t('热榜')}</h2>
          {batchInfo && (
            <span className="text-xs text-muted-foreground">
              {new Date(batchInfo.fetchedAt).toLocaleString('zh-CN', { month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit' })}
            </span>
          )}
        </div>
        <div className="flex items-center gap-2">
          {/* 平台筛选 */}
          <select
            value={selectedPlatform}
            onChange={(e) => setSelectedPlatform(e.target.value)}
            className="h-6 rounded border border-border/60 bg-background px-1.5 text-xs text-foreground outline-none focus:border-foreground/30"
          >
            {PLATFORMS.map((p) => (
              <option key={p.id} value={p.id}>{t(p.name)}</option>
            ))}
          </select>
          {/* 刷新按钮 */}
          <button
            type="button"
            onClick={handleFetch}
            disabled={fetching}
            className="inline-flex items-center gap-1 rounded-md text-xs text-muted-foreground hover:text-foreground transition-colors disabled:opacity-50"
          >
            <svg className={`w-3.5 h-3.5 ${fetching ? 'animate-spin' : ''}`} fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0 3.181 3.183a8.25 8.25 0 0 0 13.803-3.7M4.031 9.865a8.25 8.25 0 0 1 13.803-3.7l3.181 3.182" />
            </svg>
            {t('刷新')}
          </button>
        </div>
      </div>

      {/* 内容 */}
      <div className="flex-1 overflow-auto max-h-[840px] px-4 py-2">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <p className="text-xs text-muted-foreground">{t('加载中...')}</p>
          </div>
        ) : topics.length === 0 ? (
          <div className="flex flex-col items-center justify-center py-8 text-center">
            <p className="text-sm text-muted-foreground">{t('暂无热榜数据')}</p>
            <p className="mt-1 text-xs text-muted-foreground/70">{t('点击刷新按钮拉取最新热榜')}</p>
          </div>
        ) : selectedPlatform === 'all' ? (
          /* 全部平台：按平台分组显示 */
          Object.entries(grouped).map(([platformId, items]) => (
            <div key={platformId} className="mb-3 last:mb-0">
              <p className="text-xs font-medium text-muted-foreground mb-1.5">
                {t(PLATFORMS.find(p => p.id === platformId)?.name ?? items[0]?.platform_name ?? platformId)}
              </p>
              <div className="space-y-0.5">
                {items.slice(0, 10).map((item) => (
                  <TopicRow key={item.id} item={item} />
                ))}
              </div>
            </div>
          ))
        ) : (
          /* 单平台列表 */
          <div className="space-y-0.5">
            {topics.map((item) => (
              <TopicRow key={item.id} item={item} />
            ))}
          </div>
        )}
      </div>
    </div>
  )
}

/** 单条热榜行 */
function TopicRow({ item }: { item: HotTopic }) {
  const displayTitle = item.title || item.platform_name || '-'

  const inner = (
    <>
      {/* 排名 */}
      <span className={`w-5 text-center text-xs font-semibold shrink-0 ${
        item.rank != null && item.rank <= 3 ? 'text-orange-500' : 'text-muted-foreground/60'
      }`}>
        {item.rank ?? '-'}
      </span>
      {/* 标题 */}
      <span className="flex-1 text-sm text-foreground truncate" title={displayTitle}>
        {displayTitle}
      </span>
      {/* 热度 */}
      {item.heat_score != null && (
        <span className="shrink-0 text-xs text-muted-foreground/70">
          {formatHeat(item.heat_score)}
        </span>
      )}
    </>
  )

  if (item.url) {
    return (
      <a
        href={item.url}
        target="_blank"
        rel="noopener noreferrer"
        className="flex items-center gap-2 rounded px-2 py-1.5 hover:bg-muted/30 transition-colors cursor-pointer no-underline"
      >
        {inner}
      </a>
    )
  }

  return (
    <div className="flex items-center gap-2 rounded px-2 py-1.5 hover:bg-muted/30 transition-colors cursor-default">
      {inner}
    </div>
  )
}
