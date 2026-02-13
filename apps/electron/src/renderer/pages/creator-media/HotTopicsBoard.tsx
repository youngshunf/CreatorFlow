import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { TopicRecommendPanel } from './components/TopicRecommendPanel'
import { HotTopicsPanel } from './components/HotTopicsPanel'
import type { HotTopic } from '@sprouty-ai/shared/db/types'
import { navigate, routes } from '@/lib/navigate'

/** 获取当天日期 */
function getTodayDate(): string {
  return new Date().toISOString().split('T')[0]!
}

/** 格式化热度值 */
function formatHeat(score: number | null): string {
  if (score == null) return '-'
  if (score >= 10000) return `${(score / 10000).toFixed(1)}万`
  return String(score)
}

/**
 * 热点看板 — 推荐选题 + AI 侦察热点 + 全网热榜
 */
export default function HotTopicsBoard() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const {
    projects, activeProject, loading, switchProject,
  } = useCreatorMedia()

  // AI 侦察热点
  const [scoutTopics, setScoutTopics] = useState<HotTopic[]>([])
  const [scoutLoading, setScoutLoading] = useState(false)

  /** 加载 AI 侦察热点 */
  const loadScoutTopics = useCallback(async () => {
    if (!workspaceId) return
    setScoutLoading(true)
    try {
      const list = await window.electronAPI.creatorMedia.hotTopics.list(workspaceId, {
        fetchSource: 'ai-scout',
        batchDate: getTodayDate(),
      })
      setScoutTopics(list)
    } catch {
      setScoutTopics([])
    } finally {
      setScoutLoading(false)
    }
  }, [workspaceId])

  useEffect(() => { loadScoutTopics() }, [loadScoutTopics])

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* 头部 */}
      <div className="relative z-panel flex items-center justify-between px-6 py-4 border-b border-border/40">
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('热点看板')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('实时热点流与智能选题')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-6">
        {/* AI 侦察热点 */}
        <ScoutPanel
          workspaceId={workspaceId}
          topics={scoutTopics}
          loading={scoutLoading}
        />

        {/* 全网热榜 */}
        <HotTopicsPanel />
      </div>
    </div>
  )
}

/** AI 侦察热点面板 — 与 HotTopicsPanel 同风格 */
function ScoutPanel({ workspaceId, topics, loading }: {
  workspaceId: string
  topics: HotTopic[]
  loading: boolean
}) {
  const t = useT()

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
        <h2 className="text-base font-medium text-foreground">{t('AI 侦察热点')}</h2>
        <button
          type="button"
          onClick={() => navigate(routes.action.newSession({
            input: `[skill:${workspaceId}:hot-topic-scout] 为当前活跃项目侦察热点`,
            send: true,
          }))}
          className="inline-flex items-center gap-1 rounded-md text-sm text-muted-foreground hover:text-foreground transition-colors"
        >
          <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456Z" />
          </svg>
          {t('AI 侦察')}
        </button>
      </div>

      {/* 内容 */}
      <div className="flex-1 overflow-auto max-h-[840px] px-4 py-2">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
          </div>
        ) : topics.length === 0 ? (
          <div className="flex flex-col items-center justify-center py-8 text-center">
            <p className="text-sm text-muted-foreground">{t('暂无 AI 侦察热点')}</p>
            <p className="mt-1 text-xs text-muted-foreground">{t('点击"AI 侦察"按钮，基于账号画像发现相关热点')}</p>
          </div>
        ) : Object.keys(grouped).length > 1 ? (
          /* 多平台：网格布局，按平台分组显示 */
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
            {Object.entries(grouped).map(([platformId, items]) => (
              <ScoutPlatformCard
                key={platformId}
                platformId={platformId}
                platformName={t(items[0]?.platform_name ?? platformId)}
                items={items}
              />
            ))}
          </div>
        ) : (
          /* 单平台或全部同源 */
          <div className="space-y-0.5">
            {topics.map((item) => (
              <ScoutTopicRow key={item.id} item={item} />
            ))}
          </div>
        )}
      </div>
    </div>
  )
}

/** AI 侦察平台卡片 */
function ScoutPlatformCard({ platformId, platformName, items }: { platformId: string; platformName: string; items: HotTopic[] }) {
  return (
    <div className="rounded-lg border border-border/60 bg-background/40 overflow-hidden">
      {/* 平台标题 */}
      <div className="px-3 py-2 border-b border-border/40 bg-muted/20">
        <h3 className="text-sm font-medium text-foreground">{platformName}</h3>
      </div>
      {/* 热榜列表 */}
      <div className="p-2 space-y-0.5">
        {items.slice(0, 10).map((item) => (
          <CompactScoutTopicRow key={item.id} item={item} />
        ))}
      </div>
    </div>
  )
}

/** 紧凑型 AI 侦察热点行（用于卡片内） */
function CompactScoutTopicRow({ item }: { item: HotTopic }) {
  const displayTitle = item.title || item.platform_name || '-'

  const inner = (
    <>
      {/* 排名 */}
      <span className={`w-5 text-center text-sm font-semibold shrink-0 ${
        item.rank != null && item.rank <= 3 ? 'text-orange-500' : 'text-muted-foreground'
      }`}>
        {item.rank ?? '-'}
      </span>
      {/* 标题 */}
      <span className="flex-1 text-sm text-foreground truncate" title={displayTitle}>
        {displayTitle}
      </span>
      {/* 热度 */}
      {item.heat_score != null && (
        <span className="shrink-0 text-xs text-foreground/60">
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
        className="flex items-center gap-1.5 rounded px-1.5 py-1 hover:bg-muted/30 transition-colors cursor-pointer no-underline"
      >
        {inner}
      </a>
    )
  }

  return (
    <div className="flex items-center gap-1.5 rounded px-1.5 py-1 hover:bg-muted/30 transition-colors cursor-default">
      {inner}
    </div>
  )
}

/** AI 侦察热点行 — 与 HotTopicsPanel.TopicRow 同风格 */
function ScoutTopicRow({ item }: { item: HotTopic }) {
  const displayTitle = item.title || item.platform_name || '-'

  const inner = (
    <>
      <span className={`w-5 text-center text-sm font-semibold shrink-0 ${
        item.rank != null && item.rank <= 3 ? 'text-orange-500' : 'text-muted-foreground'
      }`}>
        {item.rank ?? '-'}
      </span>
      <span className="flex-1 text-sm text-foreground truncate" title={displayTitle}>
        {displayTitle}
      </span>
      {item.heat_score != null && (
        <span className="shrink-0 text-xs text-foreground/60">
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
