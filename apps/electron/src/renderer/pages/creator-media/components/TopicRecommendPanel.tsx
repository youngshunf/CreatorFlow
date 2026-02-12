import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { navigate, routes } from '@/lib/navigate'
import { Info_Markdown } from '@/components/info'
import type { RecommendedTopic } from '@sprouty-ai/shared/db/types'

/** 状态标签 */
const STATUS_MAP: Record<number, { label: string; color: string }> = {
  0: { label: '待选', color: 'bg-blue-600 text-white dark:bg-blue-500 dark:text-white' },
  1: { label: '已采纳', color: 'bg-green-600 text-white dark:bg-green-500 dark:text-white' },
  2: { label: '已忽略', color: 'bg-gray-500 text-white dark:bg-gray-500 dark:text-white' },
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
  const [adoptedTopic, setAdoptedTopic] = useState<RecommendedTopic | null>(null)
  const [expandedId, setExpandedId] = useState<string | null>(null)
  const [mdContent, setMdContent] = useState<string | null>(null)
  const [mdLoading, setMdLoading] = useState(false)

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
      const result = await window.electronAPI.creatorMedia.topics.adopt(workspaceId, topicId, projectId)
      await loadTopics()
      if (result?.topic) {
        setAdoptedTopic(result.topic)
      }
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

  /** 开始创作 — 跳转到脚本创作技能 */
  const handleStartCreate = useCallback((topic: RecommendedTopic) => {
    if (!workspaceId) return
    const mdPath = topic.md_file_path ? ` 选题详情文件路径: ${topic.md_file_path}` : ''
    navigate(routes.action.newSession({
      input: `[skill:${workspaceId}:script-create] 为选题「${topic.title}」创作脚本。${mdPath}`,
      send: true,
    }))
    setAdoptedTopic(null)
  }, [workspaceId])

  /** 展开/收起 Markdown 预览 */
  const handleToggleExpand = useCallback(async (topic: RecommendedTopic) => {
    if (expandedId === topic.id) {
      setExpandedId(null)
      setMdContent(null)
      return
    }
    setExpandedId(topic.id)
    setMdContent(null)
    if (!topic.md_file_path || !workspaceId) {
      return
    }
    setMdLoading(true)
    try {
      const content = await window.electronAPI.creatorMedia.topics.readMd(workspaceId, topic.md_file_path)
      setMdContent(content)
    } catch {
      setMdContent(null)
    } finally {
      setMdLoading(false)
    }
  }, [expandedId, workspaceId])

  return (
    <div className="flex flex-col rounded-lg border border-border/60 bg-background/40">
      {/* 头部 */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-border/40">
        <h2 className="text-sm font-medium text-foreground">{t('选题推荐')}</h2>
        <div className="flex items-center gap-2">
          <button
            type="button"
            onClick={() => navigate(routes.action.newSession({
              input: `[skill:${workspaceId}:topic-generator] 为当前活跃项目生成选题推荐`,
              send: true,
            }))}
            className="inline-flex items-center gap-1 rounded-md text-xs text-muted-foreground hover:text-foreground transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456Z" />
            </svg>
            {t('AI 生成')}
          </button>
          {/* 筛选 tab */}
          <div className="flex items-center gap-0.5 rounded-md bg-muted/40 p-0.5">
          {FILTER_TABS.map((tab) => (
            <button
              key={tab.label}
              type="button"
              onClick={() => setFilterStatus(tab.status)}
              className={`rounded px-2 py-0.5 text-xs transition-colors ${
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
      </div>

      {/* 内容 */}
      <div className="px-4 py-3">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <p className="text-xs text-muted-foreground">{t('加载中...')}</p>
          </div>
        ) : topics.length === 0 ? (
          <div className="flex flex-col items-center justify-center py-8 text-center">
            <p className="text-sm text-muted-foreground">{t('暂无选题推荐')}</p>
            <p className="mt-1 text-xs text-muted-foreground/70">{t('点击"AI 生成"按钮，使用 topic-generator 技能生成选题')}</p>
          </div>
        ) : (
          <>
            {/* 卡片网格：一行3个 */}
            <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-3">
              {topics.map((topic) => (
                <TopicCard
                  key={topic.id}
                  topic={topic}
                  adopting={adoptingId === topic.id}
                  isExpanded={expandedId === topic.id}
                  onAdopt={handleAdopt}
                  onIgnore={handleIgnore}
                  onToggleExpand={() => handleToggleExpand(topic)}
                />
              ))}
            </div>

            {/* Markdown 预览面板 — 展开时显示在卡片网格下方 */}
            {expandedId && (
              <div className="mt-3 rounded-lg border border-border/40 bg-background/60 px-4 py-3">
                {mdLoading ? (
                  <p className="text-xs text-muted-foreground py-2">{t('加载详情...')}</p>
                ) : mdContent ? (
                  <Info_Markdown maxHeight={400} fullscreen>
                    {mdContent}
                  </Info_Markdown>
                ) : (
                  <p className="text-xs text-muted-foreground/70 py-2">{t('暂无详情')}</p>
                )}
              </div>
            )}
          </>
        )}
      </div>

      {/* 采纳成功提示 */}
      {adoptedTopic && (
        <div className="flex items-center justify-between px-4 py-3 border-t border-border/40 bg-green-500/5">
          <div className="flex items-center gap-2">
            <svg className="w-4 h-4 text-green-600 dark:text-green-400" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="m4.5 12.75 6 6 9-13.5" />
            </svg>
            <span className="text-sm text-foreground">
              {t('已创建内容')}：{adoptedTopic.title}
            </span>
          </div>
          <div className="flex items-center gap-2">
            <button
              type="button"
              onClick={() => handleStartCreate(adoptedTopic)}
              className="inline-flex items-center gap-1 rounded-md bg-green-600 px-3 py-1.5 text-xs font-medium text-white hover:bg-green-700 transition-colors"
            >
              <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                <path strokeLinecap="round" strokeLinejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L10.582 16.07a4.5 4.5 0 0 1-1.897 1.13L6 18l.8-2.685a4.5 4.5 0 0 1 1.13-1.897l8.932-8.931Zm0 0L19.5 7.125M18 14v4.75A2.25 2.25 0 0 1 15.75 21H5.25A2.25 2.25 0 0 1 3 18.75V8.25A2.25 2.25 0 0 1 5.25 6H10" />
              </svg>
              {t('开始创作')}
            </button>
            <button
              type="button"
              onClick={() => setAdoptedTopic(null)}
              className="rounded-md px-2 py-1.5 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors"
            >
              {t('稍后')}
            </button>
          </div>
        </div>
      )}
    </div>
  )
}

/** 单条选题卡片 */
function TopicCard({ topic, adopting, isExpanded, onAdopt, onIgnore, onToggleExpand }: {
  topic: RecommendedTopic
  adopting: boolean
  isExpanded: boolean
  onAdopt: (id: string) => void
  onIgnore: (id: string) => void
  onToggleExpand: () => void
}) {
  const t = useT()
  const statusInfo = STATUS_MAP[topic.status] ?? STATUS_MAP[0]

  return (
    <div
      className={`flex flex-col rounded-lg border bg-background/60 px-3 py-3 cursor-pointer transition-colors hover:bg-muted/20 ${
        isExpanded ? 'border-foreground/30 ring-1 ring-foreground/10' : 'border-border/40'
      }`}
      onClick={onToggleExpand}
    >
      {/* 标题 + 状态 */}
      <div className="flex items-start justify-between gap-2 mb-2">
        <p className="text-sm font-medium text-foreground leading-snug line-clamp-2 flex-1">{topic.title}</p>
        <span className={`shrink-0 inline-flex items-center rounded-full px-1.5 py-0.5 text-xs font-medium ${statusInfo.color}`}>
          {t(statusInfo.label)}
        </span>
      </div>

      {/* 推荐理由 */}
      {topic.reason && (
        <p className="text-xs text-muted-foreground leading-relaxed line-clamp-2 mb-2">{topic.reason}</p>
      )}

      {/* 分数指标 */}
      <div className="flex items-center gap-3 text-xs text-muted-foreground mb-2">
        <span className="inline-flex items-center gap-1">
          <svg className="w-3 h-3 text-orange-400" fill="currentColor" viewBox="0 0 20 20">
            <path d="M10 15l-5.878 3.09 1.123-6.545L.489 6.91l6.572-.955L10 0l2.939 5.955 6.572.955-4.756 4.635 1.123 6.545z" />
          </svg>
          {topic.potential_score}
        </span>
        <span className="inline-flex items-center gap-1">
          <svg className="w-3 h-3 text-red-400" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M15.362 5.214A8.252 8.252 0 0 1 12 21 8.25 8.25 0 0 1 6.038 7.047 8.287 8.287 0 0 0 9 9.601a8.983 8.983 0 0 1 3.361-6.867 8.21 8.21 0 0 0 3 2.48Z" />
          </svg>
          {topic.heat_index}
        </span>
        {topic.batch_date && (
          <span className="ml-auto">{topic.batch_date}</span>
        )}
      </div>

      {/* 操作按钮 */}
      {topic.status === 0 && (
        <div className="flex items-center gap-2 pt-1 border-t border-border/30">
          <button
            type="button"
            onClick={(e) => { e.stopPropagation(); onAdopt(topic.id) }}
            disabled={adopting}
            className="inline-flex items-center gap-1 rounded px-2 py-1 text-xs font-medium text-green-600 hover:bg-green-50 dark:text-green-400 dark:hover:bg-green-900/20 transition-colors disabled:opacity-50"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="m4.5 12.75 6 6 9-13.5" />
            </svg>
            {adopting ? t('采纳中...') : t('采纳')}
          </button>
          <button
            type="button"
            onClick={(e) => { e.stopPropagation(); onIgnore(topic.id) }}
            className="inline-flex items-center gap-1 rounded px-2 py-1 text-xs font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M6 18 18 6M6 6l12 12" />
            </svg>
            {t('忽略')}
          </button>
        </div>
      )}

      {/* 已采纳提示 */}
      {topic.status === 1 && topic.content_id && (
        <p className="text-xs text-green-600 dark:text-green-400 pt-1 border-t border-border/30">
          {t('已创建内容')}
        </p>
      )}
    </div>
  )
}
