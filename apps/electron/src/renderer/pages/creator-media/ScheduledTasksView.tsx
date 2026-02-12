/**
 * 定时任务管理视图（Hooks 卡片视图 + 执行记录）
 *
 * 直接读取并展示工作区 hooks.json 的内容，以卡片视图呈现
 * 同时支持查看 events.jsonl 中的 hook 执行记录
 */

import { useState, useEffect, useCallback, useRef } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { EditPopover, EditButton, getEditConfig } from '@/components/ui/EditPopover'
import type { HookEventRecord } from '../../../shared/types'

/** hooks.json 中单个 hook 定义 */
interface HookDefinition {
  type: 'command' | 'prompt'
  command?: string
  prompt?: string
  timeout?: number
}

/** hooks.json 中的 matcher 条目 */
interface HookMatcher {
  matcher?: string
  cron?: string
  timezone?: string
  permissionMode?: string
  labels?: string[]
  enabled?: boolean
  hooks: HookDefinition[]
}

/** hooks.json 顶层结构 */
interface HooksConfig {
  hooks: Record<string, HookMatcher[]>
}

/** 扁平化后的卡片数据 */
interface HookCard {
  eventName: string
  matcherIndex: number
  matcher: HookMatcher
}

/** 事件名称的可读标签 */
const EVENT_LABELS: Record<string, string> = {
  SchedulerTick: '定时调度',
  UserPromptSubmit: '用户提交',
  SessionStart: '会话开始',
  SessionEnd: '会话结束',
  LabelAdd: '标签添加',
  LabelRemove: '标签移除',
  PreToolUse: '工具调用前',
  PostToolUse: '工具调用后',
  Stop: '会话停止',
}

/** 事件名称的颜色配置 */
const EVENT_COLORS: Record<string, string> = {
  SchedulerTick: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400',
  UserPromptSubmit: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  SessionStart: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400',
  SessionEnd: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400',
  LabelAdd: 'bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-400',
  LabelRemove: 'bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-400',
}

const DEFAULT_EVENT_COLOR = 'bg-gray-100 text-gray-700 dark:bg-gray-800 dark:text-gray-400'

/** 执行结果状态 */
function getResultStatus(event: HookEventRecord): 'success' | 'fail' | 'none' {
  if (event.type !== 'HookResult') return 'none'
  if (event.data.blocked) return 'fail'
  return event.data.success ? 'success' : 'fail'
}

const RESULT_COLORS = {
  success: 'bg-emerald-100 text-emerald-700 dark:bg-emerald-900/30 dark:text-emerald-400',
  fail: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
  none: 'bg-gray-100 text-gray-600 dark:bg-gray-800 dark:text-gray-400',
}

const RESULT_LABELS = { success: '成功', fail: '失败', none: '事件' }

/** 格式化为相对时间 */
function formatRelativeTime(isoTime: string): string {
  const diff = Date.now() - new Date(isoTime).getTime()
  if (diff < 0) return '刚刚'
  const seconds = Math.floor(diff / 1000)
  if (seconds < 60) return '刚刚'
  const minutes = Math.floor(seconds / 60)
  if (minutes < 60) return `${minutes}分钟前`
  const hours = Math.floor(minutes / 60)
  if (hours < 24) return `${hours}小时前`
  const days = Math.floor(hours / 24)
  if (days < 30) return `${days}天前`
  return new Date(isoTime).toLocaleDateString()
}

/** Hook type 标签 */
const TYPE_COLORS: Record<string, string> = {
  command: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  prompt: 'bg-cyan-100 text-cyan-700 dark:bg-cyan-900/30 dark:text-cyan-400',
}

/** 获取 hook 内容摘要 */
function getHookSummary(hooks: HookDefinition[]): string {
  if (!hooks.length) return '-'
  const first = hooks[0]
  const text = first.type === 'command' ? first.command : first.prompt
  const summary = text && text.length > 60 ? text.slice(0, 60) + '...' : text || '-'
  return hooks.length > 1 ? `${summary} (+${hooks.length - 1})` : summary
}

export default function ScheduledTasksView() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''

  const [config, setConfig] = useState<HooksConfig>({ hooks: {} })
  const [loading, setLoading] = useState(true)

  // Tab 切换 + 执行记录状态
  const [activeTab, setActiveTab] = useState<'config' | 'logs'>('config')
  const [events, setEvents] = useState<HookEventRecord[]>([])
  const [eventsLoading, setEventsLoading] = useState(false)
  const [eventTypeFilter, setEventTypeFilter] = useState('')

  /** 加载 hooks.json */
  const loadHooks = useCallback(async () => {
    if (!workspaceId) return
    setLoading(true)
    try {
      const data = await window.electronAPI.creatorMedia.hooks.read(workspaceId)
      setConfig(data as HooksConfig)
    } catch {
      setConfig({ hooks: {} })
    } finally {
      setLoading(false)
    }
  }, [workspaceId])

  useEffect(() => { loadHooks() }, [loadHooks])

  /** 加载执行记录 */
  const loadEvents = useCallback(async () => {
    if (!workspaceId) return
    setEventsLoading(true)
    try {
      const data = await window.electronAPI.creatorMedia.hookEvents.list(
        workspaceId,
        { limit: 50, eventType: eventTypeFilter || undefined }
      )
      setEvents(data as HookEventRecord[])
    } catch {
      setEvents([])
    } finally {
      setEventsLoading(false)
    }
  }, [workspaceId, eventTypeFilter])

  // 切换到执行记录 Tab 时自动加载
  useEffect(() => {
    if (activeTab === 'logs') loadEvents()
  }, [activeTab, loadEvents])

  /** 保存 hooks.json */
  const saveHooks = useCallback(async (newConfig: HooksConfig) => {
    if (!workspaceId) return
    await window.electronAPI.creatorMedia.hooks.write(workspaceId, newConfig)
    setConfig(newConfig)
  }, [workspaceId])

  /** 编辑中的卡片 */
  const [editingCard, setEditingCard] = useState<HookCard | null>(null)
  const [editJson, setEditJson] = useState('')
  const [editError, setEditError] = useState('')
  const textareaRef = useRef<HTMLTextAreaElement>(null)

  /** 打开编辑弹窗 */
  const handleEdit = useCallback((card: HookCard) => {
    setEditingCard(card)
    setEditJson(JSON.stringify(card.matcher, null, 2))
    setEditError('')
  }, [])

  /** 保存编辑 */
  const handleEditSave = useCallback(async () => {
    if (!editingCard) return
    let parsed: HookMatcher
    try {
      parsed = JSON.parse(editJson)
    } catch {
      setEditError('JSON 格式错误')
      return
    }
    if (!parsed.hooks || !Array.isArray(parsed.hooks)) {
      setEditError('缺少 hooks 数组')
      return
    }
    const newConfig = structuredClone(config)
    const matchers = newConfig.hooks[editingCard.eventName]
    if (!matchers?.[editingCard.matcherIndex]) return
    matchers[editingCard.matcherIndex] = parsed
    await saveHooks(newConfig)
    setEditingCard(null)
  }, [editingCard, editJson, config, saveHooks])

  /** 切换 enabled 状态 */
  const handleToggle = useCallback(async (eventName: string, matcherIndex: number) => {
    const newConfig = structuredClone(config)
    const matchers = newConfig.hooks[eventName]
    if (!matchers?.[matcherIndex]) return
    const current = matchers[matcherIndex].enabled
    matchers[matcherIndex].enabled = current === false ? true : false
    await saveHooks(newConfig)
  }, [config, saveHooks])

  /** 删除 matcher */
  const handleDelete = useCallback(async (eventName: string, matcherIndex: number) => {
    const confirmed = window.confirm(t('确定要删除这个 Hook 吗？此操作不可撤销。'))
    if (!confirmed) return
    const newConfig = structuredClone(config)
    const matchers = newConfig.hooks[eventName]
    if (!matchers) return
    matchers.splice(matcherIndex, 1)
    // 如果该事件下没有 matcher 了，移除整个事件 key
    if (matchers.length === 0) {
      delete newConfig.hooks[eventName]
    }
    await saveHooks(newConfig)
  }, [config, saveHooks, t])

  /** 扁平化为卡片列表 */
  const cards: HookCard[] = []
  for (const [eventName, matchers] of Object.entries(config.hooks)) {
    if (!Array.isArray(matchers)) continue
    matchers.forEach((matcher, index) => {
      cards.push({ eventName, matcherIndex: index, matcher })
    })
  }

  // EditPopover 配置
  const editConfig = activeWorkspace
    ? getEditConfig('creator-media-edit-scheduled-task', activeWorkspace.rootPath)
    : null

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
          <h1 className="text-base font-semibold text-foreground">{t('定时任务')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {t('管理工作区 Hooks 配置')}
          </p>
        </div>
        {/* Tab 切换（居中） */}
        <div className="titlebar-no-drag absolute left-1/2 -translate-x-1/2 flex items-center gap-0.5 rounded-lg bg-muted/40 p-0.5">
          <button
            type="button"
            onClick={() => setActiveTab('config')}
            className={`rounded-md px-4 py-1 text-sm font-medium transition-colors ${
              activeTab === 'config'
                ? 'bg-background text-foreground shadow-sm'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            {t('配置')}
          </button>
          <button
            type="button"
            onClick={() => setActiveTab('logs')}
            className={`rounded-md px-4 py-1 text-sm font-medium transition-colors ${
              activeTab === 'logs'
                ? 'bg-background text-foreground shadow-sm'
                : 'text-muted-foreground hover:text-foreground'
            }`}
          >
            {t('执行记录')}
          </button>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          {activeTab === 'config' && editConfig && (
            <EditPopover
              trigger={<EditButton />}
              context={editConfig.context}
              example={editConfig.example}
              model={editConfig.model}
              systemPromptPreset={editConfig.systemPromptPreset}
              inlineExecution={editConfig.inlineExecution}
              overridePlaceholder={editConfig.overridePlaceholder}
            />
          )}
          {/* 刷新按钮 */}
          <button
            type="button"
            onClick={activeTab === 'config' ? loadHooks : loadEvents}
            className="rounded-md px-2.5 py-1.5 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors"
            title={t('刷新')}
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0 3.181 3.183a8.25 8.25 0 0 0 13.803-3.7M4.031 9.865a8.25 8.25 0 0 1 13.803-3.7l3.181 3.182" />
            </svg>
          </button>
        </div>
      </div>

      {/* 配置 Tab */}
      {activeTab === 'config' && (
        <div className="flex-1 overflow-auto px-6 py-4">
          {cards.length > 0 ? (
            <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
              {cards.map((card) => {
                const { eventName, matcherIndex, matcher } = card
                const isEnabled = matcher.enabled !== false
                const hookType = matcher.hooks?.[0]?.type || 'command'

                return (
                  <div
                    key={`${eventName}-${matcherIndex}`}
                    className={`rounded-lg border bg-background/60 p-4 transition-colors cursor-pointer hover:border-accent/50 ${
                      isEnabled ? 'border-border/60' : 'border-border/30 opacity-60'
                    }`}
                    onClick={() => handleEdit(card)}
                  >
                    {/* 卡片头部：事件标签 + type 标签 + 开关 */}
                    <div className="flex items-center justify-between mb-3">
                      <div className="flex items-center gap-1.5">
                        <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${
                          EVENT_COLORS[eventName] || DEFAULT_EVENT_COLOR
                        }`}>
                          {t(EVENT_LABELS[eventName] || eventName)}
                        </span>
                        <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${
                          TYPE_COLORS[hookType] || 'bg-gray-100 text-gray-600 dark:bg-gray-800 dark:text-gray-400'
                        }`}>
                          {hookType}
                        </span>
                      </div>
                      <div className="flex items-center gap-1.5" onClick={(e) => e.stopPropagation()}>
                        {/* 启用/禁用开关 */}
                        <button
                          type="button"
                          onClick={() => handleToggle(eventName, matcherIndex)}
                          className={`relative inline-flex h-4 w-7 items-center rounded-full transition-colors ${
                            isEnabled ? 'bg-accent' : 'bg-muted'
                          }`}
                          title={isEnabled ? t('禁用') : t('启用')}
                        >
                          <span className={`inline-block h-2.5 w-2.5 rounded-full bg-white transition-transform ${
                            isEnabled ? 'translate-x-3.5' : 'translate-x-0.5'
                          }`} />
                        </button>
                        {/* 删除 */}
                        <button
                          type="button"
                          onClick={() => handleDelete(eventName, matcherIndex)}
                          className="rounded px-1 py-0.5 text-red-500/70 hover:text-red-500 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors"
                          title={t('删除')}
                        >
                          <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                            <path strokeLinecap="round" strokeLinejoin="round" d="m14.74 9-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 0 1-2.244 2.077H8.084a2.25 2.25 0 0 1-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 0 0-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 0 1 3.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 0 0-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 0 0-7.5 0" />
                          </svg>
                        </button>
                      </div>
                    </div>

                    {/* 调度信息 */}
                    <div className="space-y-1.5 text-xs">
                      {matcher.cron && (
                        <div className="flex items-center gap-2">
                          <span className="text-muted-foreground w-12 shrink-0">cron</span>
                          <code className="font-mono text-foreground/80 bg-muted/50 px-1.5 py-0.5 rounded text-[11px]">
                            {matcher.cron}
                          </code>
                          {matcher.timezone && (
                            <span className="text-muted-foreground/60 text-[10px]">{matcher.timezone}</span>
                          )}
                        </div>
                      )}
                      {matcher.matcher && (
                        <div className="flex items-center gap-2">
                          <span className="text-muted-foreground w-12 shrink-0">match</span>
                          <code className="font-mono text-foreground/80 bg-muted/50 px-1.5 py-0.5 rounded text-[11px]">
                            {matcher.matcher}
                          </code>
                        </div>
                      )}
                      <div className="flex items-start gap-2">
                        <span className="text-muted-foreground w-12 shrink-0">{hookType}</span>
                        <span className="text-foreground/70 break-all leading-relaxed">
                          {getHookSummary(matcher.hooks || [])}
                        </span>
                      </div>
                      {matcher.labels && matcher.labels.length > 0 && (
                        <div className="flex items-center gap-2">
                          <span className="text-muted-foreground w-12 shrink-0">labels</span>
                          <div className="flex flex-wrap gap-1">
                            {matcher.labels.map((label) => (
                              <span key={label} className="text-[10px] bg-muted/60 px-1.5 py-0.5 rounded">
                                {label}
                              </span>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>
                  </div>
                )
              })}
            </div>
          ) : (
            /* 空状态 */
            <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-12">
              <div className="text-center space-y-3">
                <svg className="mx-auto w-10 h-10 text-muted-foreground/40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M12 6v6h4.5m4.5 0a9 9 0 1 1-18 0 9 9 0 0 1 18 0Z" />
                </svg>
                <div>
                  <p className="text-sm font-medium text-foreground">{t('还没有配置 Hook')}</p>
                  <p className="mt-1 text-xs text-muted-foreground">
                    {t('点击右上角编辑按钮，通过 AI 助手快速创建')}
                  </p>
                </div>
              </div>
            </div>
          )}
        </div>
      )}

      {/* 执行记录 Tab */}
      {activeTab === 'logs' && (
        <div className="flex-1 overflow-auto px-6 py-4">
          {/* 筛选栏 */}
          <div className="flex items-center gap-2 mb-4">
            <select
              value={eventTypeFilter}
              onChange={(e) => setEventTypeFilter(e.target.value)}
              className="rounded-md border border-border/60 bg-background px-2.5 py-1.5 text-xs text-foreground focus:outline-none focus:ring-1 focus:ring-accent/50"
            >
              <option value="">{t('全部事件')}</option>
              <option value="HookResult">{t('执行结果')}</option>
              <option value="SchedulerTick">{t('定时调度')}</option>
              <option value="UserPromptSubmit">{t('用户提交')}</option>
              <option value="SessionStart">{t('会话开始')}</option>
              <option value="SessionEnd">{t('会话结束')}</option>
              <option value="LabelAdd">{t('标签添加')}</option>
              <option value="LabelRemove">{t('标签移除')}</option>
            </select>
          </div>

          {eventsLoading ? (
            <div className="flex items-center justify-center py-12">
              <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
            </div>
          ) : events.length > 0 ? (
            <div className="space-y-2">
              {events.map((evt) => {
                const resultStatus = getResultStatus(evt)
                const eventLabel = EVENT_LABELS[evt.data.event as string] || EVENT_LABELS[evt.type] || evt.type
                const eventColor = EVENT_COLORS[evt.data.event as string] || EVENT_COLORS[evt.type] || DEFAULT_EVENT_COLOR
                const commandOrPrompt = (evt.data.command as string) || (evt.data.prompt as string) || ''

                return (
                  <div
                    key={evt.id}
                    className="rounded-lg border border-border/60 bg-background/60 px-4 py-3"
                  >
                    <div className="flex items-center justify-between">
                      <div className="flex items-center gap-1.5">
                        <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${eventColor}`}>
                          {t(eventLabel)}
                        </span>
                        {evt.type === 'HookResult' && (
                          <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${RESULT_COLORS[resultStatus]}`}>
                            {t(RESULT_LABELS[resultStatus])}
                          </span>
                        )}
                        {evt.data.hookType && (
                          <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${
                            TYPE_COLORS[evt.data.hookType as string] || 'bg-gray-100 text-gray-600 dark:bg-gray-800 dark:text-gray-400'
                          }`}>
                            {evt.data.hookType as string}
                          </span>
                        )}
                      </div>
                      <div className="flex items-center gap-2 text-[11px] text-muted-foreground">
                        {evt.durationMs > 0 && <span>{evt.durationMs}ms</span>}
                        <span>{formatRelativeTime(evt.time)}</span>
                      </div>
                    </div>
                    {commandOrPrompt && (
                      <div className="mt-2">
                        <code className="text-[11px] font-mono text-foreground/70 bg-muted/40 px-2 py-1 rounded block break-all">
                          {evt.data.hookType === 'command' ? '$ ' : ''}{commandOrPrompt.length > 120 ? commandOrPrompt.slice(0, 120) + '...' : commandOrPrompt}
                        </code>
                      </div>
                    )}
                  </div>
                )
              })}
            </div>
          ) : (
            <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-12">
              <div className="text-center space-y-3">
                <svg className="mx-auto w-10 h-10 text-muted-foreground/40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M19.5 14.25v-2.625a3.375 3.375 0 0 0-3.375-3.375h-1.5A1.125 1.125 0 0 1 13.5 7.125v-1.5a3.375 3.375 0 0 0-3.375-3.375H8.25m0 12.75h7.5m-7.5 3H12M10.5 2.25H5.625c-.621 0-1.125.504-1.125 1.125v17.25c0 .621.504 1.125 1.125 1.125h12.75c.621 0 1.125-.504 1.125-1.125V11.25a9 9 0 0 0-9-9Z" />
                </svg>
                <div>
                  <p className="text-sm font-medium text-foreground">{t('暂无执行记录')}</p>
                  <p className="mt-1 text-xs text-muted-foreground">
                    {t('Hook 触发后，执行记录将显示在这里')}
                  </p>
                </div>
              </div>
            </div>
          )}
        </div>
      )}

      {/* 编辑弹窗 */}
      {editingCard && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40" onClick={() => setEditingCard(null)}>
          <div
            className="bg-background border border-border rounded-xl shadow-xl w-full max-w-lg mx-4 flex flex-col max-h-[80vh]"
            onClick={(e) => e.stopPropagation()}
          >
            {/* 弹窗头部 */}
            <div className="flex items-center justify-between px-5 py-3.5 border-b border-border/40">
              <div className="flex items-center gap-2">
                <h2 className="text-sm font-semibold text-foreground">{t('编辑 Hook')}</h2>
                <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${
                  EVENT_COLORS[editingCard.eventName] || DEFAULT_EVENT_COLOR
                }`}>
                  {t(EVENT_LABELS[editingCard.eventName] || editingCard.eventName)}
                </span>
              </div>
              <button
                type="button"
                onClick={() => setEditingCard(null)}
                className="text-muted-foreground hover:text-foreground transition-colors"
              >
                <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M6 18 18 6M6 6l12 12" />
                </svg>
              </button>
            </div>
            {/* JSON 编辑区 */}
            <div className="flex-1 overflow-auto px-5 py-4">
              <textarea
                ref={textareaRef}
                value={editJson}
                onChange={(e) => { setEditJson(e.target.value); setEditError('') }}
                className="w-full h-64 font-mono text-xs bg-muted/30 border border-border/60 rounded-lg p-3 resize-y focus:outline-none focus:ring-1 focus:ring-accent/50 text-foreground"
                spellCheck={false}
              />
              {editError && (
                <p className="mt-2 text-xs text-red-500">{editError}</p>
              )}
            </div>
            {/* 弹窗底部 */}
            <div className="flex items-center justify-end gap-2 px-5 py-3 border-t border-border/40">
              <button
                type="button"
                onClick={() => setEditingCard(null)}
                className="rounded-md px-3 py-1.5 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors"
              >
                {t('取消')}
              </button>
              <button
                type="button"
                onClick={handleEditSave}
                className="rounded-md px-3 py-1.5 text-xs font-medium bg-accent text-accent-foreground hover:bg-accent/90 transition-colors"
              >
                {t('保存')}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
