/**
 * 定时任务管理视图
 *
 * 支持查看/管理 scheduled_tasks（通用定时任务）和 review_tasks（采集调度任务，只读）
 */

import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { ScheduledTaskDialog } from './components/ScheduledTaskDialog'
import { EditPopover, EditButton, getEditConfig } from '@/components/ui/EditPopover'
import type { ScheduledTask, ScheduledTaskType, ScheduledTaskStatus, ReviewTask } from '@sprouty-ai/shared/db/types'

/** Tab 定义 */
type TabId = 'all' | 'review' | 'publish' | 'custom'
const TABS: { id: TabId; label: string }[] = [
  { id: 'all', label: '全部' },
  { id: 'review', label: '采集任务' },
  { id: 'publish', label: '发布任务' },
  { id: 'custom', label: '自定义' },
]

/** 任务类型 badge 配置 */
const TASK_TYPE_CONFIG: Record<string, { label: string; className: string }> = {
  review: { label: '采集', className: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400' },
  publish: { label: '发布', className: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400' },
  collect: { label: '采集', className: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400' },
  custom: { label: '自定义', className: 'bg-gray-100 text-gray-700 dark:bg-gray-800 dark:text-gray-400' },
}

/** 状态指示器配置 */
const STATUS_CONFIG: Record<string, { label: string; dotClass: string }> = {
  active: { label: '运行中', dotClass: 'bg-green-500' },
  paused: { label: '已暂停', dotClass: 'bg-yellow-500' },
  error: { label: '错误', dotClass: 'bg-red-500' },
  completed: { label: '已完成', dotClass: 'bg-gray-400' },
  // review_tasks 状态
  pending: { label: '待执行', dotClass: 'bg-yellow-500' },
  executing: { label: '执行中', dotClass: 'bg-blue-500' },
  failed: { label: '失败', dotClass: 'bg-red-500' },
  cancelled: { label: '已取消', dotClass: 'bg-gray-400' },
}

/** 调度描述可读化 */
function describeSchedule(task: ScheduledTask): string {
  if (task.schedule_mode === 'cron' && task.cron_expression) {
    return `cron: ${task.cron_expression}`
  }
  if (task.schedule_mode === 'interval' && task.interval_seconds) {
    const secs = task.interval_seconds
    if (secs >= 3600) return `每 ${secs / 3600} 小时`
    if (secs >= 60) return `每 ${secs / 60} 分钟`
    return `每 ${secs} 秒`
  }
  if (task.schedule_mode === 'once' && task.scheduled_at) {
    return `单次: ${formatTime(task.scheduled_at)}`
  }
  return '-'
}

/** 格式化时间 */
function formatTime(dateStr: string | null): string {
  if (!dateStr) return '-'
  try {
    const d = new Date(dateStr)
    return d.toLocaleString('zh-CN', { month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit' })
  } catch {
    return dateStr
  }
}

/** 统一行数据类型 */
interface TaskRow {
  id: string
  name: string
  taskType: string
  schedule: string
  status: string
  enabled: boolean
  lastRunAt: string | null
  nextRunAt: string | null
  isReviewTask: boolean
  raw: ScheduledTask | ReviewTask
}

export default function ScheduledTasksView() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const { projects, activeProject, loading, switchProject } = useCreatorMedia()

  const [activeTab, setActiveTab] = useState<TabId>('all')
  const [scheduledTasks, setScheduledTasks] = useState<ScheduledTask[]>([])
  const [reviewTasks, setReviewTasks] = useState<ReviewTask[]>([])
  const [dataLoading, setDataLoading] = useState(false)
  const [showDialog, setShowDialog] = useState(false)
  const [editingTask, setEditingTask] = useState<ScheduledTask | null>(null)

  /** 加载数据 */
  const loadData = useCallback(async () => {
    if (!workspaceId) return
    setDataLoading(true)
    try {
      const [tasks, reviews] = await Promise.all([
        window.electronAPI.creatorMedia.scheduledTasks.list(workspaceId),
        window.electronAPI.creatorMedia.reviewTasksAll.list(workspaceId),
      ])
      setScheduledTasks(tasks as ScheduledTask[])
      setReviewTasks(reviews as ReviewTask[])
    } catch {
      setScheduledTasks([])
      setReviewTasks([])
    } finally {
      setDataLoading(false)
    }
  }, [workspaceId])

  useEffect(() => { loadData() }, [loadData])

  /** 合并为统一行数据 */
  const allRows: TaskRow[] = [
    ...scheduledTasks.map((task): TaskRow => ({
      id: task.id,
      name: task.name,
      taskType: task.task_type,
      schedule: describeSchedule(task),
      status: task.status,
      enabled: task.enabled === 1,
      lastRunAt: task.last_run_at,
      nextRunAt: task.next_run_at,
      isReviewTask: false,
      raw: task,
    })),
    ...reviewTasks.map((rt): TaskRow => ({
      id: rt.id,
      name: `采集 #${rt.publish_record_id.slice(0, 8)} (${rt.review_type})`,
      taskType: 'review',
      schedule: formatTime(rt.scheduled_at),
      status: rt.status,
      enabled: rt.status === 'pending',
      lastRunAt: rt.executed_at,
      nextRunAt: rt.status === 'pending' ? rt.scheduled_at : null,
      isReviewTask: true,
      raw: rt,
    })),
  ]

  /** Tab 过滤 */
  const filteredRows = allRows.filter((row) => {
    if (activeTab === 'all') return true
    if (activeTab === 'review') return row.taskType === 'review'
    if (activeTab === 'publish') return row.taskType === 'publish'
    if (activeTab === 'custom') return row.taskType === 'custom' || row.taskType === 'collect'
    return true
  })

  /** 创建/更新任务 */
  const handleSave = useCallback(async (data: Partial<ScheduledTask> & { id: string; name: string }) => {
    if (!workspaceId) return
    if (editingTask) {
      await window.electronAPI.creatorMedia.scheduledTasks.update(workspaceId, data.id, data)
    } else {
      await window.electronAPI.creatorMedia.scheduledTasks.create(workspaceId, {
        ...data,
        project_id: activeProject?.id || null,
      })
    }
    await loadData()
  }, [workspaceId, editingTask, activeProject, loadData])

  /** 切换启用/禁用 */
  const handleToggle = useCallback(async (task: ScheduledTask) => {
    if (!workspaceId) return
    await window.electronAPI.creatorMedia.scheduledTasks.toggle(workspaceId, task.id, task.enabled !== 1)
    await loadData()
  }, [workspaceId, loadData])

  /** 删除任务 */
  const handleDelete = useCallback(async (task: ScheduledTask) => {
    if (!workspaceId) return
    const confirmed = window.confirm(t(`确定要删除任务 "${task.name}" 吗？此操作不可撤销。`))
    if (!confirmed) return
    await window.electronAPI.creatorMedia.scheduledTasks.delete(workspaceId, task.id)
    await loadData()
  }, [workspaceId, loadData, t])

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
            {t('管理定时发布、采集和自定义任务')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          {/* AI 辅助编辑 */}
          {editConfig && (
            <EditPopover
              trigger={<EditButton />}
              context={editConfig.context}
              example={editConfig.example}
              model={editConfig.model}
              systemPromptPreset={editConfig.systemPromptPreset}
              inlineExecution={editConfig.inlineExecution}
            />
          )}
          <button
            type="button"
            onClick={() => {
              setEditingTask(null)
              setShowDialog(true)
            }}
            className="inline-flex items-center gap-1.5 rounded-md bg-foreground px-3 py-1.5 text-xs font-medium text-background hover:bg-foreground/90 transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
            </svg>
            {t('新建任务')}
          </button>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* Tab 栏 */}
      <div className="px-6 pt-3">
        <div className="flex gap-1 border-b border-border/30">
          {TABS.map((tab) => (
            <button
              key={tab.id}
              type="button"
              onClick={() => setActiveTab(tab.id)}
              className={`px-3 py-1.5 text-xs font-medium transition-colors border-b-2 -mb-px ${
                activeTab === tab.id
                  ? 'border-foreground text-foreground'
                  : 'border-transparent text-muted-foreground hover:text-foreground'
              }`}
            >
              {t(tab.label)}
            </button>
          ))}
        </div>
      </div>

      {/* 数据表 */}
      <div className="flex-1 overflow-auto px-6 py-4">
        {dataLoading ? (
          <div className="flex items-center justify-center py-12">
            <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
          </div>
        ) : filteredRows.length > 0 ? (
          <table className="w-full text-xs">
            <thead>
              <tr className="border-b border-border/40 text-muted-foreground">
                <th className="text-left py-2 pr-3 font-medium">{t('名称')}</th>
                <th className="text-left py-2 pr-3 font-medium">{t('类型')}</th>
                <th className="text-left py-2 pr-3 font-medium">{t('调度')}</th>
                <th className="text-left py-2 pr-3 font-medium">{t('状态')}</th>
                <th className="text-left py-2 pr-3 font-medium">{t('上次运行')}</th>
                <th className="text-left py-2 pr-3 font-medium">{t('下次运行')}</th>
                <th className="text-right py-2 font-medium">{t('操作')}</th>
              </tr>
            </thead>
            <tbody>
              {filteredRows.map((row) => (
                <tr key={row.id} className="border-b border-border/20 hover:bg-muted/30 transition-colors">
                  {/* 名称 */}
                  <td className="py-2.5 pr-3">
                    <span className="text-sm text-foreground">{row.name}</span>
                  </td>
                  {/* 类型 badge */}
                  <td className="py-2.5 pr-3">
                    <span className={`inline-flex rounded-full px-2 py-0.5 text-[10px] font-medium ${
                      TASK_TYPE_CONFIG[row.taskType]?.className || TASK_TYPE_CONFIG.custom.className
                    }`}>
                      {t(TASK_TYPE_CONFIG[row.taskType]?.label || row.taskType)}
                    </span>
                  </td>
                  {/* 调度 */}
                  <td className="py-2.5 pr-3 text-muted-foreground font-mono text-[11px]">
                    {row.schedule}
                  </td>
                  {/* 状态 */}
                  <td className="py-2.5 pr-3">
                    <div className="flex items-center gap-1.5">
                      <span className={`w-1.5 h-1.5 rounded-full ${
                        STATUS_CONFIG[row.status]?.dotClass || 'bg-gray-400'
                      }`} />
                      <span className="text-muted-foreground">
                        {t(STATUS_CONFIG[row.status]?.label || row.status)}
                      </span>
                    </div>
                  </td>
                  {/* 上次运行 */}
                  <td className="py-2.5 pr-3 text-muted-foreground">
                    {formatTime(row.lastRunAt)}
                  </td>
                  {/* 下次运行 */}
                  <td className="py-2.5 pr-3 text-muted-foreground">
                    {formatTime(row.nextRunAt)}
                  </td>
                  {/* 操作 */}
                  <td className="py-2.5 text-right">
                    {row.isReviewTask ? (
                      <span className="text-[10px] text-muted-foreground/60">{t('只读')}</span>
                    ) : (
                      <div className="flex items-center justify-end gap-1">
                        {/* 启用/禁用 */}
                        <button
                          type="button"
                          onClick={() => handleToggle(row.raw as ScheduledTask)}
                          className={`relative inline-flex h-4 w-7 items-center rounded-full transition-colors ${
                            row.enabled ? 'bg-accent' : 'bg-muted'
                          }`}
                          title={row.enabled ? t('禁用') : t('启用')}
                        >
                          <span className={`inline-block h-2.5 w-2.5 rounded-full bg-white transition-transform ${
                            row.enabled ? 'translate-x-3.5' : 'translate-x-0.5'
                          }`} />
                        </button>
                        {/* 编辑 */}
                        <button
                          type="button"
                          onClick={() => {
                            setEditingTask(row.raw as ScheduledTask)
                            setShowDialog(true)
                          }}
                          className="rounded px-1.5 py-0.5 text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors"
                          title={t('编辑')}
                        >
                          <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                            <path strokeLinecap="round" strokeLinejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L10.582 16.07a4.5 4.5 0 0 1-1.897 1.13L6 18l.8-2.685a4.5 4.5 0 0 1 1.13-1.897l8.932-8.931Zm0 0L19.5 7.125" />
                          </svg>
                        </button>
                        {/* 删除 */}
                        <button
                          type="button"
                          onClick={() => handleDelete(row.raw as ScheduledTask)}
                          className="rounded px-1.5 py-0.5 text-red-500/70 hover:text-red-500 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors"
                          title={t('删除')}
                        >
                          <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                            <path strokeLinecap="round" strokeLinejoin="round" d="m14.74 9-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 0 1-2.244 2.077H8.084a2.25 2.25 0 0 1-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 0 0-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 0 1 3.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 0 0-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 0 0-7.5 0" />
                          </svg>
                        </button>
                      </div>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        ) : (
          /* 空状态 */
          <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-12">
            <div className="text-center space-y-3">
              <svg className="mx-auto w-10 h-10 text-muted-foreground/40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M12 6v6h4.5m4.5 0a9 9 0 1 1-18 0 9 9 0 0 1 18 0Z" />
              </svg>
              <div>
                <p className="text-sm font-medium text-foreground">{t('还没有定时任务')}</p>
                <p className="mt-1 text-xs text-muted-foreground">
                  {t('创建定时任务来自动执行发布、采集等操作')}
                </p>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* 新建/编辑 Dialog */}
      <ScheduledTaskDialog
        open={showDialog}
        onClose={() => {
          setShowDialog(false)
          setEditingTask(null)
        }}
        onSave={handleSave}
        editingTask={editingTask}
      />
    </div>
  )
}
