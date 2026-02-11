/**
 * 定时任务新建/编辑 Dialog
 */

import { useState, useEffect } from 'react'
import { useT } from '@/context/LocaleContext'
import type { ScheduledTask, ScheduledTaskType, ScheduleMode, ScheduledTaskStatus } from '@sprouty-ai/shared/db/types'

/** 任务类型选项 */
const TASK_TYPE_OPTIONS: { id: ScheduledTaskType; label: string }[] = [
  { id: 'publish', label: '定时发布' },
  { id: 'collect', label: '定时采集' },
  { id: 'custom', label: '自定义' },
]

/** 调度模式选项 */
const SCHEDULE_MODE_OPTIONS: { id: ScheduleMode; label: string; desc: string }[] = [
  { id: 'cron', label: 'Cron 表达式', desc: '灵活的定时规则' },
  { id: 'interval', label: '固定间隔', desc: '每隔一段时间执行' },
  { id: 'once', label: '单次执行', desc: '指定时间执行一次' },
]

/** 间隔单位 */
const INTERVAL_UNITS = [
  { id: 'seconds', label: '秒', multiplier: 1 },
  { id: 'minutes', label: '分钟', multiplier: 60 },
  { id: 'hours', label: '小时', multiplier: 3600 },
]

interface ScheduledTaskDialogProps {
  open: boolean
  onClose: () => void
  onSave: (data: Partial<ScheduledTask> & { id: string; name: string }) => Promise<void>
  editingTask?: ScheduledTask | null
}

export function ScheduledTaskDialog({ open, onClose, onSave, editingTask }: ScheduledTaskDialogProps) {
  const t = useT()
  const isEditing = !!editingTask

  // 表单状态
  const [name, setName] = useState('')
  const [description, setDescription] = useState('')
  const [taskType, setTaskType] = useState<ScheduledTaskType>('custom')
  const [scheduleMode, setScheduleMode] = useState<ScheduleMode>('cron')
  const [cronExpression, setCronExpression] = useState('0 */6 * * *')
  const [intervalValue, setIntervalValue] = useState(60)
  const [intervalUnit, setIntervalUnit] = useState('minutes')
  const [scheduledAt, setScheduledAt] = useState('')
  const [enabled, setEnabled] = useState(true)
  const [payload, setPayload] = useState('')
  const [saving, setSaving] = useState(false)

  // 编辑时填充表单
  useEffect(() => {
    if (editingTask) {
      setName(editingTask.name)
      setDescription(editingTask.description || '')
      setTaskType(editingTask.task_type)
      setScheduleMode(editingTask.schedule_mode)
      setCronExpression(editingTask.cron_expression || '0 */6 * * *')
      setEnabled(editingTask.enabled === 1)
      setPayload(editingTask.payload || '')
      setScheduledAt(editingTask.scheduled_at || '')
      // 解析 interval_seconds 到合适的单位
      if (editingTask.interval_seconds) {
        const secs = editingTask.interval_seconds
        if (secs >= 3600 && secs % 3600 === 0) {
          setIntervalValue(secs / 3600)
          setIntervalUnit('hours')
        } else if (secs >= 60 && secs % 60 === 0) {
          setIntervalValue(secs / 60)
          setIntervalUnit('minutes')
        } else {
          setIntervalValue(secs)
          setIntervalUnit('seconds')
        }
      }
    } else {
      // 重置表单
      setName('')
      setDescription('')
      setTaskType('custom')
      setScheduleMode('cron')
      setCronExpression('0 */6 * * *')
      setIntervalValue(60)
      setIntervalUnit('minutes')
      setScheduledAt('')
      setEnabled(true)
      setPayload('')
    }
  }, [editingTask, open])

  const handleSubmit = async () => {
    if (!name.trim()) return
    setSaving(true)
    try {
      const unit = INTERVAL_UNITS.find(u => u.id === intervalUnit)
      const intervalSeconds = scheduleMode === 'interval' ? intervalValue * (unit?.multiplier ?? 1) : null

      await onSave({
        id: editingTask?.id || crypto.randomUUID(),
        name: name.trim(),
        description: description.trim() || null,
        task_type: taskType,
        schedule_mode: scheduleMode,
        cron_expression: scheduleMode === 'cron' ? cronExpression : null,
        interval_seconds: intervalSeconds,
        scheduled_at: scheduleMode === 'once' ? scheduledAt || null : null,
        enabled: enabled ? 1 : 0,
        status: editingTask?.status || ('active' as ScheduledTaskStatus),
        payload: payload.trim() || null,
      })
      onClose()
    } finally {
      setSaving(false)
    }
  }

  if (!open) return null

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center">
      {/* 遮罩 */}
      <div className="absolute inset-0 bg-black/40" onClick={onClose} />

      {/* 对话框 */}
      <div className="relative z-10 w-full max-w-lg rounded-xl border border-border/60 bg-background shadow-xl mx-4">
        {/* 标题 */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-border/40">
          <h2 className="text-sm font-semibold text-foreground">
            {isEditing ? t('编辑定时任务') : t('新建定时任务')}
          </h2>
          <button type="button" onClick={onClose} className="text-muted-foreground hover:text-foreground">
            <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* 表单 */}
        <div className="px-6 py-4 space-y-4 max-h-[60vh] overflow-y-auto">
          {/* 名称 */}
          <div>
            <label className="block text-xs font-medium text-foreground mb-1">{t('任务名称')}</label>
            <input
              type="text"
              value={name}
              onChange={(e) => setName(e.target.value)}
              placeholder={t('输入任务名称')}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-accent"
            />
          </div>

          {/* 描述 */}
          <div>
            <label className="block text-xs font-medium text-foreground mb-1">{t('描述')}</label>
            <textarea
              value={description}
              onChange={(e) => setDescription(e.target.value)}
              placeholder={t('可选：描述任务用途')}
              rows={2}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-accent resize-none"
            />
          </div>

          {/* 任务类型 */}
          <div>
            <label className="block text-xs font-medium text-foreground mb-1">{t('任务类型')}</label>
            <select
              value={taskType}
              onChange={(e) => setTaskType(e.target.value as ScheduledTaskType)}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-accent"
            >
              {TASK_TYPE_OPTIONS.map((opt) => (
                <option key={opt.id} value={opt.id}>{t(opt.label)}</option>
              ))}
            </select>
          </div>

          {/* 调度模式 */}
          <div>
            <label className="block text-xs font-medium text-foreground mb-2">{t('调度模式')}</label>
            <div className="flex gap-2">
              {SCHEDULE_MODE_OPTIONS.map((opt) => (
                <button
                  key={opt.id}
                  type="button"
                  onClick={() => setScheduleMode(opt.id)}
                  className={`flex-1 rounded-md border px-3 py-2 text-left transition-colors ${
                    scheduleMode === opt.id
                      ? 'border-accent bg-accent/10 text-foreground'
                      : 'border-border/60 bg-background text-muted-foreground hover:border-border'
                  }`}
                >
                  <div className="text-xs font-medium">{t(opt.label)}</div>
                  <div className="text-[10px] mt-0.5 opacity-70">{t(opt.desc)}</div>
                </button>
              ))}
            </div>
          </div>

          {/* 调度配置 — 根据模式显示不同输入 */}
          {scheduleMode === 'cron' && (
            <div>
              <label className="block text-xs font-medium text-foreground mb-1">{t('Cron 表达式')}</label>
              <input
                type="text"
                value={cronExpression}
                onChange={(e) => setCronExpression(e.target.value)}
                placeholder="0 */6 * * *"
                className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm font-mono text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-accent"
              />
              <p className="mt-1 text-[10px] text-muted-foreground">
                {t('格式：分 时 日 月 周，例如 "0 */6 * * *" 表示每6小时执行')}
              </p>
            </div>
          )}

          {scheduleMode === 'interval' && (
            <div>
              <label className="block text-xs font-medium text-foreground mb-1">{t('执行间隔')}</label>
              <div className="flex gap-2">
                <input
                  type="number"
                  min={1}
                  value={intervalValue}
                  onChange={(e) => setIntervalValue(Math.max(1, parseInt(e.target.value) || 1))}
                  className="flex-1 rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-accent"
                />
                <select
                  value={intervalUnit}
                  onChange={(e) => setIntervalUnit(e.target.value)}
                  className="rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-accent"
                >
                  {INTERVAL_UNITS.map((u) => (
                    <option key={u.id} value={u.id}>{t(u.label)}</option>
                  ))}
                </select>
              </div>
            </div>
          )}

          {scheduleMode === 'once' && (
            <div>
              <label className="block text-xs font-medium text-foreground mb-1">{t('执行时间')}</label>
              <input
                type="datetime-local"
                value={scheduledAt}
                onChange={(e) => setScheduledAt(e.target.value)}
                className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-sm text-foreground focus:outline-none focus:ring-1 focus:ring-accent"
              />
            </div>
          )}

          {/* 启用/禁用 */}
          <div className="flex items-center justify-between">
            <label className="text-xs font-medium text-foreground">{t('启用任务')}</label>
            <button
              type="button"
              onClick={() => setEnabled(!enabled)}
              className={`relative inline-flex h-5 w-9 items-center rounded-full transition-colors ${
                enabled ? 'bg-accent' : 'bg-muted'
              }`}
            >
              <span
                className={`inline-block h-3.5 w-3.5 rounded-full bg-white transition-transform ${
                  enabled ? 'translate-x-4.5' : 'translate-x-0.5'
                }`}
              />
            </button>
          </div>

          {/* 高级配置 (payload JSON) */}
          <div>
            <label className="block text-xs font-medium text-foreground mb-1">
              {t('高级配置')}
              <span className="ml-1 text-muted-foreground font-normal">{t('(JSON，可选)')}</span>
            </label>
            <textarea
              value={payload}
              onChange={(e) => setPayload(e.target.value)}
              placeholder='{"key": "value"}'
              rows={3}
              className="w-full rounded-md border border-border/60 bg-background px-3 py-1.5 text-xs font-mono text-foreground placeholder:text-muted-foreground/50 focus:outline-none focus:ring-1 focus:ring-accent resize-none"
            />
          </div>
        </div>

        {/* 底部按钮 */}
        <div className="flex items-center justify-end gap-2 px-6 py-3 border-t border-border/40">
          <button
            type="button"
            onClick={onClose}
            className="rounded-md px-3 py-1.5 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors"
          >
            {t('取消')}
          </button>
          <button
            type="button"
            onClick={handleSubmit}
            disabled={!name.trim() || saving}
            className="rounded-md bg-foreground px-4 py-1.5 text-xs font-medium text-background hover:bg-foreground/90 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {saving ? t('保存中...') : isEditing ? t('保存') : t('创建')}
          </button>
        </div>
      </div>
    </div>
  )
}
