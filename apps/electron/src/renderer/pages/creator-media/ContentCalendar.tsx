import { useState, useMemo } from 'react'
import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { ContentTable } from './components/ContentTable'

/** 获取某月的天数 */
function getDaysInMonth(year: number, month: number) {
  return new Date(year, month + 1, 0).getDate()
}

/** 获取某月第一天是星期几（0=周日） */
function getFirstDayOfMonth(year: number, month: number) {
  return new Date(year, month, 1).getDay()
}

/** 格式化日期为 YYYY-MM-DD */
function formatDate(year: number, month: number, day: number) {
  return `${year}-${String(month + 1).padStart(2, '0')}-${String(day).padStart(2, '0')}`
}

const WEEKDAYS = ['日', '一', '二', '三', '四', '五', '六']

/** 状态颜色点 */
const STATUS_DOT_COLORS: Record<string, string> = {
  idea: 'bg-blue-400',
  researching: 'bg-cyan-400',
  scripting: 'bg-amber-400',
  creating: 'bg-orange-400',
  reviewing: 'bg-purple-400',
  scheduled: 'bg-indigo-400',
  published: 'bg-green-400',
  archived: 'bg-gray-400',
}

/**
 * 内容日历 — 日历视图展示创作计划和发布记录
 */
export default function ContentCalendar() {
  const t = useT()
  const {
    projects, activeProject, contents, publishRecords, loading,
    switchProject, deleteContent, updateContentStatus,
  } = useCreatorMedia()

  const now = new Date()
  const [year, setYear] = useState(now.getFullYear())
  const [month, setMonth] = useState(now.getMonth())
  const [selectedDate, setSelectedDate] = useState<string | null>(null)
  const [viewMode, setViewMode] = useState<'month' | 'week'>('month')

  /** 按日期分组内容 */
  const contentsByDate = useMemo(() => {
    const map: Record<string, typeof contents> = {}
    for (const c of contents) {
      // 优先使用 scheduled_at，其次 created_at
      const dateStr = c.scheduled_at || c.created_at
      if (!dateStr) continue
      const day = dateStr.slice(0, 10) // YYYY-MM-DD
      if (!map[day]) map[day] = []
      map[day].push(c)
    }
    // 也加入发布记录的日期
    for (const r of publishRecords) {
      if (!r.published_at) continue
      const day = r.published_at.slice(0, 10)
      if (!map[day]) map[day] = []
      // 发布记录不重复添加内容
    }
    return map
  }, [contents, publishRecords])

  /** 当前月份的日历网格 */
  const calendarDays = useMemo(() => {
    const daysInMonth = getDaysInMonth(year, month)
    const firstDay = getFirstDayOfMonth(year, month)
    const days: { day: number; date: string; isCurrentMonth: boolean }[] = []

    // 上月填充
    const prevMonthDays = getDaysInMonth(year, month - 1)
    for (let i = firstDay - 1; i >= 0; i--) {
      const d = prevMonthDays - i
      const m = month === 0 ? 11 : month - 1
      const y = month === 0 ? year - 1 : year
      days.push({ day: d, date: formatDate(y, m, d), isCurrentMonth: false })
    }

    // 当月
    for (let d = 1; d <= daysInMonth; d++) {
      days.push({ day: d, date: formatDate(year, month, d), isCurrentMonth: true })
    }

    // 下月填充到 42 天（6 行）
    const remaining = 42 - days.length
    for (let d = 1; d <= remaining; d++) {
      const m = month === 11 ? 0 : month + 1
      const y = month === 11 ? year + 1 : year
      days.push({ day: d, date: formatDate(y, m, d), isCurrentMonth: false })
    }

    return days
  }, [year, month])

  /** 选中日期的内容 */
  const selectedContents = useMemo(() => {
    if (!selectedDate) return []
    return contentsByDate[selectedDate] || []
  }, [selectedDate, contentsByDate])

  const todayStr = formatDate(now.getFullYear(), now.getMonth(), now.getDate())

  /** 月份导航 */
  const goToPrevMonth = () => {
    if (month === 0) { setYear(year - 1); setMonth(11) }
    else setMonth(month - 1)
  }
  const goToNextMonth = () => {
    if (month === 11) { setYear(year + 1); setMonth(0) }
    else setMonth(month + 1)
  }

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
        <div className="flex items-center gap-4">
          <div>
            <h1 className="text-base font-semibold text-foreground">{t('内容日历')}</h1>
            <p className="mt-0.5 text-xs text-muted-foreground">
              {activeProject ? activeProject.name : t('创作计划与发布记录')}
            </p>
          </div>
          {/* 月份导航 */}
          <div className="flex items-center gap-1">
            <button
              type="button"
              onClick={goToPrevMonth}
              className="inline-flex items-center justify-center rounded-md w-7 h-7 text-muted-foreground hover:text-foreground hover:bg-muted/40 transition-colors"
            >
              <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M15.75 19.5 8.25 12l7.5-7.5" />
              </svg>
            </button>
            <span className="text-sm font-medium text-foreground min-w-[100px] text-center">
              {year}{t('年')}{month + 1}{t('月')}
            </span>
            <button
              type="button"
              onClick={goToNextMonth}
              className="inline-flex items-center justify-center rounded-md w-7 h-7 text-muted-foreground hover:text-foreground hover:bg-muted/40 transition-colors"
            >
              <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                <path strokeLinecap="round" strokeLinejoin="round" d="m8.25 4.5 7.5 7.5-7.5 7.5" />
              </svg>
            </button>
          </div>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          {/* 视图切换 */}
          <div className="flex items-center rounded-md border border-border/60 overflow-hidden">
            <button
              type="button"
              onClick={() => setViewMode('month')}
              className={`px-2.5 py-1 text-xs transition-colors ${viewMode === 'month' ? 'bg-foreground text-background' : 'text-muted-foreground hover:bg-muted/40'}`}
            >
              {t('月')}
            </button>
            <button
              type="button"
              onClick={() => setViewMode('week')}
              className={`px-2.5 py-1 text-xs transition-colors ${viewMode === 'week' ? 'bg-foreground text-background' : 'text-muted-foreground hover:bg-muted/40'}`}
            >
              {t('周')}
            </button>
          </div>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-4 space-y-4">
        {/* 日历网格 */}
        <div className="rounded-lg border border-border/60 overflow-hidden">
          {/* 星期头 */}
          <div className="grid grid-cols-7 bg-muted/30 border-b border-border/40">
            {WEEKDAYS.map((d) => (
              <div key={d} className="px-2 py-2 text-center text-xs font-medium text-muted-foreground">
                {t(d)}
              </div>
            ))}
          </div>
          {/* 日期格子 */}
          <div className="grid grid-cols-7">
            {calendarDays.map((cell, idx) => {
              const dayContents = contentsByDate[cell.date] || []
              const isToday = cell.date === todayStr
              const isSelected = cell.date === selectedDate

              return (
                <button
                  key={idx}
                  type="button"
                  onClick={() => setSelectedDate(cell.date === selectedDate ? null : cell.date)}
                  className={`relative min-h-[72px] px-2 py-1.5 text-left border-b border-r border-border/20 transition-colors
                    ${!cell.isCurrentMonth ? 'bg-muted/10 text-muted-foreground/40' : 'hover:bg-muted/20'}
                    ${isSelected ? 'bg-blue-50 dark:bg-blue-900/20 ring-1 ring-blue-400 ring-inset' : ''}
                  `}
                >
                  <span className={`text-xs font-medium ${isToday ? 'inline-flex items-center justify-center w-5 h-5 rounded-full bg-foreground text-background' : ''}`}>
                    {cell.day}
                  </span>
                  {/* 内容标记 */}
                  {dayContents.length > 0 && (
                    <div className="mt-1 space-y-0.5">
                      {dayContents.slice(0, 3).map((c) => (
                        <div key={c.id} className="flex items-center gap-1">
                          <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${STATUS_DOT_COLORS[c.status] || 'bg-gray-400'}`} />
                          <span className="text-[10px] text-foreground/70 truncate">{c.title || t('无标题')}</span>
                        </div>
                      ))}
                      {dayContents.length > 3 && (
                        <span className="text-[10px] text-muted-foreground">+{dayContents.length - 3}</span>
                      )}
                    </div>
                  )}
                </button>
              )
            })}
          </div>
        </div>

        {/* 选中日期的详情 */}
        {selectedDate && (
          <div>
            <h2 className="text-sm font-medium text-foreground mb-2">
              {selectedDate} {t('的内容')} ({selectedContents.length})
            </h2>
            {selectedContents.length > 0 ? (
              <ContentTable
                contents={selectedContents}
                onStatusChange={updateContentStatus}
                onDelete={deleteContent}
              />
            ) : (
              <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-6 text-center">
                <p className="text-sm text-muted-foreground">{t('当日暂无内容')}</p>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  )
}
