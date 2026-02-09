import { useState, useMemo } from 'react'
import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { PublishRecordTable } from './components/PublishRecordTable'
import { ViralPatternList } from './components/ViralPatternList'
import { ViralPatternDialog } from './components/ViralPatternDialog'
import type { ViralPattern } from '@sprouty-ai/shared/db/types'

type TimeRange = 'week' | 'month' | 'all'

/**
 * 数据复盘 — 发布数据汇总 + 爆款模式库管理
 */
export default function AnalyticsDashboard() {
  const t = useT()
  const {
    projects, activeProject, publishRecords, viralPatterns, publishStats, loading,
    switchProject, createViralPattern, updateViralPattern, deleteViralPattern,
  } = useCreatorMedia()

  const [timeRange, setTimeRange] = useState<TimeRange>('all')
  const [showPatternDialog, setShowPatternDialog] = useState(false)
  const [editingPattern, setEditingPattern] = useState<ViralPattern | null>(null)

  /** 按时间范围筛选发布记录 */
  const filteredRecords = useMemo(() => {
    if (timeRange === 'all') return publishRecords
    const now = Date.now()
    const cutoff = timeRange === 'week' ? now - 7 * 86400000 : now - 30 * 86400000
    return publishRecords.filter((r) => {
      const ts = r.published_at ? new Date(r.published_at).getTime() : new Date(r.created_at).getTime()
      return ts >= cutoff
    })
  }, [publishRecords, timeRange])

  /** 筛选后的统计 */
  const filteredStats = useMemo(() => {
    const base = { totalViews: 0, totalLikes: 0, totalComments: 0, totalFavorites: 0 }
    for (const r of filteredRecords) {
      base.totalViews += r.views || 0
      base.totalLikes += r.likes || 0
      base.totalComments += r.comments || 0
      base.totalFavorites += r.favorites || 0
    }
    return base
  }, [filteredRecords])

  const handleAddPattern = () => {
    setEditingPattern(null)
    setShowPatternDialog(true)
  }

  const handleEditPattern = (pattern: ViralPattern) => {
    setEditingPattern(pattern)
    setShowPatternDialog(true)
  }

  const handleDeletePattern = async (id: string) => {
    await deleteViralPattern(id)
  }

  const handleSavePattern = async (data: Parameters<typeof createViralPattern>[0]) => {
    if (editingPattern) {
      await updateViralPattern(editingPattern.id, data)
    } else {
      await createViralPattern({ ...data, usage_count: 0, success_rate: null, tags: null, examples: null })
    }
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
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('数据复盘')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('发布数据与爆款模式')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          {/* 时间范围筛选 */}
          <div className="flex items-center rounded-md border border-border/60 overflow-hidden">
            {([
              { value: 'week' as TimeRange, label: '本周' },
              { value: 'month' as TimeRange, label: '本月' },
              { value: 'all' as TimeRange, label: '全部' },
            ]).map((item) => (
              <button
                key={item.value}
                type="button"
                onClick={() => setTimeRange(item.value)}
                className={`px-2.5 py-1 text-xs transition-colors ${timeRange === item.value ? 'bg-foreground text-background' : 'text-muted-foreground hover:bg-muted/40'}`}
              >
                {t(item.label)}
              </button>
            ))}
          </div>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-6">
        {/* 统计卡片行 */}
        <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
          <div className="flex items-center gap-3 rounded-lg border border-border/60 bg-background/40 px-4 py-3">
            <div className="flex h-9 w-9 items-center justify-center rounded-md bg-muted/50 text-blue-500">
              <svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M2 12s3-7 10-7 10 7 10 7-3 7-10 7-10-7-10-7Z" />
                <circle cx="12" cy="12" r="3" />
              </svg>
            </div>
            <div>
              <p className="text-xl font-semibold text-foreground tabular-nums">{filteredStats.totalViews.toLocaleString()}</p>
              <p className="text-xs text-muted-foreground">{t('总阅读')}</p>
            </div>
          </div>
          <div className="flex items-center gap-3 rounded-lg border border-border/60 bg-background/40 px-4 py-3">
            <div className="flex h-9 w-9 items-center justify-center rounded-md bg-muted/50 text-red-500">
              <svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M19 14c1.49-1.46 3-3.21 3-5.5A5.5 5.5 0 0 0 16.5 3c-1.76 0-3 .5-4.5 2-1.5-1.5-2.74-2-4.5-2A5.5 5.5 0 0 0 2 8.5c0 2.3 1.5 4.05 3 5.5l7 7Z" />
              </svg>
            </div>
            <div>
              <p className="text-xl font-semibold text-foreground tabular-nums">{filteredStats.totalLikes.toLocaleString()}</p>
              <p className="text-xs text-muted-foreground">{t('总点赞')}</p>
            </div>
          </div>
          <div className="flex items-center gap-3 rounded-lg border border-border/60 bg-background/40 px-4 py-3">
            <div className="flex h-9 w-9 items-center justify-center rounded-md bg-muted/50 text-amber-500">
              <svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M7.9 20A9 9 0 1 0 4 16.1L2 22Z" />
              </svg>
            </div>
            <div>
              <p className="text-xl font-semibold text-foreground tabular-nums">{filteredStats.totalComments.toLocaleString()}</p>
              <p className="text-xs text-muted-foreground">{t('总评论')}</p>
            </div>
          </div>
          <div className="flex items-center gap-3 rounded-lg border border-border/60 bg-background/40 px-4 py-3">
            <div className="flex h-9 w-9 items-center justify-center rounded-md bg-muted/50 text-purple-500">
              <svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="m19 21-7-4-7 4V5a2 2 0 0 1 2-2h10a2 2 0 0 1 2 2v16z" />
              </svg>
            </div>
            <div>
              <p className="text-xl font-semibold text-foreground tabular-nums">{filteredStats.totalFavorites.toLocaleString()}</p>
              <p className="text-xs text-muted-foreground">{t('总收藏')}</p>
            </div>
          </div>
        </div>

        {/* 内容排行区 */}
        <div>
          <h2 className="text-sm font-medium text-foreground mb-3">{t('发布记录排行')}</h2>
          <PublishRecordTable records={filteredRecords} />
        </div>

        {/* 爆款模式库区 */}
        <ViralPatternList
          patterns={viralPatterns}
          onEdit={handleEditPattern}
          onDelete={handleDeletePattern}
          onAdd={handleAddPattern}
        />
      </div>

      {/* 爆款模式对话框 */}
      <ViralPatternDialog
        open={showPatternDialog}
        onOpenChange={setShowPatternDialog}
        pattern={editingPattern}
        projectId={activeProject?.id || null}
        onSave={handleSavePattern}
      />
    </div>
  )
}
