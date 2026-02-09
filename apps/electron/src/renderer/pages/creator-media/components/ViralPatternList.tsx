import { useT } from '@/context/LocaleContext'
import type { ViralPattern } from '@sprouty-ai/shared/db/types'

/** 分类颜色 */
const CATEGORY_COLORS: Record<string, string> = {
  hook: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
  structure: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400',
  title: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  cta: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  visual: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400',
  rhythm: 'bg-cyan-100 text-cyan-700 dark:bg-cyan-900/30 dark:text-cyan-400',
}

const CATEGORY_LABELS: Record<string, string> = {
  hook: '开头钩子',
  structure: '内容结构',
  title: '标题模式',
  cta: '行动号召',
  visual: '视觉风格',
  rhythm: '节奏韵律',
}

const SOURCE_LABELS: Record<string, string> = {
  competitor_analysis: '竞品分析',
  manual: '手动添加',
  ai_discovered: 'AI 发现',
}

interface ViralPatternListProps {
  patterns: ViralPattern[]
  onEdit: (pattern: ViralPattern) => void
  onDelete: (id: string) => void
  onAdd: () => void
}

/**
 * 爆款模式列表 — 展示模式库 + CRUD 操作
 */
export function ViralPatternList({ patterns, onEdit, onDelete, onAdd }: ViralPatternListProps) {
  const t = useT()

  return (
    <div>
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-sm font-medium text-foreground">{t('爆款模式库')}</h2>
        <button
          type="button"
          onClick={onAdd}
          className="inline-flex items-center gap-1 rounded-md text-xs text-muted-foreground hover:text-foreground transition-colors"
        >
          <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
          </svg>
          {t('新增模式')}
        </button>
      </div>

      {patterns.length === 0 ? (
        <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-6 text-center">
          <p className="text-sm text-muted-foreground">{t('暂无爆款模式')}</p>
          <p className="mt-1 text-xs text-muted-foreground/70">{t('通过竞品分析或手动添加积累爆款模式')}</p>
        </div>
      ) : (
        <div className="space-y-2">
          {patterns.map((p) => (
            <div
              key={p.id}
              className="rounded-lg border border-border/60 bg-background/40 px-4 py-3 hover:bg-muted/20 transition-colors"
            >
              <div className="flex items-start justify-between gap-3">
                <div className="flex-1 min-w-0">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="text-sm font-medium text-foreground">{p.name}</span>
                    <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-medium ${CATEGORY_COLORS[p.category] || 'bg-gray-100 text-gray-700'}`}>
                      {t(CATEGORY_LABELS[p.category] || p.category)}
                    </span>
                    {p.platform && (
                      <span className="text-[10px] text-muted-foreground">{p.platform}</span>
                    )}
                  </div>
                  {p.description && (
                    <p className="text-xs text-muted-foreground line-clamp-2">{p.description}</p>
                  )}
                  <div className="flex items-center gap-4 mt-1.5 text-[10px] text-muted-foreground">
                    <span>{t('使用')} {p.usage_count} {t('次')}</span>
                    {p.success_rate != null && (
                      <span>{t('成功率')} {Math.round(p.success_rate * 100)}%</span>
                    )}
                    {p.source && (
                      <span>{t(SOURCE_LABELS[p.source] || p.source)}</span>
                    )}
                  </div>
                </div>
                <div className="flex items-center gap-1 flex-shrink-0">
                  <button
                    type="button"
                    onClick={() => onEdit(p)}
                    className="inline-flex items-center justify-center rounded-md w-7 h-7 text-muted-foreground hover:text-foreground hover:bg-muted/40 transition-colors"
                    title={t('编辑')}
                  >
                    <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L6.832 19.82a4.5 4.5 0 0 1-1.897 1.13l-2.685.8.8-2.685a4.5 4.5 0 0 1 1.13-1.897L16.863 4.487Z" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => onDelete(p.id)}
                    className="inline-flex items-center justify-center rounded-md w-7 h-7 text-red-500 hover:text-red-600 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors"
                    title={t('删除')}
                  >
                    <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="m14.74 9-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 0 1-2.244 2.077H8.084a2.25 2.25 0 0 1-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 0 0-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 0 1 3.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 0 0-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 0 0-7.5 0" />
                    </svg>
                  </button>
                </div>
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}
