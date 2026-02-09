import { useT } from '@/context/LocaleContext'
import type { Content } from '@sprouty-ai/shared/db/types'

/** 状态徽章颜色映射 */
const STATUS_COLORS: Record<string, string> = {
  idea: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400',
  researching: 'bg-cyan-100 text-cyan-700 dark:bg-cyan-900/30 dark:text-cyan-400',
  scripting: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  creating: 'bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-400',
  reviewing: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400',
  scheduled: 'bg-indigo-100 text-indigo-700 dark:bg-indigo-900/30 dark:text-indigo-400',
  published: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  archived: 'bg-gray-100 text-gray-700 dark:bg-gray-900/30 dark:text-gray-400',
}

/** 状态中文名映射 */
const STATUS_LABELS: Record<string, string> = {
  idea: '选题',
  researching: '研究中',
  scripting: '写脚本',
  creating: '创作中',
  reviewing: '审核中',
  scheduled: '待发布',
  published: '已发布',
  archived: '已归档',
}

/** 状态流转：每个状态可以转到的下一个状态 */
const STATUS_TRANSITIONS: Record<string, string[]> = {
  idea: ['researching', 'archived'],
  researching: ['scripting', 'archived'],
  scripting: ['creating', 'archived'],
  creating: ['reviewing', 'archived'],
  reviewing: ['scheduled', 'creating'],
  scheduled: ['published'],
  published: ['archived'],
  archived: ['idea'],
}

function StatusBadge({ status }: { status: string }) {
  const t = useT()
  const color = STATUS_COLORS[status] || STATUS_COLORS.idea
  const label = STATUS_LABELS[status] || status
  return (
    <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-xs font-medium ${color}`}>
      {t(label)}
    </span>
  )
}

interface ContentTableProps {
  contents: Content[]
  maxItems?: number
  onRowClick?: (content: Content) => void
  onStatusChange?: (contentId: string, status: string) => void
  onDelete?: (contentId: string) => void
}

export function ContentTable({ contents, maxItems = 10, onRowClick, onStatusChange, onDelete }: ContentTableProps) {
  const t = useT()
  const items = maxItems ? contents.slice(0, maxItems) : contents

  if (items.length === 0) {
    return (
      <div className="flex items-center justify-center py-8 text-muted-foreground">
        <p className="text-sm">{t('暂无内容，开始创作吧')}</p>
      </div>
    )
  }

  return (
    <div className="overflow-hidden rounded-lg border border-border/60">
      <table className="w-full text-sm">
        <thead>
          <tr className="border-b border-border/40 bg-muted/30">
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('标题')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('状态')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('平台')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('更新时间')}</th>
            {(onStatusChange || onDelete) && (
              <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('操作')}</th>
            )}
          </tr>
        </thead>
        <tbody>
          {items.map((item) => {
            const nextStatuses = STATUS_TRANSITIONS[item.status] || []
            return (
              <tr
                key={item.id}
                className={`border-b border-border/20 last:border-0 hover:bg-muted/20 transition-colors ${onRowClick ? 'cursor-pointer' : ''}`}
                onClick={() => onRowClick?.(item)}
              >
                <td className="px-4 py-2.5 font-medium text-foreground truncate max-w-[240px]">{item.title || t('无标题')}</td>
                <td className="px-4 py-2.5"><StatusBadge status={item.status} /></td>
                <td className="px-4 py-2.5 text-muted-foreground">{item.target_platforms || '-'}</td>
                <td className="px-4 py-2.5 text-muted-foreground text-xs">{item.updated_at ? new Date(item.updated_at).toLocaleDateString('zh-CN') : '-'}</td>
                {(onStatusChange || onDelete) && (
                  <td className="px-4 py-2.5 text-right" onClick={(e) => e.stopPropagation()}>
                    <div className="flex items-center justify-end gap-1">
                      {onStatusChange && nextStatuses.length > 0 && (
                        nextStatuses.map((ns) => (
                          <button
                            key={ns}
                            type="button"
                            onClick={() => onStatusChange(item.id, ns)}
                            className="inline-flex items-center rounded px-1.5 py-0.5 text-[10px] font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors"
                          >
                            {ns === 'archived' ? t('归档') : `→ ${t(STATUS_LABELS[ns] || ns)}`}
                          </button>
                        ))
                      )}
                      {onDelete && (
                        <button
                          type="button"
                          onClick={() => onDelete(item.id)}
                          className="inline-flex items-center rounded px-1.5 py-0.5 text-[10px] font-medium text-red-500 hover:text-red-600 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors"
                        >
                          {t('删除')}
                        </button>
                      )}
                    </div>
                  </td>
                )}
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}
