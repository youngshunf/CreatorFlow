import { useT } from '@/context/LocaleContext'
import type { Content } from '@sprouty-ai/shared/db/types'
import { Video } from 'lucide-react'

/** 内容类型标签映射 */
const CONTENT_TYPE_LABELS: Record<string, string> = {
  'image-text': '图文',
  'video': '视频',
  'short-video': '短视频',
  'article': '文章',
  'live': '直播',
}

/** 从 metadata JSON 解析视频渲染状态 */
function getVideoRenderStatus(content: Content): string | null {
  if (!content.metadata) return null
  try {
    const meta = JSON.parse(content.metadata)
    return meta.video_render_status || null
  } catch {
    return null
  }
}

/** 渲染状态中文名映射 */
const RENDER_STATUS_LABELS: Record<string, string> = {
  not_started: '未渲染',
  rendering: '渲染中',
  completed: '已渲染',
  failed: '渲染失败',
}

/** 渲染状态颜色映射 */
const RENDER_STATUS_COLORS: Record<string, string> = {
  not_started: 'bg-gray-100 text-gray-600 dark:bg-gray-800/30 dark:text-gray-400',
  rendering: 'bg-amber-100 text-amber-600 dark:bg-amber-900/30 dark:text-amber-400',
  completed: 'bg-green-100 text-green-600 dark:bg-green-900/30 dark:text-green-400',
  failed: 'bg-red-100 text-red-600 dark:bg-red-900/30 dark:text-red-400',
}

/** 判断是否为视频类型内容 */
function isVideoContent(content: Content): boolean {
  return content.content_type === 'video' || content.content_type === 'short-video'
}

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
  onVersionHistory?: (content: Content) => void
  onOpenVideoStudio?: (content: Content) => void
}

export function ContentTable({ contents, maxItems = 10, onRowClick, onStatusChange, onDelete, onVersionHistory, onOpenVideoStudio }: ContentTableProps) {
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
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('类型')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('平台')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('更新时间')}</th>
            {(onStatusChange || onDelete || onVersionHistory || onOpenVideoStudio) && (
              <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('操作')}</th>
            )}
          </tr>
        </thead>
        <tbody>
          {items.map((item) => {
            const nextStatuses = STATUS_TRANSITIONS[item.status] || []
            const isVideo = isVideoContent(item)
            const renderStatus = isVideo ? getVideoRenderStatus(item) : null
            return (
              <tr
                key={item.id}
                className={`border-b border-border/20 last:border-0 hover:bg-muted/20 transition-colors ${onRowClick ? 'cursor-pointer' : ''}`}
                onClick={() => onRowClick?.(item)}
              >
                <td className="px-4 py-2.5 font-medium text-foreground truncate max-w-[240px]">
                  {isVideo && <Video className="h-3.5 w-3.5 inline mr-1 text-muted-foreground" />}
                  {item.title || t('无标题')}
                </td>
                <td className="px-4 py-2.5">
                  <div className="flex items-center gap-1.5">
                    <StatusBadge status={item.status} />
                    {/* 视频内容在 creating 状态时显示渲染状态 */}
                    {isVideo && item.status === 'creating' && renderStatus && (
                      <span className={`inline-flex items-center rounded-full px-1.5 py-0.5 text-[10px] font-medium ${RENDER_STATUS_COLORS[renderStatus] || ''}`}>
                        {t(RENDER_STATUS_LABELS[renderStatus] || renderStatus)}
                      </span>
                    )}
                  </div>
                </td>
                <td className="px-4 py-2.5 text-muted-foreground text-xs">
                  {t(CONTENT_TYPE_LABELS[item.content_type] || item.content_type || '-')}
                </td>
                <td className="px-4 py-2.5 text-muted-foreground">{item.target_platforms || '-'}</td>
                <td className="px-4 py-2.5 text-muted-foreground text-xs">{item.updated_at ? new Date(item.updated_at).toLocaleDateString('zh-CN') : '-'}</td>
                {(onStatusChange || onDelete || onVersionHistory || onOpenVideoStudio) && (
                  <td className="px-4 py-2.5 text-right" onClick={(e) => e.stopPropagation()}>
                    <div className="flex items-center justify-end gap-1">
                      {/* 视频内容显示"视频工作台"按钮 */}
                      {onOpenVideoStudio && isVideo && (
                        <button
                          type="button"
                          onClick={() => onOpenVideoStudio(item)}
                          className="inline-flex items-center rounded px-1.5 py-0.5 text-[10px] font-medium text-blue-500 hover:text-blue-600 hover:bg-blue-50 dark:hover:bg-blue-900/20 transition-colors"
                        >
                          <Video className="h-3 w-3 mr-0.5" />
                          {t('视频工作台')}
                        </button>
                      )}
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
                      {onVersionHistory && (
                        <button
                          type="button"
                          onClick={() => onVersionHistory(item)}
                          className="inline-flex items-center rounded px-1.5 py-0.5 text-[10px] font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors"
                        >
                          {t('版本')}
                        </button>
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
