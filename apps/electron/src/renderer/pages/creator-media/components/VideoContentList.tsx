/**
 * VideoContentList — 视频内容列表组件
 *
 * 从 contents 表筛选视频类型内容，展示在视频工作台左面板
 */

import { useT } from '@/context/LocaleContext'
import { Film, Video } from 'lucide-react'
import { cn } from '@/lib/utils'
import type { Content, ContentVideoMetadata } from '@sprouty-ai/shared/db/types'

/** 状态徽章颜色 */
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

const RENDER_STATUS_LABELS: Record<string, string> = {
  not_started: '未开始',
  rendering: '渲染中',
  completed: '已完成',
  failed: '失败',
}

const RENDER_STATUS_COLORS: Record<string, string> = {
  not_started: 'text-muted-foreground',
  rendering: 'text-amber-500',
  completed: 'text-green-500',
  failed: 'text-red-500',
}

/** 解析 metadata */
function parseVideoMeta(content: Content): ContentVideoMetadata | null {
  if (!content.metadata) return null
  try {
    return JSON.parse(content.metadata) as ContentVideoMetadata
  } catch {
    return null
  }
}

interface VideoContentListProps {
  contents: Content[]
  selected?: string
  onSelect: (content: Content) => void
  className?: string
}

export function VideoContentList({ contents, selected, onSelect, className }: VideoContentListProps) {
  const t = useT()

  if (contents.length === 0) {
    return (
      <div className={cn('flex flex-col items-center justify-center py-12 text-muted-foreground', className)}>
        <Film className="h-10 w-10 mb-3 opacity-30" />
        <p className="text-xs">{t('暂无视频内容')}</p>
        <p className="text-[10px] mt-1">{t('在创作工作台新建视频类型内容')}</p>
      </div>
    )
  }

  return (
    <div className={cn('p-1.5 space-y-1', className)}>
      {contents.map((content) => {
        const meta = parseVideoMeta(content)
        const isSelected = selected === content.id
        const isCreating = content.status === 'creating'
        const renderStatus = meta?.video_render_status
        const resolution = meta?.video_resolution

        return (
          <button
            key={content.id}
            type="button"
            onClick={() => onSelect(content)}
            className={cn(
              'w-full text-left rounded-md px-2.5 py-2 transition-colors',
              isSelected
                ? 'bg-accent text-accent-foreground'
                : 'hover:bg-muted/60',
              isCreating && !isSelected && 'ring-1 ring-orange-400/40',
            )}
          >
            {/* 标题行 */}
            <div className="flex items-center gap-1.5">
              <Video className="h-3.5 w-3.5 shrink-0 text-muted-foreground" />
              <span className="text-sm font-medium truncate flex-1">
                {content.title || t('无标题')}
              </span>
            </div>

            {/* 状态行 */}
            <div className="flex items-center gap-1.5 mt-1">
              <span className={cn(
                'inline-flex items-center rounded-full px-1.5 py-0 text-[10px] font-medium',
                STATUS_COLORS[content.status] || STATUS_COLORS.idea,
              )}>
                {t(STATUS_LABELS[content.status] || content.status)}
              </span>

              {renderStatus && (
                <span className={cn(
                  'text-[10px]',
                  RENDER_STATUS_COLORS[renderStatus] || 'text-muted-foreground',
                )}>
                  {t(RENDER_STATUS_LABELS[renderStatus] || renderStatus)}
                </span>
              )}
            </div>

            {/* 元信息行 */}
            <div className="flex items-center gap-2 mt-1 text-[10px] text-muted-foreground">
              {meta?.video_template_id && (
                <span className="truncate">{meta.video_template_id}</span>
              )}
              {resolution && (
                <span>{resolution.width}×{resolution.height}</span>
              )}
              {content.updated_at && (
                <span className="ml-auto shrink-0">
                  {new Date(content.updated_at).toLocaleDateString('zh-CN')}
                </span>
              )}
            </div>
          </button>
        )
      })}
    </div>
  )
}
