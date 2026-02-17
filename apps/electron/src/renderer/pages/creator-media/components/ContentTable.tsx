import { useT } from '@/context/LocaleContext'
import type { Content, ContentStageRecord } from '@sprouty-ai/shared/db/types'
import { PLATFORM_MAP } from '@sprouty-ai/shared/db/types'
import { Video, FileText, ChevronDown, ChevronRight } from 'lucide-react'
import { useState, useEffect } from 'react'
import { useActiveWorkspace } from '@/context/AppShellContext'

/** Stage 类型标签映射 */
const STAGE_TYPE_LABELS: Record<string, string> = {
  'topic_recommend': '选题推荐',
  'research': '灵感调研',
  'script_article': '图文脚本',
  'script_video': '视频脚本',
  'draft_article': '图文原稿',
  'draft_video': '视频项目',
  'platform_adapt_article': '图文平台适配',
  'platform_adapt_video': '视频平台适配',
}

/** Stage 类型图标映射 */
const STAGE_TYPE_ICONS: Record<string, any> = {
  'script_article': FileText,
  'draft_article': FileText,
  'platform_adapt_article': FileText,
  'script_video': Video,
  'draft_video': Video,
  'platform_adapt_video': Video,
}

/** Stage 操作按钮配置 */
const STAGE_ACTIONS: Record<string, { skill: string; label: string } | null> = {
  'script_article': { skill: 'content-creator', label: '图文创作' },
  'script_video': { skill: 'video-creator', label: '视频创作' },
  'draft_article': { skill: 'platform-adapter', label: '平台适配' },
  'draft_video': { skill: 'platform-adapter', label: '平台适配' },
  'platform_adapt_article': null,
  'platform_adapt_video': null,
}

/** 状态徽章颜色映射 */
const STATUS_COLORS: Record<string, string> = {
  researching: 'bg-cyan-100 text-cyan-700 dark:bg-cyan-900/30 dark:text-cyan-400',
  scripting: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  creating: 'bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-400',
  adapting: 'bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-400',
  scheduled: 'bg-indigo-100 text-indigo-700 dark:bg-indigo-900/30 dark:text-indigo-400',
  published: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  archived: 'bg-gray-100 text-gray-700 dark:bg-gray-900/30 dark:text-gray-400',
}

/** 状态中文名映射 */
const STATUS_LABELS: Record<string, string> = {
  researching: '调研中',
  scripting: '脚本创作中',
  creating: '内容创作中',
  adapting: '平台适配中',
  scheduled: '已排期',
  published: '已发布',
  archived: '已归档',
}

/** 状态流转：每个状态可以转到的下一个状态 */
const STATUS_TRANSITIONS: Record<string, string[]> = {
  researching: ['scripting', 'archived'],
  scripting: ['creating', 'archived'],
  creating: ['adapting', 'archived'],
  adapting: ['scheduled', 'creating'],
  scheduled: ['published'],
  published: ['archived'],
  archived: ['researching'],
}

/** 下一阶段操作：每个状态对应的 skill 和操作文案 */
const NEXT_STAGE_ACTIONS: Record<string, { skill: string; label: string; icon: string } | null> = {
  researching: { skill: 'idea-researcher', label: '灵感调研', icon: '🔍' },
  scripting: null, // scripting 状态显示两个按钮：图文脚本、视频脚本
  creating: null, // creating 状态下的操作在 content_stages 子表中显示
  adapting: { skill: 'platform-adapter', label: '平台适配', icon: '📱' },
  scheduled: null,
  published: null,
  archived: null,
}

/** scripting 状态的脚本创建按钮配置 */
const SCRIPTING_ACTIONS = [
  { skill: 'article-script-create', label: '图文脚本', icon: FileText },
  { skill: 'video-script-create', label: '视频脚本', icon: Video },
]

/** 解析 target_platforms JSON 字符串为平台 id 数组 */
function parsePlatforms(raw: string | null): string[] {
  if (!raw) return []
  try {
    const parsed = JSON.parse(raw)
    return Array.isArray(parsed) ? parsed : []
  } catch {
    // 兼容逗号分隔的旧格式
    return raw.split(',').map(s => s.trim()).filter(Boolean)
  }
}

/** 平台标签组件 */
function PlatformTags({ platforms }: { platforms: string }) {
  const t = useT()
  const ids = parsePlatforms(platforms)
  if (ids.length === 0) return <span className="text-muted-foreground">-</span>
  return (
    <div className="flex flex-wrap gap-1">
      {ids.map((id) => {
        const meta = PLATFORM_MAP[id as keyof typeof PLATFORM_MAP]
        return (
          <span
            key={id}
            className={`inline-flex items-center rounded-full border px-1.5 py-0 text-[10px] font-medium whitespace-nowrap ${meta?.color || 'text-muted-foreground border-border'}`}
          >
            {meta?.shortLabel || id}
          </span>
        )
      })}
    </div>
  )
}

function StatusBadge({ status }: { status: string }) {
  const t = useT()
  const color = STATUS_COLORS[status] || STATUS_COLORS.idea
  const label = STATUS_LABELS[status] || status
  return (
    <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-xs font-medium whitespace-nowrap ${color}`}>
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
  onNextStage?: (content: Content, skillId: string) => void
  onScriptAction?: (content: Content, skillId: string) => void
  onStageAction?: (content: Content, stage: ContentStageRecord, skillId: string) => void
}

export function ContentTable({
  contents,
  maxItems = 10,
  onRowClick,
  onStatusChange,
  onDelete,
  onVersionHistory,
  onOpenVideoStudio,
  onNextStage,
  onScriptAction,
  onStageAction
}: ContentTableProps) {
  const t = useT()
  const workspace = useActiveWorkspace()
  const items = maxItems ? contents.slice(0, maxItems) : contents

  // 存储每个内容的展开状态
  const [expandedRows, setExpandedRows] = useState<Set<string>>(new Set())

  // 存储每个内容的 stages
  const [contentStages, setContentStages] = useState<Record<string, ContentStageRecord[]>>({})

  // 加载内容的 stages
  useEffect(() => {
    if (!workspace) return

    const loadStages = async () => {
      const stages: Record<string, ContentStageRecord[]> = {}

      for (const item of items) {
        try {
          const itemStages = await window.electronAPI.creatorMedia.contentStages.list(workspace.id, item.id)
          stages[item.id] = itemStages
        } catch {
          stages[item.id] = []
        }
      }

      setContentStages(stages)
    }

    loadStages()
  }, [workspace, items])

  const toggleRow = (contentId: string) => {
    setExpandedRows(prev => {
      const next = new Set(prev)
      if (next.has(contentId)) {
        next.delete(contentId)
      } else {
        next.add(contentId)
      }
      return next
    })
  }

  if (items.length === 0) {
    return (
      <div className="flex items-center justify-center py-8 text-muted-foreground">
        <p className="text-sm">{t('暂无内容，开始创作吧')}</p>
      </div>
    )
  }

  return (
    <div className="overflow-hidden rounded-lg border border-border/60">
      <div className="overflow-x-auto">
      <table className="w-full text-sm table-fixed">
        <thead>
          <tr className="border-b border-border/40 bg-muted/30">
            <th className="px-3 py-2.5 text-left font-medium text-muted-foreground w-9"></th>
            <th className="px-3 py-2.5 text-left font-medium text-muted-foreground w-[30%]">{t('标题')}</th>
            <th className="px-3 py-2.5 text-left font-medium text-muted-foreground w-[12%]">{t('状态')}</th>
            <th className="px-3 py-2.5 text-left font-medium text-muted-foreground w-[18%]">{t('平台')}</th>
            <th className="px-3 py-2.5 text-left font-medium text-muted-foreground w-[10%]">{t('更新时间')}</th>
            {(onStatusChange || onDelete || onVersionHistory || onNextStage) && (
              <th className="px-3 py-2.5 text-right font-medium text-muted-foreground">{t('操作')}</th>
            )}
          </tr>
        </thead>
        <tbody>
          {items.map((item) => {
            const nextStatuses = STATUS_TRANSITIONS[item.status] || []
            const nextStageAction = NEXT_STAGE_ACTIONS[item.status]
            const stages = contentStages[item.id] || []
            const isExpanded = expandedRows.has(item.id)
            const hasStages = stages.length > 0

            return (
              <>
                {/* 主行：Content */}
                <tr
                  key={item.id}
                  className={`border-b border-border/20 hover:bg-muted/20 transition-colors ${hasStages ? 'cursor-pointer' : ''}`}
                  onClick={() => hasStages && toggleRow(item.id)}
                >
                  <td className="px-3 py-3">
                    {hasStages && (
                      <button
                        type="button"
                        onClick={(e) => {
                          e.stopPropagation()
                          toggleRow(item.id)
                        }}
                        className="text-muted-foreground hover:text-foreground"
                      >
                        {isExpanded ? <ChevronDown className="h-4 w-4" /> : <ChevronRight className="h-4 w-4" />}
                      </button>
                    )}
                  </td>
                  <td className="px-3 py-3 font-medium text-foreground">
                    <div className="min-w-0">
                      <span className="block truncate">{item.title || t('无标题')}</span>
                    </div>
                  </td>
                  <td className="px-3 py-3">
                    <StatusBadge status={item.status} />
                  </td>
                  <td className="px-3 py-3 text-xs"><PlatformTags platforms={item.target_platforms || ''} /></td>
                  <td className="px-3 py-3 text-muted-foreground text-xs whitespace-nowrap">{item.updated_at ? new Date(item.updated_at).toLocaleDateString('zh-CN') : '-'}</td>
                  {(onStatusChange || onDelete || onVersionHistory || onNextStage || onScriptAction) && (
                    <td className="px-3 py-3 text-right" onClick={(e) => e.stopPropagation()}>
                      <div className="flex items-center justify-end gap-1 flex-wrap">
                        {/* scripting 状态显示两个脚本创建按钮 */}
                        {onScriptAction && item.status === 'scripting' && (
                          SCRIPTING_ACTIONS.map((action) => {
                            const Icon = action.icon
                            return (
                              <button
                                key={action.skill}
                                type="button"
                                onClick={() => onScriptAction(item, action.skill)}
                                className="inline-flex items-center gap-1 rounded-md px-2 py-1 text-xs font-medium text-white bg-blue-500 hover:bg-blue-600 transition-colors whitespace-nowrap"
                              >
                                <Icon className="h-3 w-3" />
                                <span>{t(action.label)}</span>
                              </button>
                            )
                          })
                        )}
                        {/* 其他状态的下一阶段操作按钮 */}
                        {onNextStage && nextStageAction && (
                          <button
                            type="button"
                            onClick={() => onNextStage(item, nextStageAction.skill)}
                            className="inline-flex items-center gap-1 rounded-md px-2 py-1 text-xs font-medium text-white bg-blue-500 hover:bg-blue-600 transition-colors whitespace-nowrap"
                          >
                            <span>{nextStageAction.icon}</span>
                            <span>{t(nextStageAction.label)}</span>
                          </button>
                        )}
                        {onStatusChange && nextStatuses.length > 0 && (
                          nextStatuses.map((ns) => (
                            <button
                              key={ns}
                              type="button"
                              onClick={() => onStatusChange(item.id, ns)}
                              className="inline-flex items-center rounded px-2 py-1 text-xs font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors whitespace-nowrap"
                            >
                              {ns === 'archived' ? t('归档') : `→ ${t(STATUS_LABELS[ns] || ns)}`}
                            </button>
                          ))
                        )}
                        {onVersionHistory && (
                          <button
                            type="button"
                            onClick={() => onVersionHistory(item)}
                            className="inline-flex items-center rounded px-2 py-1 text-xs font-medium text-muted-foreground hover:text-foreground hover:bg-muted/60 transition-colors whitespace-nowrap"
                          >
                            {t('版本')}
                          </button>
                        )}
                        {onDelete && (
                          <button
                            type="button"
                            onClick={() => onDelete(item.id)}
                            className="inline-flex items-center rounded px-2 py-1 text-xs font-medium text-red-500 hover:text-red-600 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors whitespace-nowrap"
                          >
                            {t('删除')}
                          </button>
                        )}
                      </div>
                    </td>
                  )}
                </tr>

                {/* 子行：ContentStages */}
                {isExpanded && stages.map((stage) => {
                  const Icon = STAGE_TYPE_ICONS[stage.stage]
                  const stageAction = STAGE_ACTIONS[stage.stage]

                  return (
                    <tr
                      key={stage.id}
                      className="border-b border-border/10 bg-muted/10 hover:bg-muted/20 transition-colors"
                    >
                      <td className="px-3 py-2"></td>
                      <td className="px-3 py-2 text-sm text-muted-foreground">
                        <div className="flex items-center gap-2 pl-4">
                          {Icon && <Icon className="h-3.5 w-3.5" />}
                          <span>{t(STAGE_TYPE_LABELS[stage.stage] || stage.stage)}</span>
                        </div>
                      </td>
                      <td className="px-3 py-2 text-xs text-muted-foreground"></td>
                      <td className="px-3 py-2 text-xs text-muted-foreground">
                        {stage.file_path ? (
                          <span className="block truncate" title={stage.file_path}>
                            {stage.file_path.split('/').pop()}
                          </span>
                        ) : '-'}
                      </td>
                      <td className="px-3 py-2 text-xs text-muted-foreground whitespace-nowrap">
                        {stage.updated_at ? new Date(stage.updated_at).toLocaleDateString('zh-CN') : '-'}
                      </td>
                      <td className="px-3 py-2 text-right">
                        <div className="flex items-center justify-end gap-1">
                          {/* Stage 操作按钮 */}
                          {onStageAction && stageAction && (
                            <button
                              type="button"
                              onClick={() => onStageAction(item, stage, stageAction.skill)}
                              className="inline-flex items-center gap-1 rounded-md px-2 py-1 text-xs font-medium text-white bg-green-500 hover:bg-green-600 transition-colors whitespace-nowrap"
                            >
                              <span>{t(stageAction.label)}</span>
                            </button>
                          )}
                          {/* 视频工作台按钮 */}
                          {onOpenVideoStudio && stage.stage === 'draft_video' && (
                            <button
                              type="button"
                              onClick={() => onOpenVideoStudio(item)}
                              className="inline-flex items-center gap-0.5 rounded px-2 py-1 text-xs font-medium text-blue-600 hover:text-blue-700 hover:bg-blue-50 dark:text-blue-400 dark:hover:bg-blue-900/20 transition-colors whitespace-nowrap"
                            >
                              <Video className="h-3 w-3" />
                              <span>{t('视频工作台')}</span>
                            </button>
                          )}
                        </div>
                      </td>
                    </tr>
                  )
                })}
              </>
            )
          })}
        </tbody>
      </table>
      </div>
    </div>
  )
}
