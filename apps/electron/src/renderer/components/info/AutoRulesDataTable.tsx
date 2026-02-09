/**
 * AutoRulesDataTable
 *
 * Flat data table displaying all auto-label rules across all labels.
 * Each row shows which label a rule belongs to, the regex pattern, flags,
 * value template, and description.
 *
 * Rules are collected by recursively traversing the label tree and flattening
 * all autoRules into a single list.
 */

import * as React from 'react'
import { useState, useMemo } from 'react'
import type { ColumnDef } from '@tanstack/react-table'
import { Maximize2 } from 'lucide-react'
import { Info_DataTable, SortableHeader } from './Info_DataTable'
import { Info_Badge } from './Info_Badge'
import { Tooltip, TooltipTrigger, TooltipContent } from '@sprouty-ai/ui'
import { DataTableOverlay } from '@sprouty-ai/ui'
import { LabelIcon } from '@/components/ui/label-icon'
import { cn } from '@/lib/utils'
import { toast } from 'sonner'
import type { LabelConfig, AutoLabelRule } from '@sprouty-ai/shared/labels'
import { useT } from '@/context/LocaleContext'

/**
 * Flattened auto-rule row: associates a rule with its parent label
 */
interface AutoRuleRow {
  /** The label this rule belongs to */
  label: LabelConfig
  /** The auto-label rule */
  rule: AutoLabelRule
}

interface AutoRulesDataTableProps {
  /** Label tree (root-level nodes with nested children) */
  data: LabelConfig[]
  /** Show search input */
  searchable?: boolean
  /** Max height with scroll */
  maxHeight?: number
  /** Enable fullscreen button */
  fullscreen?: boolean
  /** Title for fullscreen overlay */
  fullscreenTitle?: string
  className?: string
}

/**
 * PatternBadge - Monospace regex pattern with click-to-copy and tooltip.
 * Mirrors the PatternBadge from PermissionsDataTable for consistency.
 */
function PatternBadge({ pattern, t }: { pattern: string; t: (key: string) => string }) {
  const handleClick = async () => {
    try {
      await navigator.clipboard.writeText(pattern)
      toast.success(t('已复制到剪贴板'))
    } catch {
      toast.error(t('复制失败'))
    }
  }

  const badge = (
    <button type="button" onClick={handleClick} className="text-left">
      <Info_Badge color="muted" className="font-mono select-none">
        <span className="block overflow-hidden whitespace-nowrap text-ellipsis max-w-[200px]">
          {pattern}
        </span>
      </Info_Badge>
    </button>
  )

  if (pattern.length >= 25) {
    return (
      <Tooltip>
        <TooltipTrigger asChild>{badge}</TooltipTrigger>
        <TooltipContent className="font-mono max-w-md break-all">{pattern}</TooltipContent>
      </Tooltip>
    )
  }

  return badge
}

// Column definitions for the auto-rules flat table - created dynamically to support i18n
function createColumns(t: (key: string) => string): ColumnDef<AutoRuleRow>[] {
  return [
    {
      id: 'label',
      header: ({ column }) => <SortableHeader column={column} title={t('标签')} />,
      accessorFn: (row) => row.label.name,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5 flex items-center gap-1.5">
          <LabelIcon label={row.original.label} size="xs" />
          <span className="text-sm truncate">{row.original.label.name}</span>
        </div>
      ),
      minSize: 100,
    },
    {
      id: 'pattern',
      header: ({ column }) => <SortableHeader column={column} title={t('匹配模式')} />,
      accessorFn: (row) => row.rule.pattern,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <PatternBadge pattern={row.original.rule.pattern} t={t} />
        </div>
      ),
      minSize: 120,
    },
    {
      id: 'flags',
      header: () => <span className="p-1.5 pl-2.5">{t('标志')}</span>,
      accessorFn: (row) => row.rule.flags ?? 'gi',
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <span className="text-xs text-muted-foreground font-mono">
            {row.original.rule.flags ?? 'gi'}
          </span>
        </div>
      ),
      minSize: 50,
    },
    {
      id: 'template',
      header: () => <span className="p-1.5 pl-2.5">{t('模板')}</span>,
      accessorFn: (row) => row.rule.valueTemplate ?? '',
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          {row.original.rule.valueTemplate ? (
            <Info_Badge color="muted" className="font-mono whitespace-nowrap">
              {row.original.rule.valueTemplate}
            </Info_Badge>
          ) : (
            <span className="text-muted-foreground/50 text-sm">—</span>
          )}
        </div>
      ),
      minSize: 80,
    },
    {
      id: 'description',
      header: () => <span className="p-1.5 pl-2.5">{t('描述')}</span>,
      accessorFn: (row) => row.rule.description ?? '',
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5 min-w-0">
          <span className="truncate block text-sm">
            {row.original.rule.description || '—'}
          </span>
        </div>
      ),
      meta: { fillWidth: true, truncate: true },
    },
  ]
}

/**
 * Recursively collect all auto-rules from the label tree,
 * associating each rule with its parent label.
 */
function collectAutoRules(labels: LabelConfig[]): AutoRuleRow[] {
  const rows: AutoRuleRow[] = []

  function traverse(nodes: LabelConfig[]) {
    for (const label of nodes) {
      if (label.autoRules?.length) {
        for (const rule of label.autoRules) {
          rows.push({ label, rule })
        }
      }
      if (label.children?.length) {
        traverse(label.children)
      }
    }
  }

  traverse(labels)
  return rows
}

export function AutoRulesDataTable({
  data,
  searchable = false,
  maxHeight = 400,
  fullscreen = false,
  fullscreenTitle,
  className,
}: AutoRulesDataTableProps) {
  const t = useT()
  const [isFullscreen, setIsFullscreen] = useState(false)

  // Create columns with i18n support
  const columns = useMemo(() => createColumns(t), [t])

  // Flatten label tree into auto-rule rows
  const rows = useMemo(() => collectAutoRules(data), [data])

  // Fullscreen button (shown on hover)
  const fullscreenButton = fullscreen ? (
    <button
      onClick={() => setIsFullscreen(true)}
      className={cn(
        'p-1 rounded-[6px] transition-all',
        'opacity-0 group-hover:opacity-100',
        'bg-background/80 backdrop-blur-sm shadow-minimal',
        'text-muted-foreground/50 hover:text-foreground',
        'focus:outline-none focus-visible:ring-1 focus-visible:ring-ring focus-visible:opacity-100'
      )}
      title={t('全屏查看')}
    >
      <Maximize2 className="w-3.5 h-3.5" />
    </button>
  ) : undefined

  return (
    <>
      <Info_DataTable
        columns={columns}
        data={rows}
        searchable={searchable ? { placeholder: t('搜索规则...') } : false}
        maxHeight={maxHeight}
        emptyContent={t('未配置自动应用规则')}
        floatingAction={fullscreenButton}
        className={cn(fullscreen && 'group', className)}
      />

      {/* Fullscreen overlay */}
      {fullscreen && (
        <DataTableOverlay
          isOpen={isFullscreen}
          onClose={() => setIsFullscreen(false)}
          title={fullscreenTitle || t('自动应用规则')}
          subtitle={`${rows.length} ${t('条规则')}`}
        >
          <Info_DataTable
            columns={columns}
            data={rows}
            searchable={searchable ? { placeholder: t('搜索规则...') } : false}
            emptyContent={t('未配置自动应用规则')}
          />
        </DataTableOverlay>
      )}
    </>
  )
}
