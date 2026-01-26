/**
 * PermissionsDataTable
 *
 * Typed Data Table for displaying source permissions.
 * Features: searchable patterns, sortable columns, max-height scroll, fullscreen view.
 */

import * as React from 'react'
import { useState, useMemo } from 'react'
import type { ColumnDef } from '@tanstack/react-table'
import { Maximize2 } from 'lucide-react'
import { Info_DataTable, SortableHeader } from './Info_DataTable'
import { Info_Badge } from './Info_Badge'
import { Info_StatusBadge } from './Info_StatusBadge'
import { Tooltip, TooltipTrigger, TooltipContent } from '@creator-flow/ui'
import { DataTableOverlay } from '@creator-flow/ui'
import { cn } from '@/lib/utils'
import { useTheme } from '@/hooks/useTheme'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'

export type PermissionAccess = 'allowed' | 'blocked'
export type PermissionType = 'tool' | 'bash' | 'api' | 'mcp'

export interface PermissionRow {
  access: PermissionAccess
  type: PermissionType
  pattern: string
  comment?: string | null
}

interface PermissionsDataTableProps {
  data: PermissionRow[]
  /** Hide the type column (for MCP sources that only show pattern/comment) */
  hideTypeColumn?: boolean
  /** Show search input */
  searchable?: boolean
  /** Max height with scroll */
  maxHeight?: number
  /** Enable fullscreen button (shows Maximize2 icon on hover) */
  fullscreen?: boolean
  /** Title for the fullscreen overlay header */
  fullscreenTitle?: string
  className?: string
}

/**
 * PatternBadge - Clickable pattern badge with truncation and tooltip
 * - Dynamic width with max-width of 240px
 * - CSS truncation via text-ellipsis
 * - Tooltip shows full pattern on hover (only for patterns 30+ chars)
 * - Click to copy pattern to clipboard with toast notification
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
        <span className="block overflow-hidden whitespace-nowrap text-ellipsis max-w-[240px]">
          {pattern}
        </span>
      </Info_Badge>
    </button>
  )

  // Only show tooltip for longer patterns (30+ chars)
  if (pattern.length >= 30) {
    return (
      <Tooltip>
        <TooltipTrigger asChild>{badge}</TooltipTrigger>
        <TooltipContent className="font-mono max-w-md break-all">{pattern}</TooltipContent>
      </Tooltip>
    )
  }

  return badge
}

// Column definitions with sorting - created dynamically for i18n support
function createColumnsWithType(t: (key: string) => string): ColumnDef<PermissionRow>[] {
  return [
    {
      accessorKey: 'access',
      header: ({ column }) => <SortableHeader column={column} title={t('访问')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <Info_StatusBadge status={row.original.access} className="whitespace-nowrap" />
        </div>
      ),
      minSize: 80,
    },
    {
      accessorKey: 'type',
      header: ({ column }) => <SortableHeader column={column} title={t('类型')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <Info_Badge color="muted" className="capitalize whitespace-nowrap">
            {row.original.type}
          </Info_Badge>
        </div>
      ),
      minSize: 80,
    },
    {
      accessorKey: 'pattern',
      header: ({ column }) => <SortableHeader column={column} title={t('模式')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <PatternBadge pattern={row.original.pattern} t={t} />
        </div>
      ),
      minSize: 100,
    },
    {
      id: 'comment',
      accessorKey: 'comment',
      header: () => <span className="p-1.5 pl-2.5">{t('注释')}</span>,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5 min-w-0">
          <span className="truncate block">
            {row.original.comment ? t(row.original.comment) : '—'}
          </span>
        </div>
      ),
      meta: { fillWidth: true, truncate: true },
    },
  ]
}

function createColumnsWithoutType(t: (key: string) => string): ColumnDef<PermissionRow>[] {
  return [
    {
      accessorKey: 'access',
      header: ({ column }) => <SortableHeader column={column} title={t('访问')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <Info_StatusBadge status={row.original.access} className="whitespace-nowrap" />
        </div>
      ),
      minSize: 80,
    },
    {
      accessorKey: 'pattern',
      header: ({ column }) => <SortableHeader column={column} title={t('模式')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <PatternBadge pattern={row.original.pattern} t={t} />
        </div>
      ),
      minSize: 100,
    },
    {
      id: 'comment',
      accessorKey: 'comment',
      header: () => <span className="p-1.5 pl-2.5">{t('注释')}</span>,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5 min-w-0">
          <span className="truncate block">
            {row.original.comment ? t(row.original.comment) : '—'}
          </span>
        </div>
      ),
      meta: { fillWidth: true, truncate: true },
    },
  ]
}

export function PermissionsDataTable({
  data,
  hideTypeColumn = false,
  searchable = false,
  maxHeight = 400,
  fullscreen = false,
  fullscreenTitle,
  className,
}: PermissionsDataTableProps) {
  const t = useT()
  const [isFullscreen, setIsFullscreen] = useState(false)
  const { isDark } = useTheme()
  
  // Create columns with i18n support
  const columns = useMemo(
    () => hideTypeColumn ? createColumnsWithoutType(t) : createColumnsWithType(t),
    [hideTypeColumn, t]
  )

  // Fullscreen button for toolbar - shown on hover
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
        data={data}
        searchable={searchable ? { placeholder: t('搜索模式...') } : false}
        maxHeight={maxHeight}
        emptyContent={t('未配置权限')}
        floatingAction={fullscreenButton}
        className={cn(fullscreen && 'group', className)}
      />

      {/* Fullscreen overlay - renders the table without scroll constraints */}
      {fullscreen && (
        <DataTableOverlay
          isOpen={isFullscreen}
          onClose={() => setIsFullscreen(false)}
          title={fullscreenTitle || t('权限')}
          subtitle={`${data.length} ${t('条规则')}`}
          theme={isDark ? 'dark' : 'light'}
        >
          <Info_DataTable
            columns={columns}
            data={data}
            searchable={searchable ? { placeholder: t('搜索模式...') } : false}
            emptyContent={t('未配置权限')}
          />
        </DataTableOverlay>
      )}
    </>
  )
}
