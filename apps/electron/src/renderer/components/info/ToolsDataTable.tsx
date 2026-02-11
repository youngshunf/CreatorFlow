/**
 * ToolsDataTable
 *
 * Typed Data Table for displaying MCP tools.
 * Features: searchable tools, sortable columns, max-height scroll.
 */

import * as React from 'react'
import type { ColumnDef } from '@tanstack/react-table'
import { Info_DataTable, SortableHeader } from './Info_DataTable'
import { Info_Badge } from './Info_Badge'
import { Info_StatusBadge } from './Info_StatusBadge'
import { useT } from '@/context/LocaleContext'

export type ToolPermission = 'allowed' | 'requires-permission'

export interface ToolRow {
  name: string
  description: string
  permission: ToolPermission
}

interface ToolsDataTableProps {
  data: ToolRow[]
  /** Show loading spinner */
  loading?: boolean
  /** Show error message */
  error?: string
  /** Max height with scroll (default: 400) */
  maxHeight?: number
  className?: string
}

export function ToolsDataTable({
  data,
  loading,
  error,
  maxHeight = 400,
  className,
}: ToolsDataTableProps) {
  const t = useT()

  const columns = React.useMemo<ColumnDef<ToolRow>[]>(() => [
    {
      accessorKey: 'permission',
      header: ({ column }) => <SortableHeader column={column} title={t('访问')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <Info_StatusBadge status={row.original.permission} className="whitespace-nowrap" />
        </div>
      ),
      minSize: 80,
    },
    {
      accessorKey: 'name',
      header: ({ column }) => <SortableHeader column={column} title={t('工具')} />,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5">
          <Info_Badge color="muted" className="whitespace-nowrap">
            {row.original.name}
          </Info_Badge>
        </div>
      ),
      minSize: 100,
    },
    {
      id: 'description',
      accessorKey: 'description',
      header: () => <span className="p-1.5 pl-2.5">{t('描述')}</span>,
      cell: ({ row }) => (
        <div className="p-1.5 pl-2.5 min-w-0">
          <span className="truncate block">{row.original.description}</span>
        </div>
      ),
      meta: { fillWidth: true, truncate: true },
    },
  ], [t])

  return (
    <Info_DataTable
      columns={columns}
      data={data}
      loading={loading}
      error={error}
      maxHeight={maxHeight}
      emptyContent={t('没有可用的工具')}
      className={className}
    />
  )
}
