/**
 * DataTable - Table display component
 *
 * Renders tabular data with optional row selection.
 */

import * as React from 'react'
import { Check } from 'lucide-react'
import { cn } from '@/lib/utils'
import type { DataTableProps, TableColumn, TableRow } from '@sprouty-ai/shared/interactive-ui'

interface DataTableComponentProps {
  props: DataTableProps
  onSelect?: (selectedIds: string[]) => void
  selectedIds?: string[]
}

export function DataTable({ props, onSelect, selectedIds }: DataTableComponentProps) {
  const { columns, rows, selectable, title } = props
  const [selected, setSelected] = React.useState<string[]>(selectedIds ?? [])

  const handleToggleRow = (rowId: string) => {
    if (!selectable) return

    const isSelected = selected.includes(rowId)
    const newSelected = isSelected
      ? selected.filter((id) => id !== rowId)
      : [...selected, rowId]

    setSelected(newSelected)
    onSelect?.(newSelected)
  }

  const handleToggleAll = () => {
    if (!selectable) return

    const allSelected = selected.length === rows.length
    const newSelected = allSelected ? [] : rows.map((r) => r.id)

    setSelected(newSelected)
    onSelect?.(newSelected)
  }

  const allSelected = selected.length === rows.length && rows.length > 0
  const someSelected = selected.length > 0 && selected.length < rows.length

  return (
    <div className="rounded-lg border bg-card overflow-hidden">
      {/* Title */}
      {title && (
        <div className="px-4 py-3 border-b bg-muted/30">
          <h4 className="text-sm font-semibold text-foreground">{title}</h4>
        </div>
      )}

      {/* Table */}
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead className="bg-muted/50">
            <tr>
              {/* Selection column */}
              {selectable && (
                <th className="w-10 px-3 py-2">
                  <button
                    type="button"
                    onClick={handleToggleAll}
                    className={cn(
                      'w-4 h-4 rounded border-2 flex items-center justify-center transition-colors',
                      allSelected || someSelected
                        ? 'border-primary bg-primary'
                        : 'border-muted-foreground/40'
                    )}
                  >
                    {allSelected && (
                      <Check className="w-3 h-3 text-primary-foreground" strokeWidth={3} />
                    )}
                    {someSelected && !allSelected && (
                      <div className="w-2 h-0.5 bg-primary-foreground rounded" />
                    )}
                  </button>
                </th>
              )}

              {/* Column headers */}
              {columns.map((column) => (
                <th
                  key={column.id}
                  className={cn(
                    'px-3 py-2 text-left font-medium text-muted-foreground',
                    column.width && `w-[${column.width}]`,
                    column.align === 'center' && 'text-center',
                    column.align === 'right' && 'text-right'
                  )}
                >
                  {column.label}
                </th>
              ))}
            </tr>
          </thead>

          <tbody className="divide-y divide-border">
            {rows.map((row) => (
              <TableRowComponent
                key={row.id}
                row={row}
                columns={columns}
                selectable={selectable}
                isSelected={selected.includes(row.id)}
                onToggle={() => handleToggleRow(row.id)}
              />
            ))}
          </tbody>
        </table>
      </div>

      {/* Empty state */}
      {rows.length === 0 && (
        <div className="px-4 py-8 text-center text-sm text-muted-foreground">
          No data available
        </div>
      )}

      {/* Selection summary */}
      {selectable && selected.length > 0 && (
        <div className="px-4 py-2 border-t bg-muted/30 text-xs text-muted-foreground">
          {selected.length} of {rows.length} row(s) selected
        </div>
      )}
    </div>
  )
}

interface TableRowComponentProps {
  row: TableRow
  columns: TableColumn[]
  selectable?: boolean
  isSelected: boolean
  onToggle: () => void
}

function TableRowComponent({
  row,
  columns,
  selectable,
  isSelected,
  onToggle,
}: TableRowComponentProps) {
  return (
    <tr
      className={cn(
        'hover:bg-muted/30 transition-colors',
        isSelected && 'bg-primary/5',
        selectable && 'cursor-pointer'
      )}
      onClick={selectable ? onToggle : undefined}
    >
      {/* Selection cell */}
      {selectable && (
        <td className="w-10 px-3 py-2">
          <div
            className={cn(
              'w-4 h-4 rounded border-2 flex items-center justify-center transition-colors',
              isSelected ? 'border-primary bg-primary' : 'border-muted-foreground/40'
            )}
          >
            {isSelected && (
              <Check className="w-3 h-3 text-primary-foreground" strokeWidth={3} />
            )}
          </div>
        </td>
      )}

      {/* Data cells */}
      {columns.map((column) => {
        const value = row.data[column.id]

        return (
          <td
            key={column.id}
            className={cn(
              'px-3 py-2 text-foreground',
              column.align === 'center' && 'text-center',
              column.align === 'right' && 'text-right'
            )}
          >
            <CellRenderer value={value} type={column.type} />
          </td>
        )
      })}
    </tr>
  )
}

interface CellRendererProps {
  value: unknown
  type?: 'text' | 'number' | 'date' | 'boolean' | 'badge'
}

function CellRenderer({ value, type }: CellRendererProps) {
  if (value === null || value === undefined) {
    return <span className="text-muted-foreground">—</span>
  }

  switch (type) {
    case 'boolean':
      return (
        <span
          className={cn(
            'inline-flex items-center justify-center w-4 h-4 rounded-full',
            value ? 'bg-green-100 text-green-600' : 'bg-red-100 text-red-600'
          )}
        >
          {value ? '✓' : '✗'}
        </span>
      )

    case 'badge':
      return (
        <span className="inline-flex px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">
          {String(value)}
        </span>
      )

    case 'date':
      if (value instanceof Date) {
        return <span>{value.toLocaleDateString()}</span>
      }
      return <span>{String(value)}</span>

    case 'number':
      return <span className="font-mono">{String(value)}</span>

    default:
      return <span>{String(value)}</span>
  }
}
