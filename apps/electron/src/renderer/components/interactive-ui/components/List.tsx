/**
 * List - Item list component
 *
 * Renders a list of items with optional selection support.
 */

import * as React from 'react'
import { Check } from 'lucide-react'
import { cn } from '@/lib/utils'
import type { ListProps, ListItem } from '@sprouty-ai/shared/interactive-ui'

interface ListComponentProps {
  props: ListProps
  onSelect?: (selectedIds: string[]) => void
  selectedIds?: string[]
}

export function List({ props, onSelect, selectedIds }: ListComponentProps) {
  const { title, items, ordered = false, selectable = false, compact = false } = props
  const [selected, setSelected] = React.useState<string[]>(selectedIds ?? [])

  const handleToggle = (itemId: string) => {
    if (!selectable) return

    const isSelected = selected.includes(itemId)
    const newSelected = isSelected
      ? selected.filter((id) => id !== itemId)
      : [...selected, itemId]

    setSelected(newSelected)
    onSelect?.(newSelected)
  }

  const ListTag = ordered ? 'ol' : 'ul'

  return (
    <div className="rounded-lg border bg-card overflow-hidden">
      {/* Title */}
      {title && (
        <div className="px-4 py-3 border-b bg-muted/30">
          <h4 className="text-sm font-semibold text-foreground">{title}</h4>
        </div>
      )}

      {/* List items */}
      <ListTag className={cn('divide-y divide-border', ordered && 'list-decimal list-inside')}>
        {items.map((item, index) => (
          <ListItemComponent
            key={item.id}
            item={item}
            index={index}
            ordered={ordered}
            selectable={selectable}
            compact={compact}
            isSelected={selected.includes(item.id)}
            onToggle={() => handleToggle(item.id)}
          />
        ))}
      </ListTag>

      {/* Empty state */}
      {items.length === 0 && (
        <div className="px-4 py-8 text-center text-sm text-muted-foreground">
          No items
        </div>
      )}
    </div>
  )
}

interface ListItemComponentProps {
  item: ListItem
  index: number
  ordered: boolean
  selectable: boolean
  compact: boolean
  isSelected: boolean
  onToggle: () => void
}

function ListItemComponent({
  item,
  index,
  ordered,
  selectable,
  compact,
  isSelected,
  onToggle,
}: ListItemComponentProps) {
  return (
    <li
      className={cn(
        'flex items-start gap-3 transition-colors',
        compact ? 'px-3 py-2' : 'px-4 py-3',
        selectable && 'cursor-pointer hover:bg-muted/50',
        isSelected && 'bg-primary/5'
      )}
      onClick={selectable ? onToggle : undefined}
    >
      {/* Order number for ordered lists */}
      {ordered && (
        <span className="shrink-0 w-6 text-sm font-medium text-muted-foreground">
          {index + 1}.
        </span>
      )}

      {/* Selection checkbox */}
      {selectable && !ordered && (
        <div
          className={cn(
            'shrink-0 w-4 h-4 mt-0.5 rounded border-2 flex items-center justify-center transition-colors',
            isSelected ? 'border-primary bg-primary' : 'border-muted-foreground/40'
          )}
        >
          {isSelected && (
            <Check className="w-3 h-3 text-primary-foreground" strokeWidth={3} />
          )}
        </div>
      )}

      {/* Icon */}
      {item.icon && !ordered && (
        <span className="shrink-0 text-base mt-0.5">{item.icon}</span>
      )}

      {/* Content */}
      <div className="flex-1 min-w-0">
        <div className="flex items-center gap-2">
          <span className={cn('text-sm text-foreground', compact && 'text-xs')}>
            {item.label}
          </span>
          {item.badge && (
            <span className="inline-flex px-1.5 py-0.5 rounded-full bg-primary/10 text-primary text-[10px] font-medium">
              {item.badge}
            </span>
          )}
        </div>
        {item.description && (
          <p className={cn('text-muted-foreground mt-0.5', compact ? 'text-[10px]' : 'text-xs')}>
            {item.description}
          </p>
        )}
      </div>

      {/* Secondary content */}
      {item.secondary && (
        <span className="shrink-0 text-xs text-muted-foreground">
          {item.secondary}
        </span>
      )}
    </li>
  )
}
