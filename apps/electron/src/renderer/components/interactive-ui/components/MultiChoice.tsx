/**
 * MultiChoice - Checkbox selection component
 *
 * Renders a list of checkboxes for multiple selection.
 * Supports min/max selection constraints.
 * User must confirm selection with a button.
 */

import * as React from 'react'
import { Check } from 'lucide-react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import type { MultiChoiceProps, ChoiceOption } from '@sprouty-ai/shared/interactive-ui'

interface MultiChoiceComponentProps {
  props: MultiChoiceProps
  onSelect: (values: string[]) => void
  selectedValues?: string[]
  /** Whether this component is in completed state (read-only) */
  disabled?: boolean
}

export function MultiChoice({ props, onSelect, selectedValues, disabled }: MultiChoiceComponentProps) {
  const { label, options, layout = 'vertical', min, max } = props
  const [selected, setSelected] = React.useState<string[]>(selectedValues ?? props.values ?? [])
  const [error, setError] = React.useState<string | null>(null)

  const handleToggle = (optionId: string) => {
    if (disabled) return

    const isSelected = selected.includes(optionId)
    let newSelected: string[]

    if (isSelected) {
      // Deselect - check min constraint (only warn, don't block)
      newSelected = selected.filter((id) => id !== optionId)
    } else {
      // Select - check max constraint
      if (max !== undefined && selected.length >= max) {
        return // Can't select more, at maximum
      }
      newSelected = [...selected, optionId]
    }

    setSelected(newSelected)
    setError(null) // Clear error when user makes a change
  }

  const handleConfirm = () => {
    // Validate min constraint
    if (min !== undefined && selected.length < min) {
      setError(`Please select at least ${min} option${min > 1 ? 's' : ''}`)
      return
    }
    onSelect(selected)
  }

  // Check if can confirm (meets minimum requirement)
  const canConfirm = min === undefined || selected.length >= min

  // Show constraint hint if applicable
  const constraintHint = React.useMemo(() => {
    if (min !== undefined && max !== undefined) {
      if (min === max) return `Select exactly ${min}`
      return `Select ${min}-${max} options`
    }
    if (min !== undefined) return `Select at least ${min}`
    if (max !== undefined) return `Select up to ${max}`
    return null
  }, [min, max])

  // Show current selection count
  const selectionStatus = React.useMemo(() => {
    if (max !== undefined) {
      return `${selected.length}/${max} selected`
    }
    if (selected.length > 0) {
      return `${selected.length} selected`
    }
    return null
  }, [selected.length, max])

  return (
    <div className="space-y-3">
      {/* Label/Question */}
      <div className="flex items-baseline gap-2">
        <span className="text-sm font-medium text-foreground">{label}</span>
        {min !== undefined && <span className="text-destructive">*</span>}
        {constraintHint && (
          <span className="text-xs text-muted-foreground">({constraintHint})</span>
        )}
      </div>

      {/* Options */}
      <div
        className={cn(
          'gap-2',
          layout === 'horizontal' ? 'flex flex-wrap' : 'flex flex-col'
        )}
      >
        {options.map((option) => (
          <CheckboxOptionItem
            key={option.id}
            option={option}
            checked={selected.includes(option.id)}
            onClick={() => handleToggle(option.id)}
            layout={layout}
            disabled={
              disabled ||
              option.disabled ||
              // Disable unchecked items if at max
              (!selected.includes(option.id) && max !== undefined && selected.length >= max)
            }
          />
        ))}
      </div>

      {/* Error message */}
      {error && (
        <p className="text-xs text-destructive">{error}</p>
      )}

      {/* Selection status and Confirm button */}
      {!disabled && (
        <div className="flex items-center justify-between pt-2">
          <span className="text-xs text-muted-foreground">
            {selectionStatus}
          </span>
          <Button
            size="sm"
            onClick={handleConfirm}
            disabled={!canConfirm}
          >
            Confirm
          </Button>
        </div>
      )}
    </div>
  )
}

interface CheckboxOptionItemProps {
  option: ChoiceOption
  checked: boolean
  onClick: () => void
  layout: 'vertical' | 'horizontal'
  disabled?: boolean
}

function CheckboxOptionItem({ option, checked, onClick, layout, disabled }: CheckboxOptionItemProps) {
  return (
    <button
      type="button"
      disabled={disabled}
      onClick={onClick}
      className={cn(
        'flex items-start gap-3 p-3 rounded-lg border transition-all text-left',
        !disabled && 'hover:bg-muted/50',
        checked
          ? 'border-primary bg-primary/5 ring-1 ring-primary/20'
          : 'border-border',
        !disabled && !checked && 'hover:border-foreground/20',
        disabled && 'opacity-60 cursor-not-allowed',
        layout === 'horizontal' && 'flex-1 min-w-[140px]'
      )}
    >
      {/* Checkbox indicator */}
      <div
        className={cn(
          'w-4 h-4 mt-0.5 rounded border-2 flex items-center justify-center shrink-0 transition-colors',
          checked ? 'border-primary bg-primary' : 'border-muted-foreground/40'
        )}
      >
        {checked && <Check className="w-3 h-3 text-primary-foreground" strokeWidth={3} />}
      </div>

      {/* Content */}
      <div className="flex-1 min-w-0">
        <div className="flex items-center gap-2">
          {option.icon && <span className="text-base">{option.icon}</span>}
          <span className="text-sm font-medium text-foreground">{option.label}</span>
        </div>
        {option.description && (
          <p className="text-xs text-muted-foreground mt-0.5">{option.description}</p>
        )}
      </div>
    </button>
  )
}
