/**
 * SingleChoice - Radio button selection component
 *
 * Renders a list of options as radio buttons.
 * Only one option can be selected at a time.
 * User must confirm selection with a button.
 */

import * as React from 'react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import type { SingleChoiceProps, ChoiceOption } from '@creator-flow/shared/interactive-ui'

interface SingleChoiceComponentProps {
  props: SingleChoiceProps
  onSelect: (value: string) => void
  selectedValue?: string
  /** Whether this component is in completed state (read-only) */
  disabled?: boolean
}

export function SingleChoice({ props, onSelect, selectedValue, disabled }: SingleChoiceComponentProps) {
  const { label, options, layout = 'vertical', required = true } = props
  const [selected, setSelected] = React.useState<string | undefined>(selectedValue ?? props.value)
  const [error, setError] = React.useState<string | null>(null)

  const handleSelect = (optionId: string) => {
    if (disabled) return
    setSelected(optionId)
    setError(null) // Clear error when user makes a selection
  }

  const handleConfirm = () => {
    // Validate selection
    if (required && !selected) {
      setError('Please select an option')
      return
    }
    if (selected) {
      onSelect(selected)
    }
  }

  const canConfirm = !required || !!selected

  return (
    <div className="space-y-3">
      {/* Label/Question */}
      <div className="flex items-baseline gap-1">
        <span className="text-sm font-medium text-foreground">{label}</span>
        {required && <span className="text-destructive">*</span>}
      </div>

      {/* Options */}
      <div
        className={cn(
          'gap-2',
          layout === 'horizontal' ? 'flex flex-wrap' : 'flex flex-col'
        )}
      >
        {options.map((option) => (
          <OptionItem
            key={option.id}
            option={option}
            selected={selected === option.id}
            onClick={() => handleSelect(option.id)}
            layout={layout}
            disabled={disabled || option.disabled}
          />
        ))}
      </div>

      {/* Error message */}
      {error && (
        <p className="text-xs text-destructive">{error}</p>
      )}

      {/* Confirm button */}
      {!disabled && (
        <div className="flex justify-end pt-2">
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

interface OptionItemProps {
  option: ChoiceOption
  selected: boolean
  onClick: () => void
  layout: 'vertical' | 'horizontal'
  disabled?: boolean
}

function OptionItem({ option, selected, onClick, layout, disabled }: OptionItemProps) {
  return (
    <button
      type="button"
      disabled={disabled}
      onClick={onClick}
      className={cn(
        'flex items-start gap-3 p-3 rounded-lg border transition-all text-left',
        !disabled && 'hover:bg-muted/50',
        selected
          ? 'border-primary bg-primary/5 ring-1 ring-primary/20'
          : 'border-border',
        !disabled && !selected && 'hover:border-foreground/20',
        disabled && 'opacity-60 cursor-not-allowed',
        layout === 'horizontal' && 'flex-1 min-w-[140px]'
      )}
    >
      {/* Radio indicator */}
      <div
        className={cn(
          'w-4 h-4 mt-0.5 rounded-full border-2 flex items-center justify-center shrink-0 transition-colors',
          selected ? 'border-primary' : 'border-muted-foreground/40'
        )}
      >
        {selected && <div className="w-2 h-2 rounded-full bg-primary" />}
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
