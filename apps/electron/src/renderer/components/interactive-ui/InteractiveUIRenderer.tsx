/**
 * InteractiveUIRenderer - Main entry point for rendering interactive UI elements
 *
 * Supports two modes:
 * 1. Single element mode - renders one interactive component
 * 2. Multi-element form mode - combines multiple questions into a single form
 *
 * When multiple elements are provided, they are rendered as form fields
 * with a single submit button at the end.
 */

import * as React from 'react'
import { Button } from '@/components/ui/button'
import { cn } from '@/lib/utils'
import type {
  InteractiveUIElement,
  SingleChoiceProps,
  MultiChoiceProps,
  ConfirmProps,
  FormProps,
  CardProps,
  ListProps,
  DataTableProps,
  ButtonGroupProps,
  InteractiveResponse,
} from '@creator-flow/shared/interactive-ui'

// Local alias for convenience
type InteractiveElement = InteractiveUIElement

import { SingleChoice } from './components/SingleChoice'
import { MultiChoice } from './components/MultiChoice'
import { Confirm } from './components/Confirm'
import { Form } from './components/Form'
import { Card } from './components/Card'
import { List } from './components/List'
import { DataTable } from './components/DataTable'
import { ButtonGroup } from './components/ButtonGroup'
import { requiresResponse, type InteractiveResponseValue } from './registry'

// ============================================
// Multi-Element Form Renderer (Primary)
// ============================================

interface InteractiveFormRendererProps {
  /** Multiple interactive elements to render as a form */
  elements: InteractiveElement[]
  /** Request ID for tracking the response */
  requestId: string
  /** Called when user submits all responses */
  onRespond: (response: InteractiveResponse) => void
  /** Whether the interaction is already completed */
  completed?: boolean
  /** Previous response values (keyed by element key) */
  previousResponses?: Record<string, unknown>
  /** Optional prompt/title to display above the form */
  prompt?: string
}

/**
 * Renders multiple interactive elements as a unified form.
 * User fills in all fields and submits with a single button.
 */
export function InteractiveFormRenderer({
  elements,
  requestId,
  onRespond,
  completed = false,
  previousResponses,
  prompt,
}: InteractiveFormRendererProps) {
  // Track responses for each element by key
  const [responses, setResponses] = React.useState<Record<string, unknown>>(
    previousResponses ?? {}
  )
  // Track validation errors
  const [errors, setErrors] = React.useState<Record<string, string>>({})

  // Update response for a specific element
  const updateResponse = React.useCallback((key: string, value: unknown) => {
    setResponses((prev) => ({ ...prev, [key]: value }))
    // Clear error when user provides a value
    setErrors((prev) => {
      const next = { ...prev }
      delete next[key]
      return next
    })
  }, [])

  // Validate all required fields
  const validate = React.useCallback((): boolean => {
    const newErrors: Record<string, string> = {}

    for (const element of elements) {
      const value = responses[element.key]
      const props = element.props as Record<string, unknown>

      // Check required fields based on element type
      switch (element.type) {
        case 'single-choice': {
          const required = props.required !== false // Default to required
          if (required && !value) {
            newErrors[element.key] = 'Please select an option'
          }
          break
        }
        case 'multi-choice': {
          const min = props.min as number | undefined
          const values = value as string[] | undefined
          if (min !== undefined && (!values || values.length < min)) {
            newErrors[element.key] = `Please select at least ${min} option${min > 1 ? 's' : ''}`
          }
          break
        }
        case 'confirm': {
          // Confirm is always required - must choose yes or no
          if (value === undefined) {
            newErrors[element.key] = 'Please make a selection'
          }
          break
        }
        // Form has its own validation
        // Card, List, DataTable, ButtonGroup are optional by nature
      }
    }

    setErrors(newErrors)
    return Object.keys(newErrors).length === 0
  }, [elements, responses])

  // Submit all responses
  const handleSubmit = React.useCallback(() => {
    if (!validate()) return

    const response: InteractiveResponse = {
      requestId,
      type: 'submit',
      data: responses, // Form data keyed by element key
    }
    onRespond(response)
  }, [requestId, responses, onRespond, validate])

  // Check if form can be submitted (has all required values)
  const hasRequiredValues = React.useMemo(() => {
    for (const element of elements) {
      const value = responses[element.key]
      const props = element.props as Record<string, unknown>

      switch (element.type) {
        case 'single-choice': {
          const required = props.required !== false
          if (required && !value) return false
          break
        }
        case 'multi-choice': {
          const min = props.min as number | undefined
          const values = value as string[] | undefined
          if (min !== undefined && (!values || values.length < min)) return false
          break
        }
        case 'confirm': {
          if (value === undefined) return false
          break
        }
      }
    }
    return true
  }, [elements, responses])

  return (
    <div className="interactive-form-container space-y-4">
      {/* Optional prompt/title */}
      {prompt && (
        <p className="text-sm font-medium text-foreground">{prompt}</p>
      )}

      {/* Form fields */}
      <div className={cn('space-y-6', completed && 'opacity-75 pointer-events-none')}>
        {elements.map((element) => (
          <FormFieldWrapper
            key={element.key}
            element={element}
            value={responses[element.key]}
            onChange={(value) => updateResponse(element.key, value)}
            error={errors[element.key]}
            disabled={completed}
          />
        ))}
      </div>

      {/* Submit button */}
      {!completed && (
        <div className="flex justify-end pt-4 mt-2 border-t">
          <Button
            onClick={handleSubmit}
            disabled={!hasRequiredValues}
            className={cn(
              'px-6 py-2.5 font-semibold shadow-sm',
              hasRequiredValues && 'hover:shadow-md active:scale-[0.98]'
            )}
          >
            提交
          </Button>
        </div>
      )}

      {/* Completed indicator */}
      {completed && (
        <p className="text-xs text-muted-foreground italic">Response submitted</p>
      )}
    </div>
  )
}

// ============================================
// Form Field Wrapper (renders individual elements)
// ============================================

interface FormFieldWrapperProps {
  element: InteractiveElement
  value: unknown
  onChange: (value: unknown) => void
  error?: string
  disabled?: boolean
}

function FormFieldWrapper({ element, value, onChange, error, disabled }: FormFieldWrapperProps) {
  const renderField = () => {
    switch (element.type) {
      case 'single-choice':
        return (
          <SingleChoiceField
            props={element.props as SingleChoiceProps}
            value={value as string | undefined}
            onChange={onChange}
            disabled={disabled}
          />
        )

      case 'multi-choice':
        return (
          <MultiChoiceField
            props={element.props as MultiChoiceProps}
            values={value as string[] | undefined}
            onChange={onChange}
            disabled={disabled}
          />
        )

      case 'confirm':
        return (
          <ConfirmField
            props={element.props as ConfirmProps}
            value={value as boolean | undefined}
            onChange={onChange}
            disabled={disabled}
          />
        )

      case 'form':
        // Nested form - render inline
        return (
          <Form
            props={element.props as FormProps}
            onSubmit={(values) => onChange(values)}
          />
        )

      case 'card':
        return (
          <Card
            props={element.props as CardProps}
            onAction={(actionId) => onChange(actionId)}
          />
        )

      case 'list':
        return (
          <List
            props={element.props as ListProps}
            onSelect={(selectedIds) => onChange(selectedIds)}
            selectedIds={value as string[] | undefined}
          />
        )

      case 'data-table':
        return (
          <DataTable
            props={element.props as DataTableProps}
            onSelect={(selectedIds) => onChange(selectedIds)}
            selectedIds={value as string[] | undefined}
          />
        )

      case 'button-group':
        return (
          <ButtonGroup
            props={element.props as ButtonGroupProps}
            onAction={(actionId) => onChange(actionId)}
          />
        )

      default:
        return (
          <div className="p-4 rounded-lg border border-destructive/50 bg-destructive/5 text-sm text-destructive">
            Unknown component type: {element.type}
          </div>
        )
    }
  }

  return (
    <div className="form-field space-y-1">
      {renderField()}
      {error && <p className="text-xs text-destructive">{error}</p>}
    </div>
  )
}

// ============================================
// Field Components (no submit button, controlled)
// ============================================

/** SingleChoice as controlled field (no submit button) */
function SingleChoiceField({
  props,
  value,
  onChange,
  disabled,
}: {
  props: SingleChoiceProps
  value?: string
  onChange: (value: string) => void
  disabled?: boolean
}) {
  const { label, options, layout = 'vertical', required = true } = props

  return (
    <div className="space-y-3">
      <div className="flex items-baseline gap-1">
        <span className="text-sm font-medium text-foreground">{label}</span>
        {required && <span className="text-destructive">*</span>}
      </div>
      <div
        className={cn(
          'gap-2',
          layout === 'horizontal' ? 'flex flex-wrap' : 'flex flex-col'
        )}
      >
        {options.map((option) => (
          <button
            key={option.id}
            type="button"
            disabled={disabled || option.disabled}
            onClick={() => onChange(option.id)}
            className={cn(
              'flex items-start gap-3 p-3 rounded-lg border-2 transition-all text-left',
              !disabled && 'hover:bg-muted/50',
              value === option.id
                ? 'border-accent bg-accent/10 ring-2 ring-accent/40 shadow-sm'
                : 'border-border',
              !disabled && value !== option.id && 'hover:border-foreground/20',
              (disabled || option.disabled) && 'opacity-60 cursor-not-allowed',
              layout === 'horizontal' && 'flex-1 min-w-[140px]'
            )}
          >
            <div
              className={cn(
                'w-4 h-4 mt-0.5 rounded-full border-2 flex items-center justify-center shrink-0 transition-all',
                value === option.id
                  ? 'border-accent bg-accent/20'
                  : 'border-muted-foreground/40 bg-transparent'
              )}
            >
              {value === option.id && <div className="w-2 h-2 rounded-full bg-accent" />}
            </div>
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
        ))}
      </div>
    </div>
  )
}

/** MultiChoice as controlled field (no submit button) */
import { Check } from 'lucide-react'

function MultiChoiceField({
  props,
  values,
  onChange,
  disabled,
}: {
  props: MultiChoiceProps
  values?: string[]
  onChange: (values: string[]) => void
  disabled?: boolean
}) {
  const { label, options, layout = 'vertical', min, max } = props
  const selected = values ?? []

  const handleToggle = (optionId: string) => {
    if (disabled) return
    const isSelected = selected.includes(optionId)
    if (isSelected) {
      onChange(selected.filter((id) => id !== optionId))
    } else {
      if (max !== undefined && selected.length >= max) return
      onChange([...selected, optionId])
    }
  }

  const constraintHint = React.useMemo(() => {
    if (min !== undefined && max !== undefined) {
      if (min === max) return `Select exactly ${min}`
      return `Select ${min}-${max} options`
    }
    if (min !== undefined) return `Select at least ${min}`
    if (max !== undefined) return `Select up to ${max}`
    return null
  }, [min, max])

  return (
    <div className="space-y-3">
      <div className="flex items-baseline gap-2">
        <span className="text-sm font-medium text-foreground">{label}</span>
        {min !== undefined && <span className="text-destructive">*</span>}
        {constraintHint && (
          <span className="text-xs text-muted-foreground">({constraintHint})</span>
        )}
      </div>
      <div
        className={cn(
          'gap-2',
          layout === 'horizontal' ? 'flex flex-wrap' : 'flex flex-col'
        )}
      >
        {options.map((option) => {
          const checked = selected.includes(option.id)
          const isDisabled = disabled || option.disabled ||
            (!checked && max !== undefined && selected.length >= max)

          return (
            <button
              key={option.id}
              type="button"
              disabled={isDisabled}
              onClick={() => handleToggle(option.id)}
              className={cn(
                'flex items-start gap-3 p-3 rounded-lg border-2 transition-all text-left',
                !isDisabled && 'hover:bg-muted/50',
                checked
                  ? 'border-accent bg-accent/10 ring-2 ring-accent/40 shadow-sm'
                  : 'border-border',
                !isDisabled && !checked && 'hover:border-foreground/20',
                isDisabled && 'opacity-60 cursor-not-allowed',
                layout === 'horizontal' && 'flex-1 min-w-[140px]'
              )}
            >
              <div
                className={cn(
                  'w-4 h-4 mt-0.5 rounded border-2 flex items-center justify-center shrink-0 transition-all',
                  checked ? 'border-accent bg-accent' : 'border-muted-foreground/40 bg-transparent'
                )}
              >
                {checked && <Check className="w-3 h-3 text-accent-foreground" strokeWidth={3} />}
              </div>
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
        })}
      </div>
      {selected.length > 0 && (
        <p className="text-xs text-muted-foreground">
          {selected.length} selected{max !== undefined ? ` / ${max}` : ''}
        </p>
      )}
    </div>
  )
}

/** Confirm as controlled field (inline buttons) */
function ConfirmField({
  props,
  value,
  onChange,
  disabled,
}: {
  props: ConfirmProps
  value?: boolean
  onChange: (value: boolean) => void
  disabled?: boolean
}) {
  const { title, message, confirmLabel = 'Yes', cancelLabel = 'No' } = props

  return (
    <div className="space-y-3">
      <div>
        <span className="text-sm font-medium text-foreground">{title}</span>
        <span className="text-destructive">*</span>
      </div>
      {message && <p className="text-sm text-muted-foreground">{message}</p>}
      <div className="flex gap-2">
        <Button
          type="button"
          variant={value === true ? 'default' : 'outline'}
          size="sm"
          onClick={() => onChange(true)}
          disabled={disabled}
        >
          {confirmLabel}
        </Button>
        <Button
          type="button"
          variant={value === false ? 'default' : 'outline'}
          size="sm"
          onClick={() => onChange(false)}
          disabled={disabled}
        >
          {cancelLabel}
        </Button>
      </div>
    </div>
  )
}

// ============================================
// Single Element Renderer (Legacy/Simple)
// ============================================

interface InteractiveUIRendererProps {
  /** The interactive element to render */
  element: InteractiveElement
  /** Request ID for tracking the response */
  requestId: string
  /** Called when user submits a response */
  onRespond: (response: InteractiveResponse) => void
  /** Whether the interaction is already completed */
  completed?: boolean
  /** Previous response value (for displaying completed state) */
  previousResponse?: unknown
  /** Optional prompt to display above the component */
  prompt?: string
}

/**
 * Renders a single interactive element.
 * For multiple elements, use InteractiveFormRenderer instead.
 */
export function InteractiveUIRenderer({
  element,
  requestId,
  onRespond,
  completed = false,
  previousResponse,
  prompt,
}: InteractiveUIRendererProps) {
  // For single element, wrap in array and use form renderer
  return (
    <InteractiveFormRenderer
      elements={[element]}
      requestId={requestId}
      onRespond={onRespond}
      completed={completed}
      previousResponses={previousResponse ? { [element.key]: previousResponse } : undefined}
      prompt={prompt}
    />
  )
}

/**
 * Check if an element type requires user interaction
 */
export { requiresResponse }
