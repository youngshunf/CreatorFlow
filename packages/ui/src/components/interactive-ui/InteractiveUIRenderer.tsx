/**
 * InteractiveUIRenderer - Renders interactive UI elements from AI output
 *
 * Supports form components (single-choice, multi-choice, confirm) that
 * collect user input, and display components (card, list) for information.
 *
 * Form components render with a Submit button. User responses are returned
 * via the onRespond callback.
 */

import * as React from 'react'
import { Check } from 'lucide-react'
import type {
  InteractiveUIElement,
  SingleChoiceProps,
  MultiChoiceProps,
  ConfirmProps,
  FormProps,
  InteractiveResponse,
} from '@creator-flow/shared/interactive-ui'
import { cn } from '../../lib/utils'

// ============================================
// Types
// ============================================

type InteractiveElement = InteractiveUIElement

/** Default translations (Chinese) */
const DEFAULT_TRANSLATIONS = {
  submit: '提交',
  submitted: '已提交',
  selected: '已选择',
  selectAtLeast: '至少选择',
  selectExactly: '请选择',
  selectUpTo: '最多选择',
  selectRange: '选择',
  options: '个选项',
  pleaseSelect: '请选择一个选项',
  pleaseSelectAtLeast: '请至少选择',
  pleaseConfirm: '请确认',
}

export type InteractiveUITranslations = typeof DEFAULT_TRANSLATIONS

interface InteractiveFormRendererProps {
  /** Elements to render as a form */
  elements: InteractiveElement[]
  /** Request ID for tracking */
  requestId: string
  /** Called when user submits */
  onRespond: (response: InteractiveResponse) => void
  /** Whether already completed */
  completed?: boolean
  /** Previous responses (for completed state display) */
  previousResponses?: Record<string, unknown>
  /** Optional prompt above the form */
  prompt?: string
  /** Optional translations for i18n */
  translations?: Partial<InteractiveUITranslations>
}

// ============================================
// Main Form Renderer
// ============================================

export function InteractiveFormRenderer({
  elements,
  requestId,
  onRespond,
  completed = false,
  previousResponses,
  prompt,
  translations: userTranslations,
}: InteractiveFormRendererProps) {
  const t = React.useMemo(() => ({ ...DEFAULT_TRANSLATIONS, ...userTranslations }), [userTranslations])
  const [responses, setResponses] = React.useState<Record<string, unknown>>(
    previousResponses ?? {}
  )
  const [errors, setErrors] = React.useState<Record<string, string>>({})

  const updateResponse = React.useCallback((key: string, value: unknown) => {
    setResponses((prev) => ({ ...prev, [key]: value }))
    setErrors((prev) => {
      const next = { ...prev }
      delete next[key]
      return next
    })
  }, [])

  const validate = React.useCallback((): boolean => {
    const newErrors: Record<string, string> = {}

    for (const element of elements) {
      const value = responses[element.key]
      const props = element.props as Record<string, unknown>

      switch (element.type) {
        case 'single-choice': {
          const required = props.required !== false
          if (required && !value) {
            newErrors[element.key] = t.pleaseSelect
          }
          break
        }
        case 'multi-choice': {
          const min = props.min as number | undefined
          const values = value as string[] | undefined
          if (min !== undefined && (!values || values.length < min)) {
            newErrors[element.key] = `${t.pleaseSelectAtLeast} ${min} ${t.options}`
          }
          break
        }
        case 'confirm': {
          if (value === undefined) {
            newErrors[element.key] = t.pleaseConfirm
          }
          break
        }
      }
    }

    setErrors(newErrors)
    return Object.keys(newErrors).length === 0
  }, [elements, responses, t])

  const handleSubmit = React.useCallback(() => {
    if (!validate()) return

    const response: InteractiveResponse = {
      requestId,
      type: 'submit',
      data: responses,
    }
    onRespond(response)
  }, [requestId, responses, onRespond, validate])

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

  // Check if there are any form elements that need submission
  const hasFormElements = elements.some((el) =>
    ['single-choice', 'multi-choice', 'confirm', 'form'].includes(el.type)
  )

  return (
    <div className="interactive-form-container space-y-4">
      {prompt && (
        <p className="text-sm font-medium text-foreground">{prompt}</p>
      )}

      <div className={cn('space-y-6', completed && 'opacity-75 pointer-events-none')}>
        {elements.map((element) => (
          <FormFieldWrapper
            key={element.key}
            element={element}
            value={responses[element.key]}
            onChange={(value) => updateResponse(element.key, value)}
            error={errors[element.key]}
            disabled={completed}
            translations={t}
          />
        ))}
      </div>

      {/* Submit button - only for form elements, not display elements */}
      {!completed && hasFormElements && (
        <div className="flex justify-end pt-4 mt-2 border-t border-border">
          <button
            onClick={handleSubmit}
            disabled={!hasRequiredValues}
            className={cn(
              'px-6 py-2.5 text-sm font-semibold rounded-lg transition-all shadow-sm',
              hasRequiredValues
                ? 'bg-accent text-background hover:bg-accent/90 hover:shadow-md active:scale-[0.98]'
                : 'bg-muted text-muted-foreground cursor-not-allowed opacity-60'
            )}
          >
            {t.submit}
          </button>
        </div>
      )}

      {completed && (
        <p className="text-xs text-muted-foreground italic">{t.submitted}</p>
      )}
    </div>
  )
}

// ============================================
// Form Field Wrapper
// ============================================

interface FormFieldWrapperProps {
  element: InteractiveElement
  value: unknown
  onChange: (value: unknown) => void
  error?: string
  disabled?: boolean
  translations: InteractiveUITranslations
}

function FormFieldWrapper({ element, value, onChange, error, disabled, translations }: FormFieldWrapperProps) {
  const renderField = () => {
    switch (element.type) {
      case 'single-choice':
        return (
          <SingleChoiceField
            props={element.props as SingleChoiceProps}
            value={value as string | undefined}
            onChange={onChange}
            disabled={disabled}
            translations={translations}
          />
        )

      case 'multi-choice':
        return (
          <MultiChoiceField
            props={element.props as MultiChoiceProps}
            values={value as string[] | undefined}
            onChange={onChange}
            disabled={disabled}
            translations={translations}
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
        // Nested form - render form fields
        const formProps = element.props as FormProps
        return (
          <div className="space-y-4">
            {formProps.title && (
              <p className="text-sm font-medium text-foreground">{formProps.title}</p>
            )}
            {formProps.description && (
              <p className="text-xs text-muted-foreground mb-3">{formProps.description}</p>
            )}
            <div className="space-y-3">
              {formProps.fields?.map((field) => (
                <FormField
                  key={field.id}
                  field={field}
                  value={(value as Record<string, unknown>)?.[field.id]}
                  onChange={(fieldValue) => {
                    const currentValue = (value as Record<string, unknown>) || {}
                    onChange({ ...currentValue, [field.id]: fieldValue })
                  }}
                  disabled={disabled}
                />
              ))}
            </div>
          </div>
        )

      // Display components - render as-is, no interaction needed
      case 'card':
      case 'list':
      case 'data-table':
      case 'button-group':
        return (
          <div className="text-sm text-muted-foreground">
            [Display component: {element.type}]
          </div>
        )

      default:
        return (
          <div className="p-3 rounded-lg border border-destructive/30 bg-destructive/5 text-sm text-destructive">
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
// Field Components
// ============================================

function SingleChoiceField({
  props,
  value,
  onChange,
  disabled,
  translations: t,
}: {
  props: SingleChoiceProps
  value?: string
  onChange: (value: string) => void
  disabled?: boolean
  translations: InteractiveUITranslations
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
      {options.map((option: { id: string; label: string; description?: string; icon?: string; disabled?: boolean }) => (
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
                : 'border-border bg-background',
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

function MultiChoiceField({
  props,
  values,
  onChange,
  disabled,
  translations: t,
}: {
  props: MultiChoiceProps
  values?: string[]
  onChange: (values: string[]) => void
  disabled?: boolean
  translations: InteractiveUITranslations
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
      if (min === max) return `${t.selectExactly} ${min} ${t.options}`
      return `${t.selectRange} ${min}-${max} ${t.options}`
    }
    if (min !== undefined) return `${t.selectAtLeast} ${min} ${t.options}`
    if (max !== undefined) return `${t.selectUpTo} ${max} ${t.options}`
    return null
  }, [min, max, t])

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
      {options.map((option: { id: string; label: string; description?: string; icon?: string; disabled?: boolean }) => {
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
                {checked && <Check className="w-3 h-3 text-background" strokeWidth={3} />}
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
          {t.selected} {selected.length}{max !== undefined ? ` / ${max}` : ''}
        </p>
      )}
    </div>
  )
}

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
        <button
          type="button"
          onClick={() => onChange(true)}
          disabled={disabled}
          className={cn(
            'px-3 py-1.5 text-sm font-medium rounded-md transition-colors',
            value === true
              ? 'bg-accent text-background'
              : 'border border-border hover:bg-muted',
            disabled && 'opacity-60 cursor-not-allowed'
          )}
        >
          {confirmLabel}
        </button>
        <button
          type="button"
          onClick={() => onChange(false)}
          disabled={disabled}
          className={cn(
            'px-3 py-1.5 text-sm font-medium rounded-md transition-colors',
            value === false
              ? 'bg-accent text-background'
              : 'border border-border hover:bg-muted',
            disabled && 'opacity-60 cursor-not-allowed'
          )}
        >
          {cancelLabel}
        </button>
      </div>
    </div>
  )
}

// ============================================
// Form Field Component
// ============================================

import type { FormField as FormFieldType } from '@creator-flow/shared/interactive-ui'

function FormField({
  field,
  value,
  onChange,
  disabled,
}: {
  field: FormFieldType
  value: unknown
  onChange: (value: unknown) => void
  disabled?: boolean
}) {
  const renderInput = () => {
    switch (field.type) {
      case 'text':
      case 'email':
      case 'password':
        return (
          <input
            type={field.type}
            id={field.id}
            value={(value as string) || ''}
            onChange={(e) => onChange(e.target.value)}
            placeholder={field.placeholder}
            disabled={disabled}
            className={cn(
              'w-full px-3 py-2 text-sm rounded-md border border-border bg-background',
              'focus:outline-none focus:ring-2 focus:ring-accent/40 focus:border-accent',
              'disabled:opacity-60 disabled:cursor-not-allowed',
              'transition-all'
            )}
          />
        )

      case 'number':
        return (
          <input
            type="number"
            id={field.id}
            value={(value as number) || ''}
            onChange={(e) => onChange(Number(e.target.value))}
            placeholder={field.placeholder}
            disabled={disabled}
            min={field.validation?.min}
            max={field.validation?.max}
            className={cn(
              'w-full px-3 py-2 text-sm rounded-md border border-border bg-background',
              'focus:outline-none focus:ring-2 focus:ring-accent/40 focus:border-accent',
              'disabled:opacity-60 disabled:cursor-not-allowed',
              'transition-all'
            )}
          />
        )

      case 'textarea':
        return (
          <textarea
            id={field.id}
            value={(value as string) || ''}
            onChange={(e) => onChange(e.target.value)}
            placeholder={field.placeholder}
            disabled={disabled}
            rows={3}
            className={cn(
              'w-full px-3 py-2 text-sm rounded-md border border-border bg-background',
              'focus:outline-none focus:ring-2 focus:ring-accent/40 focus:border-accent',
              'disabled:opacity-60 disabled:cursor-not-allowed',
              'transition-all resize-y'
            )}
          />
        )

      case 'select':
        return (
          <select
            id={field.id}
            value={(value as string) || ''}
            onChange={(e) => onChange(e.target.value)}
            disabled={disabled}
            className={cn(
              'w-full px-3 py-2 text-sm rounded-md border border-border bg-background',
              'focus:outline-none focus:ring-2 focus:ring-accent/40 focus:border-accent',
              'disabled:opacity-60 disabled:cursor-not-allowed',
              'transition-all'
            )}
          >
            <option value="">请选择...</option>
            {field.options?.map((option) => (
              <option key={option.id} value={option.id}>
                {option.label}
              </option>
            ))}
          </select>
        )

      case 'checkbox':
        return (
          <label className="flex items-center gap-2 cursor-pointer">
            <input
              type="checkbox"
              id={field.id}
              checked={(value as boolean) || false}
              onChange={(e) => onChange(e.target.checked)}
              disabled={disabled}
              className={cn(
                'w-4 h-4 rounded border-2 border-border',
                'focus:outline-none focus:ring-2 focus:ring-accent/40',
                'disabled:opacity-60 disabled:cursor-not-allowed',
                'transition-all'
              )}
            />
            <span className="text-sm text-foreground">{field.label}</span>
          </label>
        )

      case 'date':
        return (
          <input
            type="date"
            id={field.id}
            value={(value as string) || ''}
            onChange={(e) => onChange(e.target.value)}
            disabled={disabled}
            className={cn(
              'w-full px-3 py-2 text-sm rounded-md border border-border bg-background',
              'focus:outline-none focus:ring-2 focus:ring-accent/40 focus:border-accent',
              'disabled:opacity-60 disabled:cursor-not-allowed',
              'transition-all'
            )}
          />
        )

      default:
        return null
    }
  }

  // Checkbox has its own label
  if (field.type === 'checkbox') {
    return <div className="space-y-1">{renderInput()}</div>
  }

  return (
    <div className="space-y-1.5">
      <label htmlFor={field.id} className="flex items-baseline gap-1">
        <span className="text-sm font-medium text-foreground">{field.label}</span>
        {field.required && <span className="text-destructive">*</span>}
      </label>
      {renderInput()}
      {field.helpText && (
        <p className="text-xs text-muted-foreground">{field.helpText}</p>
      )}
    </div>
  )
}

// ============================================
// Legacy single-element renderer (wraps FormRenderer)
// ============================================

interface InteractiveUIRendererProps {
  element: InteractiveElement
  requestId: string
  onRespond: (response: InteractiveResponse) => void
  completed?: boolean
  previousResponse?: unknown
  prompt?: string
}

export function InteractiveUIRenderer({
  element,
  requestId,
  onRespond,
  completed = false,
  previousResponse,
  prompt,
}: InteractiveUIRendererProps) {
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
