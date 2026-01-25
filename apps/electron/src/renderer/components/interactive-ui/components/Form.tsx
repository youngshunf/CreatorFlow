/**
 * Form - Multi-field form component
 *
 * Renders a form with multiple field types:
 * text, number, email, password, textarea, select, checkbox
 */

import * as React from 'react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Label } from '@/components/ui/label'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import type { FormProps, FormField } from '@creator-flow/shared/interactive-ui'

interface FormComponentProps {
  props: FormProps
  onSubmit: (values: Record<string, unknown>) => void
  onCancel?: () => void
}

export function Form({ props, onSubmit, onCancel }: FormComponentProps) {
  const { title, description, fields, submitLabel = 'Submit', cancelLabel } = props

  // Initialize form values with defaults
  const [values, setValues] = React.useState<Record<string, unknown>>(() => {
    const initial: Record<string, unknown> = {}
    for (const field of fields) {
      if (field.defaultValue !== undefined) {
        initial[field.id] = field.defaultValue
      } else if (field.type === 'checkbox') {
        initial[field.id] = false
      } else {
        initial[field.id] = ''
      }
    }
    return initial
  })

  // Track validation errors
  const [errors, setErrors] = React.useState<Record<string, string>>({})

  const validateField = (field: FormField, value: unknown): string | null => {
    if (field.required) {
      if (value === '' || value === null || value === undefined) {
        return `${field.label} is required`
      }
    }
    if (field.validation) {
      if (field.validation.pattern) {
        const regex = new RegExp(field.validation.pattern)
        if (typeof value === 'string' && !regex.test(value)) {
          return field.validation.message || 'Invalid format'
        }
      }
      if (field.validation.minLength !== undefined && typeof value === 'string') {
        if (value.length < field.validation.minLength) {
          return `Minimum ${field.validation.minLength} characters`
        }
      }
      if (field.validation.maxLength !== undefined && typeof value === 'string') {
        if (value.length > field.validation.maxLength) {
          return `Maximum ${field.validation.maxLength} characters`
        }
      }
    }
    return null
  }

  const handleChange = (fieldId: string, value: unknown) => {
    setValues((prev) => ({ ...prev, [fieldId]: value }))
    // Clear error on change
    if (errors[fieldId]) {
      setErrors((prev) => {
        const next = { ...prev }
        delete next[fieldId]
        return next
      })
    }
  }

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault()

    // Validate all fields
    const newErrors: Record<string, string> = {}
    for (const field of fields) {
      const error = validateField(field, values[field.id])
      if (error) {
        newErrors[field.id] = error
      }
    }

    if (Object.keys(newErrors).length > 0) {
      setErrors(newErrors)
      return
    }

    onSubmit(values)
  }

  return (
    <form onSubmit={handleSubmit} className="rounded-lg border bg-card p-4 space-y-4 max-w-lg">
      {/* Header */}
      {(title || description) && (
        <div className="space-y-1">
          {title && <h4 className="text-sm font-semibold text-foreground">{title}</h4>}
          {description && (
            <p className="text-xs text-muted-foreground">{description}</p>
          )}
        </div>
      )}

      {/* Fields */}
      <div className="space-y-4">
        {fields.map((field) => (
          <FormFieldRenderer
            key={field.id}
            field={field}
            value={values[field.id]}
            error={errors[field.id]}
            onChange={(value) => handleChange(field.id, value)}
          />
        ))}
      </div>

      {/* Actions */}
      <div className="flex justify-end gap-2 pt-2">
        {cancelLabel && onCancel && (
          <Button type="button" variant="outline" size="sm" onClick={onCancel}>
            {cancelLabel}
          </Button>
        )}
        <Button type="submit" size="sm">
          {submitLabel}
        </Button>
      </div>
    </form>
  )
}

interface FormFieldRendererProps {
  field: FormField
  value: unknown
  error?: string
  onChange: (value: unknown) => void
}

function FormFieldRenderer({ field, value, error, onChange }: FormFieldRendererProps) {
  const id = `form-field-${field.id}`

  const renderInput = () => {
    switch (field.type) {
      case 'textarea':
        return (
          <Textarea
            id={id}
            placeholder={field.placeholder}
            value={String(value ?? '')}
            onChange={(e) => onChange(e.target.value)}
            disabled={field.disabled}
            rows={4}
            className={cn(error && 'border-destructive')}
          />
        )

      case 'select':
        return (
          <Select
            value={String(value ?? '')}
            onValueChange={onChange}
            disabled={field.disabled}
          >
            <SelectTrigger className={cn(error && 'border-destructive')}>
              <SelectValue placeholder={field.placeholder || 'Select...'} />
            </SelectTrigger>
            <SelectContent>
              {field.options?.map((option) => (
                <SelectItem key={option.id} value={option.id}>
                  {option.icon && <span className="mr-2">{option.icon}</span>}
                  {option.label}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        )

      case 'checkbox':
        return (
          <div className="flex items-center gap-2">
            <input
              type="checkbox"
              id={id}
              checked={Boolean(value)}
              onChange={(e) => onChange(e.target.checked)}
              disabled={field.disabled}
              className="h-4 w-4 rounded border-gray-300 text-primary focus:ring-primary"
            />
            <Label htmlFor={id} className="text-sm font-normal cursor-pointer">
              {field.label}
              {field.required && <span className="text-destructive ml-0.5">*</span>}
            </Label>
          </div>
        )

      case 'number':
        return (
          <Input
            id={id}
            type="number"
            placeholder={field.placeholder}
            value={value === '' || value === undefined ? '' : String(value)}
            onChange={(e) => onChange(e.target.value ? Number(e.target.value) : '')}
            disabled={field.disabled}
            className={cn(error && 'border-destructive')}
          />
        )

      default:
        return (
          <Input
            id={id}
            type={field.type || 'text'}
            placeholder={field.placeholder}
            value={String(value ?? '')}
            onChange={(e) => onChange(e.target.value)}
            disabled={field.disabled}
            className={cn(error && 'border-destructive')}
          />
        )
    }
  }

  // Checkbox has integrated label
  if (field.type === 'checkbox') {
    return (
      <div className="space-y-1">
        {renderInput()}
        {error && <p className="text-xs text-destructive">{error}</p>}
      </div>
    )
  }

  return (
    <div className="space-y-1.5">
      <Label htmlFor={id} className="text-sm">
        {field.label}
        {field.required && <span className="text-destructive ml-0.5">*</span>}
      </Label>
      {renderInput()}
      {field.hint && !error && (
        <p className="text-xs text-muted-foreground">{field.hint}</p>
      )}
      {error && <p className="text-xs text-destructive">{error}</p>}
    </div>
  )
}
