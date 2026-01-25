/**
 * Confirm - Yes/No confirmation component
 *
 * Renders a confirmation dialog with title, message, and two action buttons.
 */

import * as React from 'react'
import { AlertCircle, Info, HelpCircle, CheckCircle } from 'lucide-react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import type { ConfirmProps } from '@creator-flow/shared/interactive-ui'

interface ConfirmComponentProps {
  props: ConfirmProps
  onConfirm: (confirmed: boolean) => void
}

export function Confirm({ props, onConfirm }: ConfirmComponentProps) {
  const {
    title,
    message,
    confirmLabel = 'Yes',
    cancelLabel = 'No',
    variant = 'default',
  } = props

  // Variant styling
  const variants = {
    default: {
      icon: HelpCircle,
      iconClass: 'text-primary',
      confirmClass: '',
    },
    info: {
      icon: Info,
      iconClass: 'text-blue-500',
      confirmClass: '',
    },
    warning: {
      icon: AlertCircle,
      iconClass: 'text-yellow-500',
      confirmClass: 'bg-yellow-500 hover:bg-yellow-600',
    },
    danger: {
      icon: AlertCircle,
      iconClass: 'text-destructive',
      confirmClass: 'bg-destructive hover:bg-destructive/90',
    },
    success: {
      icon: CheckCircle,
      iconClass: 'text-green-500',
      confirmClass: 'bg-green-500 hover:bg-green-600',
    },
  }

  const variantConfig = variants[variant] || variants.default
  const IconComponent = variantConfig.icon

  return (
    <div className="rounded-lg border bg-card p-4 space-y-4 max-w-md">
      {/* Header with icon and title */}
      <div className="flex items-start gap-3">
        <div className={cn('shrink-0 mt-0.5', variantConfig.iconClass)}>
          <IconComponent className="w-5 h-5" />
        </div>
        <div className="flex-1 min-w-0">
          <h4 className="text-sm font-semibold text-foreground">{title}</h4>
          {message && (
            <p className="text-sm text-muted-foreground mt-1">{message}</p>
          )}
        </div>
      </div>

      {/* Action buttons */}
      <div className="flex justify-end gap-2 pt-2">
        <Button
          variant="outline"
          size="sm"
          onClick={() => onConfirm(false)}
        >
          {cancelLabel}
        </Button>
        <Button
          size="sm"
          className={cn(variantConfig.confirmClass)}
          onClick={() => onConfirm(true)}
        >
          {confirmLabel}
        </Button>
      </div>
    </div>
  )
}
