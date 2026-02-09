/**
 * Card - Information card component
 *
 * Displays structured information with optional image and action buttons.
 */

import * as React from 'react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import type { CardProps, CardAction } from '@sprouty-ai/shared/interactive-ui'

interface CardComponentProps {
  props: CardProps
  onAction?: (actionId: string) => void
}

export function Card({ props, onAction }: CardComponentProps) {
  const { title, description, image, footer, metadata, actions } = props

  return (
    <div className="rounded-lg border bg-card overflow-hidden max-w-sm">
      {/* Image */}
      {image && (
        <div className="relative aspect-video bg-muted">
          <img
            src={image}
            alt={title}
            className="w-full h-full object-cover"
            loading="lazy"
          />
        </div>
      )}

      {/* Content */}
      <div className="p-4 space-y-3">
        {/* Title */}
        <h4 className="text-sm font-semibold text-foreground line-clamp-2">{title}</h4>

        {/* Description */}
        {description && (
          <p className="text-sm text-muted-foreground line-clamp-3">{description}</p>
        )}

        {/* Metadata */}
        {metadata && metadata.length > 0 && (
          <div className="flex flex-wrap gap-2">
            {metadata.map((item, index) => (
              <span
                key={index}
                className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full bg-muted text-xs text-muted-foreground"
              >
                {item.icon && <span>{item.icon}</span>}
                {item.label && <span className="font-medium">{item.label}:</span>}
                <span>{item.value}</span>
              </span>
            ))}
          </div>
        )}

        {/* Footer text */}
        {footer && (
          <p className="text-xs text-muted-foreground pt-2 border-t">{footer}</p>
        )}

        {/* Actions */}
        {actions && actions.length > 0 && (
          <div className="flex flex-wrap gap-2 pt-2">
            {actions.map((action) => (
              <CardActionButton
                key={action.id}
                action={action}
                onClick={() => onAction?.(action.id)}
              />
            ))}
          </div>
        )}
      </div>
    </div>
  )
}

interface CardActionButtonProps {
  action: CardAction
  onClick: () => void
}

function CardActionButton({ action, onClick }: CardActionButtonProps) {
  const variantMap: Record<string, 'default' | 'outline' | 'ghost' | 'destructive'> = {
    primary: 'default',
    secondary: 'outline',
    ghost: 'ghost',
    danger: 'destructive',
  }

  return (
    <Button
      variant={variantMap[action.variant || 'secondary'] || 'outline'}
      size="sm"
      onClick={onClick}
      disabled={action.disabled}
      className="text-xs"
    >
      {action.icon && <span className="mr-1">{action.icon}</span>}
      {action.label}
    </Button>
  )
}
