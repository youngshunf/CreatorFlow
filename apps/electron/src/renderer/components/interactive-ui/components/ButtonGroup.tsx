/**
 * ButtonGroup - Button row component
 *
 * Renders a group of action buttons in horizontal or vertical layout.
 */

import * as React from 'react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import type { ButtonGroupProps, ButtonProps } from '@sprouty-ai/shared/interactive-ui'

interface ButtonGroupComponentProps {
  props: ButtonGroupProps
  onAction: (actionId: string) => void
}

export function ButtonGroup({ props, onAction }: ButtonGroupComponentProps) {
  const { buttons, layout = 'horizontal', align = 'start' } = props

  const alignmentClass = {
    start: 'justify-start',
    center: 'justify-center',
    end: 'justify-end',
    between: 'justify-between',
    stretch: 'justify-stretch',
  }

  return (
    <div
      className={cn(
        'gap-2',
        layout === 'horizontal'
          ? `flex flex-wrap items-center ${alignmentClass[align]}`
          : 'flex flex-col items-stretch'
      )}
    >
      {buttons.map((button) => (
        <ActionButton key={button.action} button={button} onClick={onAction} />
      ))}
    </div>
  )
}

interface ActionButtonProps {
  button: ButtonProps
  onClick: (actionId: string) => void
}

function ActionButton({ button, onClick }: ActionButtonProps) {
  const { label, action, variant = 'primary', disabled, icon, size = 'default' } = button

  const variantMap: Record<string, 'default' | 'outline' | 'ghost' | 'destructive' | 'secondary'> = {
    primary: 'default',
    secondary: 'secondary',
    outline: 'outline',
    ghost: 'ghost',
    danger: 'destructive',
  }

  const sizeMap: Record<string, 'default' | 'sm' | 'lg' | 'icon'> = {
    small: 'sm',
    default: 'default',
    large: 'lg',
  }

  return (
    <Button
      variant={variantMap[variant] || 'default'}
      size={sizeMap[size] || 'default'}
      disabled={disabled}
      onClick={() => onClick(action)}
    >
      {icon && <span className="mr-1.5">{icon}</span>}
      {label}
    </Button>
  )
}

// Re-export for standalone button usage
export { ActionButton }
