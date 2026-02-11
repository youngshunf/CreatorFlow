/**
 * Info_StatusBadge
 *
 * Status badge for permission states using Info_Badge.
 */

import * as React from 'react'
import { Info_Badge, type BadgeColor } from './Info_Badge'
import { useT } from '@/context/LocaleContext'

type PermissionStatus = 'allowed' | 'blocked' | 'requires-permission'

const statusConfig: Record<PermissionStatus, { labelKey: string; color: BadgeColor }> = {
  allowed: { labelKey: '已允许', color: 'success' },
  blocked: { labelKey: '已阻止', color: 'destructive' },
  'requires-permission': { labelKey: '需询问', color: 'warning' },
}

export interface Info_StatusBadgeProps
  extends Omit<React.HTMLAttributes<HTMLSpanElement>, 'children'> {
  /** Status type */
  status?: PermissionStatus | null
  /** Override the default label */
  label?: string
}

export function Info_StatusBadge({
  status,
  label,
  ...props
}: Info_StatusBadgeProps) {
  const t = useT()
  const key: PermissionStatus = status ?? 'allowed'
  const config = statusConfig[key]
  const displayLabel = label ?? t(config.labelKey)

  return (
    <Info_Badge {...props} color={config.color}>
      {displayLabel}
    </Info_Badge>
  )
}
