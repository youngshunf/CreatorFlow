/**
 * Interactive UI Component Registry
 *
 * Maps component types to their React implementations.
 */

import type { InteractiveComponentType, InteractiveUIElement } from '@creator-flow/shared/interactive-ui'

// Re-export component type
export type { InteractiveComponentType }

// Response value types
export type InteractiveResponseValue =
  | string // SingleChoice selection
  | string[] // MultiChoice selections, List selections, DataTable selections
  | boolean // Confirm result
  | Record<string, unknown> // Form submission

// Component callback types
export interface InteractiveCallbacks {
  onSingleChoice?: (value: string) => void
  onMultiChoice?: (values: string[]) => void
  onConfirm?: (confirmed: boolean) => void
  onFormSubmit?: (values: Record<string, unknown>) => void
  onCardAction?: (actionId: string) => void
  onButtonAction?: (actionId: string) => void
  onListSelect?: (selectedIds: string[]) => void
  onTableSelect?: (selectedIds: string[]) => void
}

// Determine if a component type requires user response
export function requiresResponse(type: InteractiveComponentType): boolean {
  const responsiveTypes: InteractiveComponentType[] = [
    'single-choice',
    'multi-choice',
    'confirm',
    'form',
    'button-group',
  ]
  return responsiveTypes.includes(type)
}

// Determine if a component type is display-only
export function isDisplayOnly(type: InteractiveComponentType): boolean {
  const displayTypes: InteractiveComponentType[] = ['card', 'list', 'data-table']
  return displayTypes.includes(type)
}

// Extract response from element based on type
export function extractResponseValue(
  element: InteractiveUIElement,
  response: unknown
): InteractiveResponseValue | null {
  switch (element.type) {
    case 'single-choice':
      return typeof response === 'string' ? response : null

    case 'multi-choice':
      return Array.isArray(response) ? response : null

    case 'confirm':
      return typeof response === 'boolean' ? response : null

    case 'form':
      return response && typeof response === 'object' && !Array.isArray(response)
        ? (response as Record<string, unknown>)
        : null

    case 'button-group':
      return typeof response === 'string' ? response : null

    case 'list':
    case 'data-table':
      // These can optionally return selected IDs
      return Array.isArray(response) ? response : null

    case 'card':
      // Card action returns action ID
      return typeof response === 'string' ? response : null

    default:
      return null
  }
}

// Default props for missing optional fields
export const defaultProps = {
  singleChoice: {
    layout: 'vertical' as const,
  },
  multiChoice: {
    layout: 'vertical' as const,
  },
  confirm: {
    confirmLabel: 'Yes',
    cancelLabel: 'No',
    variant: 'default' as const,
  },
  form: {
    submitLabel: 'Submit',
  },
  buttonGroup: {
    layout: 'horizontal' as const,
    align: 'start' as const,
  },
  list: {
    ordered: false,
    selectable: false,
    compact: false,
  },
  dataTable: {
    selectable: false,
  },
}
