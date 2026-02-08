/**
 * InteractiveRequest - Structured input wrapper for interactive UI requests
 *
 * Renders the interactive UI component based on the request type
 * and handles response submission back to the agent.
 *
 * Supports both single elements and trees with multiple children.
 * When tree has children, renders them as a unified form.
 */

import * as React from 'react'
import { cn } from '@/lib/utils'
import type {
  InteractiveRequest as InteractiveRequestType,
  InteractiveResponse,
  InteractiveUIElement,
} from '@sprouty-ai/shared/interactive-ui'
import { InteractiveFormRenderer } from '../../../interactive-ui'

interface InteractiveRequestProps {
  request: InteractiveRequestType
  onResponse: (response: InteractiveResponse) => void
  /** When true, removes container styling - used when wrapped by InputContainer */
  unstyled?: boolean
}

/**
 * Flatten tree structure to array of elements.
 * If tree has children, use children. Otherwise use tree itself as single element.
 */
function flattenTree(tree: InteractiveUIElement): InteractiveUIElement[] {
  // If tree has children, use children as the elements
  if (tree.children && tree.children.length > 0) {
    return tree.children
  }
  // Otherwise treat the tree itself as a single element
  return [tree]
}

export function InteractiveRequest({ request, onResponse, unstyled = false }: InteractiveRequestProps) {
  const [completed, setCompleted] = React.useState(false)
  const [submittedResponses, setSubmittedResponses] = React.useState<Record<string, unknown> | null>(null)

  // Flatten tree to elements array
  const elements = React.useMemo(() => flattenTree(request.tree), [request.tree])

  const handleRespond = React.useCallback(
    (response: InteractiveResponse) => {
      setCompleted(true)
      // For form responses, data is Record<string, unknown>
      setSubmittedResponses(response.data as Record<string, unknown>)
      onResponse(response)
    },
    [onResponse]
  )

  return (
    <div
      className={cn(
        'w-full',
        !unstyled && 'p-4 bg-background border rounded-xl shadow-sm'
      )}
    >
      <InteractiveFormRenderer
        elements={elements}
        requestId={request.requestId}
        onRespond={handleRespond}
        completed={completed}
        previousResponses={submittedResponses ?? undefined}
        prompt={request.prompt}
      />
    </div>
  )
}
