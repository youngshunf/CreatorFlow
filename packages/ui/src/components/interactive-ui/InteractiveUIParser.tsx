/**
 * InteractiveUIParser - Parses <interactive-ui> blocks from AI output
 *
 * Extracts JSON UI definitions from markdown/text content and renders them
 * as interactive components, following json-render's approach.
 */

import * as React from 'react'
import type { InteractiveUIElement, InteractiveResponse } from '@creator-flow/shared/interactive-ui'
import { InteractiveFormRenderer } from './InteractiveUIRenderer'

// Regex to match <interactive-ui>...</interactive-ui> blocks
const INTERACTIVE_UI_REGEX = /<interactive-ui>\s*([\s\S]*?)\s*<\/interactive-ui>/g

interface ParsedSegment {
  type: 'text' | 'interactive'
  content: string
  elements?: InteractiveUIElement[]
  parseError?: string
}

/**
 * Parse content into segments of text and interactive UI
 */
function parseContent(content: string): ParsedSegment[] {
  const segments: ParsedSegment[] = []
  let lastIndex = 0

  INTERACTIVE_UI_REGEX.lastIndex = 0

  let match
  while ((match = INTERACTIVE_UI_REGEX.exec(content)) !== null) {
    // Add text before this match
    if (match.index > lastIndex) {
      const textContent = content.slice(lastIndex, match.index).trim()
      if (textContent) {
        segments.push({ type: 'text', content: textContent })
      }
    }

    // Parse the JSON inside the tag
    const jsonContent = (match[1] ?? '').trim()
    try {
      const parsed = JSON.parse(jsonContent)

      if (parsed.elements && Array.isArray(parsed.elements)) {
        const elements: InteractiveUIElement[] = parsed.elements.map((el: unknown, index: number) => {
          const element = el as Record<string, unknown>
          return {
            key: (element.key as string) || `element-${index}`,
            type: element.type as InteractiveUIElement['type'],
            props: (element.props as Record<string, unknown>) || {},
          }
        })

        segments.push({
          type: 'interactive',
          content: jsonContent,
          elements,
        })
      } else {
        segments.push({
          type: 'interactive',
          content: jsonContent,
          parseError: 'Invalid structure: missing "elements" array',
        })
      }
    } catch (err) {
      segments.push({
        type: 'interactive',
        content: jsonContent,
        parseError: err instanceof Error ? err.message : 'JSON parse error',
      })
    }

    lastIndex = match.index + match[0].length
  }

  // Add remaining text
  if (lastIndex < content.length) {
    const textContent = content.slice(lastIndex).trim()
    if (textContent) {
      segments.push({ type: 'text', content: textContent })
    }
  }

  return segments
}

export interface InteractiveUIParserProps {
  /** The raw content that may contain <interactive-ui> blocks */
  content: string
  /** Called when user submits interactive response (includes elements for label lookup) */
  onInteractiveResponse?: (response: InteractiveResponse, elements?: InteractiveUIElement[]) => void
  /** Render function for text segments (e.g., markdown renderer) */
  renderText?: (text: string, index: number) => React.ReactNode
  /** Unique ID for tracking */
  messageId?: string
  /** Additional className for wrapper */
  className?: string
}

/**
 * Renders content with interactive UI blocks
 *
 * Parses <interactive-ui> tags and renders them as interactive components.
 * Text segments are rendered via the renderText prop (allowing custom markdown rendering).
 */
export function InteractiveUIParser({
  content,
  onInteractiveResponse,
  renderText,
  messageId,
  className,
}: InteractiveUIParserProps) {
  const [submittedBlocks, setSubmittedBlocks] = React.useState<Set<number>>(new Set())
  const [responses, setResponses] = React.useState<Map<number, Record<string, unknown>>>(new Map())

  const segments = React.useMemo(() => parseContent(content), [content])
  const hasInteractiveBlocks = segments.some((s) => s.type === 'interactive' && s.elements)

  // If no interactive blocks, render all as text
  if (!hasInteractiveBlocks) {
    return (
      <div className={className}>
        {renderText ? renderText(content, 0) : <p>{content}</p>}
      </div>
    )
  }

  const handleResponse = (blockIndex: number, response: InteractiveResponse, elements?: InteractiveUIElement[]) => {
    setSubmittedBlocks((prev) => new Set(prev).add(blockIndex))
    setResponses((prev) => new Map(prev).set(blockIndex, response.data as Record<string, unknown>))
    onInteractiveResponse?.(response, elements)
  }

  return (
    <div className={className}>
      {segments.map((segment, index) => {
        if (segment.type === 'text') {
          return (
            <div key={`text-${index}`}>
              {renderText ? renderText(segment.content, index) : <p>{segment.content}</p>}
            </div>
          )
        }

        if (segment.parseError) {
          return (
            <div
              key={`error-${index}`}
              className="my-4 p-4 rounded-lg border border-destructive/30 bg-destructive/5"
            >
              <div className="text-xs text-destructive/70 mb-1">Interactive UI Parse Error</div>
              <div className="text-sm text-destructive">{segment.parseError}</div>
              <pre className="mt-2 text-xs text-muted-foreground overflow-auto">
                {segment.content.slice(0, 200)}
                {segment.content.length > 200 ? '...' : ''}
              </pre>
            </div>
          )
        }

        if (segment.elements) {
          const isCompleted = submittedBlocks.has(index)
          const previousResponses = responses.get(index)

          return (
            <div key={`interactive-${index}`} className="my-4">
              <InteractiveFormRenderer
                elements={segment.elements}
                requestId={`inline-${messageId || 'msg'}-${index}`}
                onRespond={(response) => handleResponse(index, response, segment.elements)}
                completed={isCompleted}
                previousResponses={previousResponses}
              />
            </div>
          )
        }

        return null
      })}
    </div>
  )
}

/**
 * Check if content contains interactive UI blocks
 */
export function hasInteractiveUI(content: string): boolean {
  INTERACTIVE_UI_REGEX.lastIndex = 0
  return INTERACTIVE_UI_REGEX.test(content)
}
