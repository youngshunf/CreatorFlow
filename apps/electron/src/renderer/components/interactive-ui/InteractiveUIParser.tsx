/**
 * InteractiveUIParser - Parses <interactive-ui> blocks from AI output
 *
 * Extracts JSON UI definitions from markdown content and renders them
 * as interactive components, following json-render's approach.
 *
 * Usage: Wrap assistant message content with this component to enable
 * interactive UI rendering inline with markdown.
 */

import * as React from 'react'
import { Markdown, type RenderMode } from '@creator-flow/ui'
import { InteractiveFormRenderer } from './InteractiveUIRenderer'
import type { InteractiveUIElement, InteractiveResponse } from '@creator-flow/shared/interactive-ui'

// Regex to match <interactive-ui>...</interactive-ui> blocks
const INTERACTIVE_UI_REGEX = /<interactive-ui>\s*([\s\S]*?)\s*<\/interactive-ui>/g

interface ParsedSegment {
  type: 'text' | 'interactive'
  content: string
  // For interactive segments
  elements?: InteractiveUIElement[]
  parseError?: string
}

/**
 * Parse content into segments of text and interactive UI
 */
function parseContent(content: string): ParsedSegment[] {
  const segments: ParsedSegment[] = []
  let lastIndex = 0

  // Reset regex state
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
    const jsonContent = match[1].trim()
    try {
      const parsed = JSON.parse(jsonContent)

      // Validate structure
      if (parsed.elements && Array.isArray(parsed.elements)) {
        // Convert to InteractiveUIElement format
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
        // Invalid structure
        segments.push({
          type: 'interactive',
          content: jsonContent,
          parseError: 'Invalid structure: missing "elements" array',
        })
      }
    } catch (err) {
      // JSON parse error - show error state
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
  /** Markdown render mode */
  renderMode?: RenderMode
  /** Whether content is still streaming */
  isStreaming?: boolean
  /** Called when user submits interactive response */
  onInteractiveResponse?: (response: InteractiveResponse) => void
  /** URL click handler for markdown */
  onUrlClick?: (url: string) => void
  /** File click handler for markdown */
  onFileClick?: (path: string) => void
  /** Unique ID for markdown rendering */
  messageId?: string
  /** Additional className for wrapper */
  className?: string
  /** Session ID for sending response back to agent */
  sessionId?: string
}

/**
 * Renders content with interactive UI blocks
 *
 * Parses <interactive-ui> tags and renders them as interactive components,
 * while rendering regular content as markdown.
 */
export function InteractiveUIParser({
  content,
  renderMode = 'minimal',
  isStreaming = false,
  onInteractiveResponse,
  onUrlClick,
  onFileClick,
  messageId,
  className,
  sessionId,
}: InteractiveUIParserProps) {
  // Track which interactive blocks have been submitted
  const [submittedBlocks, setSubmittedBlocks] = React.useState<Set<number>>(new Set())
  const [responses, setResponses] = React.useState<Map<number, Record<string, unknown>>>(new Map())

  // Parse content into segments
  const segments = React.useMemo(() => parseContent(content), [content])

  // Check if there are any interactive blocks
  const hasInteractiveBlocks = segments.some((s) => s.type === 'interactive' && s.elements)

  // If no interactive blocks, just render markdown
  if (!hasInteractiveBlocks) {
    return (
      <Markdown
        mode={renderMode}
        onUrlClick={onUrlClick}
        onFileClick={onFileClick}
        id={messageId}
        className={className}
      >
        {content}
      </Markdown>
    )
  }

  // Handle interactive response
  const handleResponse = (blockIndex: number, response: InteractiveResponse) => {
    setSubmittedBlocks((prev) => new Set(prev).add(blockIndex))
    setResponses((prev) => new Map(prev).set(blockIndex, response.data as Record<string, unknown>))

    // Forward to parent callback if provided
    onInteractiveResponse?.(response)

    // Dispatch event to send response as user message to agent
    // Format: JSON representation of the form data
    if (sessionId && response.data) {
      const formattedResponse = JSON.stringify(response.data, null, 2)
      window.dispatchEvent(
        new CustomEvent('craft:interactive-response', {
          detail: {
            sessionId,
            text: formattedResponse,
          },
        })
      )
    }
  }

  return (
    <div className={className}>
      {segments.map((segment, index) => {
        if (segment.type === 'text') {
          return (
            <Markdown
              key={`text-${index}`}
              mode={renderMode}
              onUrlClick={onUrlClick}
              onFileClick={onFileClick}
              id={messageId ? `${messageId}-${index}` : undefined}
            >
              {segment.content}
            </Markdown>
          )
        }

        // Interactive segment
        if (segment.parseError) {
          // Show parse error
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
                onRespond={(response) => handleResponse(index, response)}
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
