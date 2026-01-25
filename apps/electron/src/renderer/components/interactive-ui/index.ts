/**
 * Interactive UI Components
 *
 * Local implementation that wraps shared components with Markdown support.
 * The parser uses local Markdown rendering, while renderer uses shared UI.
 */

// Re-export from local implementation (supports Markdown rendering)
export { InteractiveUIParser, hasInteractiveUI, type InteractiveUIParserProps } from './InteractiveUIParser'

// Re-export from shared UI package
export { InteractiveUIRenderer, InteractiveFormRenderer } from '@creator-flow/ui'

// Re-export types from shared
export type { InteractiveUIElement, InteractiveResponse } from '@creator-flow/shared/interactive-ui'
