/**
 * Interactive UI Components
 *
 * Provides components for rendering interactive UI from AI output.
 * AI can embed <interactive-ui> blocks in text, which are parsed and
 * rendered as interactive form elements.
 */

export {
  InteractiveFormRenderer,
  InteractiveUIRenderer,
} from './InteractiveUIRenderer'

export {
  InteractiveUIParser,
  hasInteractiveUI,
  type InteractiveUIParserProps,
} from './InteractiveUIParser'

// Re-export types from shared for convenience
export type {
  InteractiveUIElement,
  InteractiveResponse,
  SingleChoiceProps,
  MultiChoiceProps,
  ConfirmProps,
} from '@sprouty-ai/shared/interactive-ui'
