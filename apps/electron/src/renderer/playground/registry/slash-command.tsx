import * as React from 'react'
import type { ComponentEntry } from './types'
import { SquareSlash } from 'lucide-react'
import { Button } from '@/components/ui/button'
import {
  SlashCommandMenu,
  InlineSlashCommand,
  useInlineSlashCommand,
  DEFAULT_SLASH_COMMANDS,
  type SlashCommandId,
} from '@/components/ui/slash-command-menu'

// ============================================================================
// SlashCommandDemo - Full interactive demo
// ============================================================================

function SlashCommandDemo() {
  const inputRef = React.useRef<{ getBoundingClientRect: () => DOMRect; value: string; selectionStart: number }>(null)
  const textareaRef = React.useRef<HTMLTextAreaElement>(null)
  const [inputValue, setInputValue] = React.useState('')
  const [activeCommands, setActiveCommands] = React.useState<SlashCommandId[]>([])
  const [buttonMenuOpen, setButtonMenuOpen] = React.useState(false)

  // Sync inputRef with textarea values for the hook
  React.useEffect(() => {
    if (textareaRef.current) {
      (inputRef as React.MutableRefObject<{ getBoundingClientRect: () => DOMRect; value: string; selectionStart: number }>).current = {
        getBoundingClientRect: () => textareaRef.current!.getBoundingClientRect(),
        get value() { return textareaRef.current?.value ?? '' },
        get selectionStart() { return textareaRef.current?.selectionStart ?? 0 },
      }
    }
  }, [])

  // Inline slash command hook
  const inlineSlash = useInlineSlashCommand({
    inputRef,
  })

  const handleInputChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    const value = e.target.value
    const cursorPosition = e.target.selectionStart
    setInputValue(value)
    inlineSlash.handleInputChange(value, cursorPosition)
  }

  const handleInlineCommandSelect = (commandName: string, hasArgs: boolean) => {
    const { value: newValue } = inlineSlash.handleSelectCommand(commandName, hasArgs)
    setInputValue(newValue)
    setTimeout(() => textareaRef.current?.focus(), 0)
  }

  const handleButtonSelect = (commandId: SlashCommandId) => {
    setActiveCommands(prev =>
      prev.includes(commandId)
        ? prev.filter(id => id !== commandId)
        : [...prev, commandId]
    )
    setButtonMenuOpen(false)
    textareaRef.current?.focus()
  }

  return (
    <div className="flex flex-col h-full">
      {/* Description */}
      <div className="shrink-0 p-4 border-b border-border/50">
        <h2 className="text-sm font-medium text-foreground/80 mb-2">
          Slash Command Menu Demo
        </h2>
        <p className="text-xs text-muted-foreground">
          Type <code className="px-1 py-0.5 bg-muted rounded">/</code> to trigger inline autocomplete, or click the button to open the menu.
        </p>
      </div>

      {/* Main Content Area */}
      <div className="flex-1 flex gap-4 p-4">
        {/* Left: Button Menu with Filter */}
        <div className="flex-1">
          <div className="text-xs font-medium text-muted-foreground mb-2">
            Button Menu (with filter input)
          </div>
          <div className="relative">
            <Button
              variant="outline"
              size="sm"
              className="gap-2"
              onClick={() => setButtonMenuOpen(!buttonMenuOpen)}
            >
              <SquareSlash className="h-4 w-4" />
              Commands
            </Button>
            {buttonMenuOpen && (
              <div className="absolute top-full left-0 mt-2 z-10">
                <SlashCommandMenu
                  commands={DEFAULT_SLASH_COMMANDS}
                  activeCommands={activeCommands}
                  onSelect={handleButtonSelect}
                  showFilter={true}
                  className="w-[240px]"
                />
              </div>
            )}
          </div>
        </div>

        {/* Right: Static Menu (no filter) */}
        <div className="flex-1">
          <div className="text-xs font-medium text-muted-foreground mb-2">
            Static Menu (no filter)
          </div>
          <SlashCommandMenu
            commands={DEFAULT_SLASH_COMMANDS}
            activeCommands={activeCommands}
            onSelect={handleButtonSelect}
            className="w-full"
          />
        </div>
      </div>

      {/* Input Area with Inline Autocomplete */}
      <div className="shrink-0 p-4 border-t border-border/50">
        <div className="text-xs font-medium text-muted-foreground mb-2">
          Inline Autocomplete (type / in the textarea)
        </div>
        <div className="relative rounded-lg border bg-background shadow-sm">
          <textarea
            ref={textareaRef}
            value={inputValue}
            onChange={handleInputChange}
            placeholder="Type / to see commands..."
            className="w-full min-h-[80px] px-4 py-3 text-sm bg-transparent outline-none resize-none"
            rows={3}
          />
          <div className="flex items-center gap-2 px-3 py-2 border-t border-border/50">
            <Button
              variant="ghost"
              size="icon"
              className="h-7 w-7"
              onClick={() => setButtonMenuOpen(!buttonMenuOpen)}
            >
              <SquareSlash className="h-4 w-4" />
            </Button>
            <span className="text-xs text-muted-foreground flex-1">
              {inlineSlash.isOpen ? `Filtering: "${inlineSlash.filter}"` : 'Type / to trigger'}
            </span>
          </div>
        </div>

        {/* Inline Slash Command Menu */}
        <InlineSlashCommand
          open={inlineSlash.isOpen}
          onOpenChange={(open) => {
            if (!open) inlineSlash.close()
          }}
          sections={inlineSlash.sections}
          onSelectCommand={handleInlineCommandSelect}
          filter={inlineSlash.filter}
          position={inlineSlash.position}
        />
      </div>
    </div>
  )
}

// ============================================================================
// Basic Menu Preview (for props-based customization)
// ============================================================================

interface SlashCommandMenuPlaygroundProps {
  showFilter?: boolean
  filterPlaceholder?: string
}

function SlashCommandMenuPlayground({
  showFilter = true,
  filterPlaceholder = 'Search commands...',
}: SlashCommandMenuPlaygroundProps) {
  const [activeCommands, setActiveCommands] = React.useState<SlashCommandId[]>([])

  const handleSelect = (commandId: SlashCommandId) => {
    setActiveCommands(prev =>
      prev.includes(commandId)
        ? prev.filter(id => id !== commandId)
        : [...prev, commandId]
    )
  }

  return (
    <SlashCommandMenu
      commands={DEFAULT_SLASH_COMMANDS}
      activeCommands={activeCommands}
      onSelect={handleSelect}
      showFilter={showFilter}
      filterPlaceholder={filterPlaceholder}
      className="w-[280px]"
    />
  )
}

// ============================================================================
// Component Registry Entries
// ============================================================================

export const slashCommandComponents: ComponentEntry[] = [
  {
    id: 'slash-command-demo',
    name: 'Slash Command Demo',
    category: 'Chat Inputs',
    description: 'Interactive demo showing both button-triggered and inline slash command menus',
    component: SlashCommandDemo,
    layout: 'full',
    props: [],
    variants: [],
    mockData: () => ({}),
  },
  {
    id: 'slash-command-menu',
    name: 'SlashCommandMenu',
    category: 'Chat Inputs',
    description: 'Command palette for slash commands using cmdk',
    component: SlashCommandMenuPlayground,
    props: [
      {
        name: 'showFilter',
        description: 'Show filter input above the menu',
        control: { type: 'boolean' },
        defaultValue: true,
      },
      {
        name: 'filterPlaceholder',
        description: 'Placeholder text for filter input',
        control: { type: 'string', placeholder: 'Search...' },
        defaultValue: 'Search commands...',
      },
    ],
    variants: [
      { name: 'With Filter', props: { showFilter: true } },
      { name: 'Without Filter', props: { showFilter: false } },
    ],
    mockData: () => ({}),
  },
]
