import * as React from 'react'
import { Command as CommandPrimitive } from 'cmdk'
import { Brain, Check } from 'lucide-react'
import { Icon_Folder } from '@creator-flow/ui'
import { cn } from '@/lib/utils'
import { PERMISSION_MODE_CONFIG, PERMISSION_MODE_ORDER, type PermissionMode } from '@creator-flow/shared/agent/modes'
import { SDK_COMMAND_TRANSLATIONS, getCommandDisplay } from '@creator-flow/shared/agent/slash-command-data'
import { sdkSlashCommandsAtom, commandTranslationsAtom } from '@/atoms/sdk-commands'
import { useAtomValue } from 'jotai'
import { useT } from '@/context/LocaleContext'

// ============================================================================
// Types
// ============================================================================

export type SlashCommandId = 'safe' | 'ask' | 'allow-all' | 'ultrathink'

/** Union type for all item types in the slash menu */
export type SlashItemType = 'command' | 'folder'

export interface SlashCommand {
  id: SlashCommandId
  label: string
  description: string
  icon: React.ReactNode
  shortcut?: string
  /** Optional color for the command (hex color string) */
  color?: string
}

/** Folder item for the slash menu */
export interface SlashFolderItem {
  id: string
  type: 'folder'
  label: string
  description: string
  path: string
}

/** Section with header for the inline slash menu */
export interface SlashSection {
  id: string
  label: string
  items: (SlashCommand | SlashFolderItem)[]
}

export interface CommandGroup {
  id: string
  commands: SlashCommand[]
}

// ============================================================================
// Permission Mode Icon Component
// ============================================================================

interface PermissionModeIconProps {
  mode: PermissionMode
  className?: string
}

function PermissionModeIcon({ mode, className }: PermissionModeIconProps) {
  const config = PERMISSION_MODE_CONFIG[mode]
  return (
    <svg
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth={2}
      strokeLinecap="round"
      strokeLinejoin="round"
      className={className}
    >
      <path d={config.svgPath} />
    </svg>
  )
}

// ============================================================================
// Default Commands
// ============================================================================

// Icon size constant
const MENU_ICON_SIZE = 'h-3.5 w-3.5'

/**
 * Permission mode translation mapping for UI display.
 * These provide localized labels and descriptions for permission modes.
 */
const PERMISSION_MODE_TRANSLATIONS: Record<PermissionMode, { label: string; description: string }> = {
  'safe': { label: '探索模式', description: '只读模式，禁止写入，不会提示' },
  'ask': { label: '询问模式', description: '编辑前先询问' },
  'allow-all': { label: '执行模式', description: '自动执行，不会提示' },
}

/**
 * Get localized permission mode commands.
 * Uses translation mapping for Chinese UI while keeping original IDs.
 */
export function getPermissionModeCommands(t: (text: string) => string): SlashCommand[] {
  return PERMISSION_MODE_ORDER.map(mode => {
    const config = PERMISSION_MODE_CONFIG[mode]
    const translation = PERMISSION_MODE_TRANSLATIONS[mode]
    return {
      id: mode as SlashCommandId,
      label: t(translation.label),
      description: t(translation.description),
      icon: <PermissionModeIcon mode={mode} className={MENU_ICON_SIZE} />,
    }
  })
}

/**
 * Get localized ultrathink command.
 */
export function getUltrathinkCommand(t: (text: string) => string): SlashCommand {
  return {
    id: 'ultrathink',
    label: t('深度思考'),
    description: t('用于复杂问题的扩展推理'),
    icon: <Brain className={MENU_ICON_SIZE} />,
  }
}

// Static versions for backward compatibility (English defaults)
const permissionModeCommands: SlashCommand[] = PERMISSION_MODE_ORDER.map(mode => {
  const config = PERMISSION_MODE_CONFIG[mode]
  return {
    id: mode as SlashCommandId,
    label: config.displayName,
    description: config.description,
    icon: <PermissionModeIcon mode={mode} className={MENU_ICON_SIZE} />,
  }
})

const ultrathinkCommand: SlashCommand = {
  id: 'ultrathink',
  label: 'Ultrathink',
  description: 'Extended reasoning for complex problems',
  icon: <Brain className={MENU_ICON_SIZE} />,
}

export const DEFAULT_SLASH_COMMANDS: SlashCommand[] = [
  ...permissionModeCommands,
  ultrathinkCommand,
]

export const DEFAULT_SLASH_COMMAND_GROUPS: CommandGroup[] = [
  { id: 'modes', commands: permissionModeCommands },
  { id: 'features', commands: [ultrathinkCommand] },
]

/**
 * Get localized command groups.
 * Use this instead of DEFAULT_SLASH_COMMAND_GROUPS for localized UI.
 */
export function getLocalizedCommandGroups(t: (text: string) => string): CommandGroup[] {
  return [
    { id: 'modes', commands: getPermissionModeCommands(t) },
    { id: 'features', commands: [getUltrathinkCommand(t)] },
  ]
}

// ============================================================================
// Shared Styles
// ============================================================================

const MENU_CONTAINER_STYLE = 'min-w-[200px] overflow-hidden rounded-[8px] bg-background text-foreground shadow-modal-small'
const MENU_LIST_STYLE = 'max-h-[260px] overflow-y-auto py-1'
const MENU_ITEM_STYLE = 'flex cursor-pointer select-none items-center gap-2 rounded-[6px] mx-1 px-2 py-1.5 text-[13px]'
const MENU_ITEM_SELECTED = 'bg-foreground/5'
const MENU_SECTION_HEADER = 'px-3 py-1.5 mb-0.5 text-[12px] font-medium text-muted-foreground border-b border-foreground/5'

// ============================================================================
// Shared: Filter utilities
// ============================================================================

function filterCommands(commands: SlashCommand[], filter: string): SlashCommand[] {
  if (!filter) return commands
  const lowerFilter = filter.toLowerCase()
  return commands.filter(
    cmd =>
      cmd.label.toLowerCase().includes(lowerFilter) ||
      cmd.id.toLowerCase().includes(lowerFilter)
  )
}

/** Check if an item is a folder */
function isFolder(item: SlashCommand | SlashFolderItem): item is SlashFolderItem {
  return 'type' in item && item.type === 'folder'
}

/** Filter sections by label/id, keeping sections grouped */
function filterSections(sections: SlashSection[], filter: string): SlashSection[] {
  if (!filter) return sections
  const lowerFilter = filter.toLowerCase()

  // Filter items within each section, keeping section structure
  return sections
    .map(section => ({
      ...section,
      items: section.items.filter(item =>
        item.label.toLowerCase().includes(lowerFilter) ||
        item.id.toLowerCase().includes(lowerFilter) ||
        item.description?.toLowerCase().includes(lowerFilter)
      ),
    }))
    .filter(section => section.items.length > 0)
}

/** Flatten sections into a single array of items */
function flattenSections(sections: SlashSection[]): (SlashCommand | SlashFolderItem)[] {
  return sections.flatMap(section => section.items)
}

// ============================================================================
// Shared: Command Item Content
// ============================================================================

function CommandItemContent({ command, isActive, showDescription = false }: { command: SlashCommand; isActive: boolean; showDescription?: boolean }) {
  return (
    <>
      <div className="shrink-0 text-muted-foreground">{command.icon}</div>
      <div className="flex-1 min-w-0">
        <div>{command.label}</div>
        {showDescription && command.description && (
          <div className="text-[11px] text-muted-foreground mt-0.5 leading-tight">{command.description}</div>
        )}
      </div>
      {isActive && (
        <div className="shrink-0 h-4 w-4 rounded-full bg-current flex items-center justify-center">
          <Check className="h-2.5 w-2.5 text-white dark:text-black" strokeWidth={3} />
        </div>
      )}
    </>
  )
}

// ============================================================================
// SlashCommandMenu Component (Button-triggered popup)
// ============================================================================

export interface SlashCommandMenuProps {
  /** Flat list of commands (use this OR commandGroups, not both) */
  commands?: SlashCommand[]
  /** Grouped commands with separators between groups */
  commandGroups?: CommandGroup[]
  activeCommands?: SlashCommandId[]
  onSelect: (commandId: SlashCommandId) => void
  showFilter?: boolean
  filterPlaceholder?: string
  /** Show description for each command item */
  showDescription?: boolean
  className?: string
}

export function SlashCommandMenu({
  commands,
  commandGroups,
  activeCommands = [],
  onSelect,
  showFilter = false,
  filterPlaceholder,
  showDescription = false,
  className,
}: SlashCommandMenuProps) {
  const t = useT()
  const [filter, setFilter] = React.useState('')
  const effectivePlaceholder = filterPlaceholder ?? t('搜索命令...')
  const inputRef = React.useRef<HTMLInputElement>(null)

  // If groups provided, filter within each group; otherwise use flat commands
  const filteredGroups = React.useMemo(() => {
    if (commandGroups) {
      return commandGroups.map(group => ({
        ...group,
        commands: filterCommands(group.commands, filter),
      })).filter(group => group.commands.length > 0)
    }
    return null
  }, [commandGroups, filter])

  const filteredCommands = React.useMemo(() => {
    if (commands && !commandGroups) {
      return filterCommands(commands, filter)
    }
    return null
  }, [commands, commandGroups, filter])

  // Get all commands for defaultValue calculation
  const allFilteredCommands = filteredGroups
    ? filteredGroups.flatMap(g => g.commands)
    : (filteredCommands ?? [])

  // Default to the first active command, or first command if none active
  const defaultValue = activeCommands[0] ?? allFilteredCommands[0]?.id

  React.useEffect(() => {
    if (showFilter && inputRef.current) {
      inputRef.current.focus()
    }
  }, [showFilter])

  if (allFilteredCommands.length === 0 && !showFilter) return null

  // Render a single command item
  const renderCommandItem = (cmd: SlashCommand) => {
    const isActive = activeCommands.includes(cmd.id)
    return (
      <CommandPrimitive.Item
        key={cmd.id}
        value={cmd.id}
        onSelect={() => onSelect(cmd.id)}
        data-tutorial={`permission-mode-${cmd.id}`}
        className={cn(
          MENU_ITEM_STYLE,
          'outline-none',
          'data-[selected=true]:bg-foreground/5',
          showDescription && 'py-2' // More padding when showing descriptions
        )}
      >
        <CommandItemContent command={cmd} isActive={isActive} showDescription={showDescription} />
      </CommandPrimitive.Item>
    )
  }

  return (
    <CommandPrimitive
      className={cn(MENU_CONTAINER_STYLE, className)}
      shouldFilter={false}
      defaultValue={defaultValue}
    >
      {showFilter && (
        <div className="border-b border-border/50 px-3 py-2">
          <CommandPrimitive.Input
            ref={inputRef}
            value={filter}
            onValueChange={setFilter}
            placeholder={effectivePlaceholder}
            className="w-full bg-transparent text-sm outline-none placeholder:text-muted-foreground"
          />
        </div>
      )}
      <CommandPrimitive.List className={MENU_LIST_STYLE}>
        {allFilteredCommands.length === 0 ? (
          <CommandPrimitive.Empty className="py-4 text-center text-sm text-muted-foreground">
            {t('未找到命令')}
          </CommandPrimitive.Empty>
        ) : filteredGroups ? (
          // Group-based rendering with smart separators
          filteredGroups.map((group, groupIndex) => (
            <React.Fragment key={group.id}>
              {group.commands.map(renderCommandItem)}
              {/* Separator: only show if there's another group after this one */}
              {groupIndex < filteredGroups.length - 1 && (
                <div className="h-px bg-border/50 my-1 mx-2" />
              )}
            </React.Fragment>
          ))
        ) : (
          // Flat list rendering
          filteredCommands?.map(renderCommandItem)
        )}
      </CommandPrimitive.List>
    </CommandPrimitive>
  )
}

// ============================================================================
// InlineSlashCommand - Autocomplete that follows cursor (shows SDK + plugin slash commands)
// ============================================================================

/** Item type for the inline slash command menu */
export interface SlashCommandItem {
  id: string
  label: string
  description: string
  hasArgs: boolean
}

/** Section for the inline slash command menu */
export interface InlineSlashSection {
  id: string
  label: string
  items: SlashCommandItem[]
}

export interface InlineSlashCommandProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  sections: InlineSlashSection[]
  onSelectCommand: (commandName: string, hasArgs: boolean) => void
  filter?: string
  position: { x: number; y: number }
  className?: string
}

/** Filter inline slash sections by label/id/description */
function filterInlineSections(sections: InlineSlashSection[], filter: string): InlineSlashSection[] {
  if (!filter) return sections
  const lowerFilter = filter.toLowerCase()

  return sections
    .map(section => ({
      ...section,
      items: section.items.filter(item =>
        item.label.toLowerCase().includes(lowerFilter) ||
        item.id.toLowerCase().includes(lowerFilter) ||
        item.description?.toLowerCase().includes(lowerFilter)
      ),
    }))
    .filter(section => section.items.length > 0)
}

/** Flatten inline slash sections into a single array */
function flattenInlineSections(sections: InlineSlashSection[]): SlashCommandItem[] {
  return sections.flatMap(section => section.items)
}

export function InlineSlashCommand({
  open,
  onOpenChange,
  sections,
  onSelectCommand,
  filter = '',
  position,
  className,
}: InlineSlashCommandProps) {
  const t = useT()
  const menuRef = React.useRef<HTMLDivElement>(null)
  const listRef = React.useRef<HTMLDivElement>(null)
  const [selectedIndex, setSelectedIndex] = React.useState(0)
  const filteredSections = filterInlineSections(sections, filter)
  const flatItems = flattenInlineSections(filteredSections)

  // Reset selection when filter changes
  React.useEffect(() => {
    setSelectedIndex(0)
  }, [filter])

  // Scroll selected item into view
  React.useEffect(() => {
    if (!listRef.current) return
    const selectedEl = listRef.current.querySelector('[data-selected="true"]')
    if (selectedEl) {
      selectedEl.scrollIntoView({ block: 'nearest' })
    }
  }, [selectedIndex])

  // Handle item selection
  const handleSelect = React.useCallback((item: SlashCommandItem) => {
    onSelectCommand(item.id, item.hasArgs)
    onOpenChange(false)
  }, [onSelectCommand, onOpenChange])

  // Keyboard navigation
  React.useEffect(() => {
    if (!open || flatItems.length === 0) return

    const handleKeyDown = (e: KeyboardEvent) => {
      switch (e.key) {
        case 'ArrowDown':
          e.preventDefault()
          setSelectedIndex(prev => (prev < flatItems.length - 1 ? prev + 1 : 0))
          break
        case 'ArrowUp':
          e.preventDefault()
          setSelectedIndex(prev => (prev > 0 ? prev - 1 : flatItems.length - 1))
          break
        case 'Enter':
        case 'Tab':
          e.preventDefault()
          if (flatItems[selectedIndex]) {
            handleSelect(flatItems[selectedIndex])
          }
          break
        case 'Escape':
          e.preventDefault()
          onOpenChange(false)
          break
      }
    }

    document.addEventListener('keydown', handleKeyDown)
    return () => document.removeEventListener('keydown', handleKeyDown)
  }, [open, flatItems, selectedIndex, handleSelect, onOpenChange])

  // Close on click outside
  React.useEffect(() => {
    if (!open) return

    const handleClickOutside = (e: MouseEvent) => {
      if (menuRef.current && !menuRef.current.contains(e.target as Node)) {
        onOpenChange(false)
      }
    }

    document.addEventListener('mousedown', handleClickOutside)
    return () => document.removeEventListener('mousedown', handleClickOutside)
  }, [open, onOpenChange])

  // Hide if no results or not open
  if (!open || flatItems.length === 0) return null

  // Calculate bottom position from window height (menu appears above cursor)
  const bottomPosition = typeof window !== 'undefined'
    ? window.innerHeight - Math.round(position.y) + 8
    : 0

  // Track current item index across all sections
  let currentItemIndex = 0

  return (
    <div
      ref={menuRef}
      className={cn('fixed z-dropdown', MENU_CONTAINER_STYLE, className)}
      style={{ left: Math.round(position.x) - 10, bottom: bottomPosition, minWidth: 260, maxWidth: 380 }}
    >
      {/* Menu header */}
      <div className={MENU_SECTION_HEADER}>
        {t('斜杠命令')}
      </div>

      <div ref={listRef} className={MENU_LIST_STYLE}>
        {filteredSections.map((section) => (
          <React.Fragment key={section.id}>
            {/* Show section label only when there are multiple sections */}
            {filteredSections.length > 1 && (
              <div className="px-3 py-1 text-[11px] font-medium text-muted-foreground/60 uppercase tracking-wider">
                {section.label}
              </div>
            )}

            {/* Section items */}
            {section.items.map((item) => {
              const itemIndex = currentItemIndex++
              const isSelected = itemIndex === selectedIndex

              return (
                <div
                  key={item.id}
                  data-selected={isSelected}
                  onClick={() => handleSelect(item)}
                  onMouseEnter={() => setSelectedIndex(itemIndex)}
                  className={cn(
                    MENU_ITEM_STYLE,
                    'items-start py-1.5',
                    isSelected && MENU_ITEM_SELECTED
                  )}
                >
                  {/* Terminal icon */}
                  <div className="shrink-0 text-muted-foreground mt-0.5">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                      <polyline points="4 17 10 11 4 5" />
                      <line x1="12" y1="19" x2="20" y2="19" />
                    </svg>
                  </div>
                  <div className="flex-1 min-w-0">
                    <span className="truncate block">{item.label}</span>
                    {item.description && (
                      <span className="truncate block text-[11px] text-muted-foreground/60 leading-tight">{item.description}</span>
                    )}
                  </div>
                  <span className="text-[11px] text-muted-foreground/50 shrink-0 mt-0.5">/{item.id}</span>
                </div>
              )
            })}
          </React.Fragment>
        ))}
      </div>
    </div>
  )
}

// ============================================================================
// Hook for managing inline slash command state (SDK + plugin commands)
// ============================================================================

/** Interface for elements that can be used with useInlineSlashCommand */
export interface SlashCommandInputElement {
  getBoundingClientRect: () => DOMRect
  getCaretRect?: () => DOMRect | null
  value: string
  selectionStart: number
}

export interface UseInlineSlashCommandOptions {
  /** Ref to input element (textarea or RichTextInput handle) */
  inputRef: React.RefObject<SlashCommandInputElement | null>
}

export interface UseInlineSlashCommandReturn {
  isOpen: boolean
  filter: string
  position: { x: number; y: number }
  sections: InlineSlashSection[]
  handleInputChange: (value: string, cursorPosition: number) => void
  close: () => void
  handleSelectCommand: (commandName: string, hasArgs: boolean) => { value: string; cursorPosition: number }
}

export function useInlineSlashCommand({
  inputRef,
}: UseInlineSlashCommandOptions): UseInlineSlashCommandReturn {
  const t = useT()
  const [isOpen, setIsOpen] = React.useState(false)
  const [filter, setFilter] = React.useState('')
  const [position, setPosition] = React.useState({ x: 0, y: 0 })
  const [slashStart, setSlashStart] = React.useState(-1)
  // Store current input state for handleSelect
  const currentInputRef = React.useRef({ value: '', cursorPosition: 0 })

  // Read SDK + plugin slash commands from atom
  const sdkSlashCommands = useAtomValue(sdkSlashCommandsAtom)
  const commandTranslations = useAtomValue(commandTranslationsAtom)

  // Build sections from SDK built-in commands + plugin commands
  const sections = React.useMemo((): InlineSlashSection[] => {
    const builtinNames = new Set(Object.keys(SDK_COMMAND_TRANSLATIONS))

    // SDK built-in commands (always available)
    const builtinItems: SlashCommandItem[] = Object.entries(SDK_COMMAND_TRANSLATIONS).map(([name, trans]) => ({
      id: name,
      label: t(trans.label),
      description: t(trans.description),
      hasArgs: false,
    }))

    // Plugin commands from atom (available after IPC event)
    const pluginItems: SlashCommandItem[] = sdkSlashCommands
      .filter(cmd => !builtinNames.has(cmd.name))
      .map(cmd => {
        const display = getCommandDisplay(cmd, commandTranslations)
        return {
          id: cmd.name,
          label: t(display.label),
          description: t(display.description),
          hasArgs: display.hasArgs,
        }
      })

    const result: InlineSlashSection[] = []

    // Built-in commands section
    if (builtinItems.length > 0) {
      result.push({
        id: 'builtin',
        label: t('内置命令'),
        items: builtinItems.sort((a, b) => a.label.localeCompare(b.label)),
      })
    }

    // Plugin commands section
    if (pluginItems.length > 0) {
      result.push({
        id: 'plugins',
        label: t('插件命令'),
        items: pluginItems.sort((a, b) => a.label.localeCompare(b.label)),
      })
    }

    return result
  }, [sdkSlashCommands, commandTranslations, t])

  const handleInputChange = React.useCallback((value: string, cursorPosition: number) => {
    // Store current state for handleSelect
    currentInputRef.current = { value, cursorPosition }

    const textBeforeCursor = value.slice(0, cursorPosition)
    const slashMatch = textBeforeCursor.match(/(?:^|\s)\/(\w[\w\-:]*)?$/)

    // Only show menu if we have sections with items
    const hasItems = sections.some(s => s.items.length > 0)

    if (slashMatch && hasItems) {
      const filterText = slashMatch[1] || ''
      // Check if there are any filtered results before opening menu
      const filtered = filterInlineSections(sections, filterText)
      const hasFilteredItems = filtered.some(s => s.items.length > 0)

      if (!hasFilteredItems) {
        setIsOpen(false)
        setFilter('')
        setSlashStart(-1)
        return
      }

      const matchStart = textBeforeCursor.lastIndexOf('/')
      setSlashStart(matchStart)
      setFilter(filterText)

      if (inputRef.current) {
        const caretRect = inputRef.current.getCaretRect?.()

        if (caretRect && caretRect.x > 0) {
          setPosition({
            x: caretRect.x,
            y: caretRect.y,
          })
        } else {
          const rect = inputRef.current.getBoundingClientRect()
          const lineHeight = 20
          const linesBeforeCursor = textBeforeCursor.split('\n').length - 1
          setPosition({
            x: rect.left,
            y: rect.top + (linesBeforeCursor + 1) * lineHeight,
          })
        }
      }

      setIsOpen(true)
    } else {
      setIsOpen(false)
      setFilter('')
      setSlashStart(-1)
    }
  }, [inputRef, sections])

  const handleSelectCommand = React.useCallback((commandName: string, hasArgs: boolean): { value: string; cursorPosition: number } => {
    let result = ''
    let newCursorPosition = 0

    if (slashStart >= 0) {
      const { value: currentValue, cursorPosition } = currentInputRef.current
      const before = currentValue.slice(0, slashStart)
      const after = currentValue.slice(cursorPosition)

      if (hasArgs) {
        // Has args: insert /{name} and let user type arguments
        const commandText = `/${commandName} `
        result = before + commandText + after
        newCursorPosition = before.length + commandText.length
      } else {
        // No args: replace with /{name} for direct execution
        result = `/${commandName}`
        newCursorPosition = result.length
      }
    }

    setIsOpen(false)
    return { value: result, cursorPosition: newCursorPosition }
  }, [slashStart])

  const close = React.useCallback(() => {
    setIsOpen(false)
    setFilter('')
    setSlashStart(-1)
  }, [])

  return {
    isOpen,
    filter,
    position,
    sections,
    handleInputChange,
    close,
    handleSelectCommand,
  }
}
