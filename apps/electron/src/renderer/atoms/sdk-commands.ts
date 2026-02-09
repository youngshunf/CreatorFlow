/**
 * SDK Slash Commands Atom
 *
 * Stores SDK slash commands received from the agent (via IPC).
 * Used by the @ mention menu to display available commands.
 */

import { atom } from 'jotai'
import type { SdkSlashCommand } from '../../shared/types'

/**
 * Atom to store SDK slash commands for the current workspace.
 * Populated when the agent sends 'slash_commands_available' event.
 * Read by useInlineMention hook to build the Commands section.
 */
export const sdkSlashCommandsAtom = atom<SdkSlashCommand[]>([])

/**
 * Atom to store plugin command translations for the current workspace.
 * Populated alongside sdkSlashCommandsAtom from the same IPC event.
 */
export const commandTranslationsAtom = atom<Record<string, { label: string; description: string }>>({})
