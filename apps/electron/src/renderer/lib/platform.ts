/**
 * Platform Detection Utilities
 *
 * Centralized platform detection for the renderer process.
 * Use these instead of accessing navigator.platform directly.
 *
 * @example
 * import { isMac, isWindows, PATH_SEP, getPathBasename } from '@/lib/platform'
 *
 * // Platform checks
 * const modifier = isMac ? '⌘' : 'Ctrl'
 *
 * // Path handling
 * const folderName = getPathBasename('/Users/alice/projects') // 'projects'
 */

/** True if running on macOS */
export const isMac =
  typeof navigator !== 'undefined' &&
  navigator.platform.toLowerCase().includes('mac')

/** True if running on Windows */
export const isWindows =
  typeof navigator !== 'undefined' &&
  navigator.platform.toLowerCase().includes('win')

/** True if running on Linux */
export const isLinux =
  typeof navigator !== 'undefined' &&
  navigator.platform.toLowerCase().includes('linux')

/** Native path separator for current OS */
export const PATH_SEP = isWindows ? '\\' : '/'

/**
 * Get the last segment of a path (folder/file name).
 * Handles both Unix (/) and Windows (\\) separators based on current OS.
 * Handles trailing slashes and empty segments.
 */
export function getPathBasename(path: string): string {
  if (!path) return ''
  
  // Remove trailing slash(es) to handle paths like "/foo/bar/"
  let normalized = path
  while (normalized.endsWith(PATH_SEP) && normalized.length > 1) {
    normalized = normalized.slice(0, -1)
  }
  
  // Handle root paths (just "/" or "C:\")
  if (normalized === PATH_SEP || /^[A-Z]:\\?$/.test(normalized)) {
    return normalized
  }
  
  // Split and get last non-empty segment
  const segments = normalized.split(PATH_SEP).filter(s => s.length > 0)
  return segments[segments.length - 1] || ''
}
