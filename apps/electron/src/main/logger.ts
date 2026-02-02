import log from 'electron-log/main'
import { app } from 'electron'

/**
 * Debug mode is enabled when running from source (not packaged) or with --debug flag.
 * - true: `bun run electron:start` or `electron .` or packaged app with `--debug`
 * - false: bundled .app/.exe release without --debug flag
 */
export const isDebugMode = !app.isPackaged || process.argv.includes('--debug')

// Helper: format timestamp in local time (no 8h UTC offset)
function formatLocalTimestamp(date: Date): string {
  const pad = (n: number, width = 2) => n.toString().padStart(width, '0')
  const year = date.getFullYear()
  const month = pad(date.getMonth() + 1)
  const day = pad(date.getDate())
  const hours = pad(date.getHours())
  const minutes = pad(date.getMinutes())
  const seconds = pad(date.getSeconds())
  const millis = pad(date.getMilliseconds(), 3)
  return `${year}-${month}-${day}T${hours}:${minutes}:${seconds}.${millis}`
}

// Configure transports based on debug mode
if (isDebugMode) {
  // JSON format for file (agent-parseable)
  // Note: format expects (params: FormatParams) => any[], where params.message has the LogMessage fields
  log.transports.file.format = ({ message }) => [
    JSON.stringify({
      // Use local time instead of UTC to avoid 8-hour offset when inspecting logs
      timestamp: formatLocalTimestamp(message.date),
      level: message.level,
      scope: message.scope,
      message: message.data,
    }),
  ]

  log.transports.file.maxSize = 5 * 1024 * 1024 // 5MB

  // Console output in debug mode with readable format
  // Note: format must return an array - electron-log's transformStyles calls .reduce() on it
  log.transports.console.format = ({ message }) => {
    const scope = message.scope ? `[${message.scope}]` : ''
    const level = message.level.toUpperCase().padEnd(5)
    const data = message.data
      .map((d: unknown) => (typeof d === 'object' ? JSON.stringify(d) : String(d)))
      .join(' ')
    return [`${formatLocalTimestamp(message.date)} ${level} ${scope} ${data}`]
  }
  log.transports.console.level = 'debug'
} else {
  // Production: still log to file for debugging, but no console
  log.transports.file.format = ({ message }) => [
    JSON.stringify({
      timestamp: formatLocalTimestamp(message.date),
      level: message.level,
      scope: message.scope,
      message: message.data,
    }),
  ]
  log.transports.file.maxSize = 5 * 1024 * 1024 // 5MB
  log.transports.file.level = 'info' // Log info and above in production
  log.transports.console.level = false
}

// Export scoped loggers for different modules
export const mainLog = log.scope('main')
export const sessionLog = log.scope('session')
export const ipcLog = log.scope('ipc')
export const windowLog = log.scope('window')
export const agentLog = log.scope('agent')
export const searchLog = log.scope('search')

/**
 * Get the path to the current log file.
 * Returns undefined if file logging is disabled.
 */
export function getLogFilePath(): string | undefined {
  if (!isDebugMode) return undefined
  return log.transports.file.getFile()?.path
}

export default log
