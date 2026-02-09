/**
 * Local File Protocol Handler
 *
 * Registers a custom `localfile://` protocol that serves local files
 * directly to the renderer process. This enables <video>, <audio>, and
 * other HTML elements to load local files without base64 encoding overhead.
 *
 * URL format: localfile://file/<encodeURIComponent(absolutePath)>
 * Examples:
 *   macOS:   localfile://file/%2FUsers%2Ffoo%2Fvideo.mp4
 *   Windows: localfile://file/C%3A%5CUsers%5Cfoo%5Cvideo.mp4
 */

import { protocol } from 'electron'
import { stat } from 'fs/promises'
import { createReadStream } from 'fs'
import { isAbsolute, extname } from 'path'
import { mainLog } from './logger'

/** MIME type mapping for common media files */
const MIME_TYPES: Record<string, string> = {
  // Video
  '.mp4': 'video/mp4',
  '.webm': 'video/webm',
  '.ogg': 'video/ogg',
  '.ogv': 'video/ogg',
  '.mov': 'video/quicktime',
  '.avi': 'video/x-msvideo',
  '.mkv': 'video/x-matroska',
  '.m4v': 'video/mp4',
  // Audio
  '.mp3': 'audio/mpeg',
  '.wav': 'audio/wav',
  '.flac': 'audio/flac',
  '.aac': 'audio/aac',
  '.m4a': 'audio/mp4',
  '.wma': 'audio/x-ms-wma',
  '.oga': 'audio/ogg',
  // Images
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.gif': 'image/gif',
  '.webp': 'image/webp',
  '.svg': 'image/svg+xml',
  '.bmp': 'image/bmp',
  '.ico': 'image/x-icon',
  // Documents
  '.pdf': 'application/pdf',
}

/**
 * Register the localfile:// custom protocol scheme.
 * MUST be called before app.whenReady().
 */
export function registerLocalFileScheme(): void {
  protocol.registerSchemesAsPrivileged([
    {
      scheme: 'localfile',
      privileges: {
        supportFetchAPI: true,
        standard: true,
        corsEnabled: true,
        stream: true,
      },
    },
  ])
}

/**
 * Register the localfile:// protocol handler.
 * Must be called after app.whenReady().
 */
export function registerLocalFileHandler(): void {
  protocol.handle('localfile', async (request) => {
    try {
      mainLog.info('Local file protocol request:', request.url)
      const url = new URL(request.url)
      mainLog.info('Local file protocol parsed - host:', url.host, 'pathname:', url.pathname)
      const filePath = decodeURIComponent(url.pathname.slice(1))
      mainLog.info('Local file protocol resolved path:', filePath)

      if (!filePath || !isAbsolute(filePath)) {
        mainLog.warn('Local file protocol: invalid path', filePath)
        return new Response(null, { status: 400 })
      }

      // Get file info
      let fileSize: number
      try {
        const fileStat = await stat(filePath)
        fileSize = fileStat.size
      } catch {
        return new Response(null, { status: 404 })
      }

      const ext = extname(filePath).toLowerCase()
      const mimeType = MIME_TYPES[ext] || 'application/octet-stream'

      // Handle Range requests (needed for video/audio seeking)
      const rangeHeader = request.headers.get('Range')
      if (rangeHeader) {
        const match = rangeHeader.match(/bytes=(\d+)-(\d*)/)
        if (match) {
          const start = parseInt(match[1], 10)
          const end = match[2] ? parseInt(match[2], 10) : fileSize - 1
          const chunkSize = end - start + 1

          const stream = createReadStream(filePath, { start, end })
          const readable = new ReadableStream({
            start(controller) {
              stream.on('data', (chunk: Buffer) => controller.enqueue(new Uint8Array(chunk)))
              stream.on('end', () => controller.close())
              stream.on('error', (err) => controller.error(err))
            },
            cancel() {
              stream.destroy()
            },
          })

          return new Response(readable, {
            status: 206,
            headers: {
              'Content-Type': mimeType,
              'Content-Range': `bytes ${start}-${end}/${fileSize}`,
              'Content-Length': String(chunkSize),
              'Accept-Ranges': 'bytes',
            },
          })
        }
      }

      // Full file response
      const stream = createReadStream(filePath)
      const readable = new ReadableStream({
        start(controller) {
          stream.on('data', (chunk: Buffer) => controller.enqueue(new Uint8Array(chunk)))
          stream.on('end', () => controller.close())
          stream.on('error', (err) => controller.error(err))
        },
        cancel() {
          stream.destroy()
        },
      })

      return new Response(readable, {
        status: 200,
        headers: {
          'Content-Type': mimeType,
          'Content-Length': String(fileSize),
          'Accept-Ranges': 'bytes',
        },
      })
    } catch (error) {
      mainLog.error('Local file protocol error:', error)
      return new Response(null, { status: 500 })
    }
  })

  mainLog.info('Registered localfile:// protocol handler')
}
