/**
 * File Manager IPC handlers
 * Provides file browsing, creation, deletion, rename, move, copy operations
 * for the workspace file manager feature.
 */
import { ipcMain, BrowserWindow } from 'electron'
import * as fs from 'fs/promises'
import { watch, type FSWatcher } from 'fs'
import * as path from 'path'
import { IPC_CHANNELS, type FMFileEntry, type FMFileInfo, type FMDirectoryChangeEvent } from '../shared/types'
import { ipcLog } from './logger'

// Active directory watchers (path -> watcher)
const watchers = new Map<string, FSWatcher>()

/**
 * Register file manager IPC handlers
 */
export function registerFileManagerIpc(): void {
  // List directory contents
  ipcMain.handle(IPC_CHANNELS.FM_LIST_DIRECTORY, async (_, dirPath: string): Promise<FMFileEntry[]> => {
    try {
      const entries = await fs.readdir(dirPath, { withFileTypes: true })
      const result: FMFileEntry[] = []

      for (const entry of entries) {
        // Skip hidden files (starting with .)
        if (entry.name.startsWith('.')) continue

        const fullPath = path.join(dirPath, entry.name)
        
        try {
          const stats = await fs.stat(fullPath)
          
          result.push({
            path: fullPath,
            name: entry.name,
            isDirectory: entry.isDirectory(),
            size: entry.isDirectory() ? undefined : stats.size,
            modifiedTime: stats.mtimeMs,
            extension: entry.isDirectory() ? undefined : path.extname(entry.name).slice(1).toLowerCase(),
          })
        } catch (statErr) {
          // Skip files we can't stat (permission denied, broken symlinks, etc.)
          ipcLog.info(`Skipping file ${fullPath}: ${statErr}`)
        }
      }

      // Sort: directories first, then by name (case-insensitive)
      return result.sort((a, b) => {
        if (a.isDirectory !== b.isDirectory) return a.isDirectory ? -1 : 1
        return a.name.localeCompare(b.name, undefined, { sensitivity: 'base' })
      })
    } catch (err) {
      ipcLog.error('fm:listDirectory error:', err)
      throw err
    }
  })

  // Create a new folder
  ipcMain.handle(IPC_CHANNELS.FM_CREATE_FOLDER, async (_, fullPath: string): Promise<string> => {
    try {
      await fs.mkdir(fullPath, { recursive: true })
      ipcLog.info(`Created folder: ${fullPath}`)
      return fullPath
    } catch (err) {
      ipcLog.error('fm:createFolder error:', err)
      throw err
    }
  })

  // Delete files/folders
  ipcMain.handle(IPC_CHANNELS.FM_DELETE, async (_, paths: string[]): Promise<void> => {
    try {
      for (const p of paths) {
        await fs.rm(p, { recursive: true, force: true })
        ipcLog.info(`Deleted: ${p}`)
      }
    } catch (err) {
      ipcLog.error('fm:delete error:', err)
      throw err
    }
  })

  // Rename a file/folder
  ipcMain.handle(IPC_CHANNELS.FM_RENAME, async (_, oldPath: string, newName: string): Promise<string> => {
    try {
      const dir = path.dirname(oldPath)
      const newPath = path.join(dir, newName)
      await fs.rename(oldPath, newPath)
      ipcLog.info(`Renamed: ${oldPath} -> ${newPath}`)
      return newPath
    } catch (err) {
      ipcLog.error('fm:rename error:', err)
      throw err
    }
  })

  // Move files/folders to a new directory
  ipcMain.handle(IPC_CHANNELS.FM_MOVE, async (_, srcPaths: string[], destDir: string): Promise<void> => {
    try {
      for (const src of srcPaths) {
        const dest = path.join(destDir, path.basename(src))
        await fs.rename(src, dest)
        ipcLog.info(`Moved: ${src} -> ${dest}`)
      }
    } catch (err) {
      ipcLog.error('fm:move error:', err)
      throw err
    }
  })

  // Copy files/folders to a new directory
  ipcMain.handle(IPC_CHANNELS.FM_COPY, async (_, srcPaths: string[], destDir: string): Promise<void> => {
    try {
      for (const src of srcPaths) {
        const dest = path.join(destDir, path.basename(src))
        await fs.cp(src, dest, { recursive: true })
        ipcLog.info(`Copied: ${src} -> ${dest}`)
      }
    } catch (err) {
      ipcLog.error('fm:copy error:', err)
      throw err
    }
  })

  // Get detailed file info
  ipcMain.handle(IPC_CHANNELS.FM_GET_FILE_INFO, async (_, filePath: string): Promise<FMFileInfo> => {
    try {
      const stats = await fs.stat(filePath)
      return {
        path: filePath,
        name: path.basename(filePath),
        size: stats.size,
        isDir: stats.isDirectory(),
        createdAt: stats.birthtime,
        modifiedAt: stats.mtime,
      }
    } catch (err) {
      ipcLog.error('fm:getFileInfo error:', err)
      throw err
    }
  })

  // Read file as base64 (for binary files like images)
  ipcMain.handle(IPC_CHANNELS.FM_READ_FILE_BASE64, async (_, filePath: string, maxSize?: number): Promise<{ data: string; mimeType: string }> => {
    try {
      const stats = await fs.stat(filePath)
      const limit = maxSize || 10 * 1024 * 1024 // Default 10MB limit

      if (stats.size > limit) {
        throw new Error(`File too large: ${stats.size} bytes (max ${limit})`)
      }

      const buffer = await fs.readFile(filePath)
      const base64 = buffer.toString('base64')

      // Determine MIME type from extension
      const ext = path.extname(filePath).toLowerCase().slice(1)
      const mimeTypes: Record<string, string> = {
        // Images
        png: 'image/png',
        jpg: 'image/jpeg',
        jpeg: 'image/jpeg',
        gif: 'image/gif',
        webp: 'image/webp',
        svg: 'image/svg+xml',
        ico: 'image/x-icon',
        bmp: 'image/bmp',
        // Documents
        pdf: 'application/pdf',
        // Audio
        mp3: 'audio/mpeg',
        wav: 'audio/wav',
        ogg: 'audio/ogg',
        // Video
        mp4: 'video/mp4',
        webm: 'video/webm',
      }

      const mimeType = mimeTypes[ext] || 'application/octet-stream'

      return { data: base64, mimeType }
    } catch (err) {
      ipcLog.error('fm:readFileBase64 error:', err)
      throw err
    }
  })

  // Write text file content
  ipcMain.handle(IPC_CHANNELS.FM_WRITE_FILE, async (_, filePath: string, content: string): Promise<void> => {
    try {
      // 允许写入的文件扩展名（文本/代码类文件）
      const ext = path.extname(filePath).toLowerCase().slice(1)
      const allowedExtensions = [
        // 文本
        'txt', 'md', 'markdown', 'mdx', 'log', 'csv', 'ini', 'conf', 'cfg',
        // 代码
        'ts', 'tsx', 'js', 'jsx', 'json', 'yaml', 'yml', 'py', 'rb', 'go', 'rs',
        'java', 'c', 'cpp', 'h', 'css', 'scss', 'less', 'html', 'xml', 'sh', 'bash',
        'zsh', 'vue', 'svelte', 'swift', 'kt', 'php', 'sql', 'graphql', 'toml', 'env',
        'dockerfile', 'makefile',
      ]

      if (!allowedExtensions.includes(ext)) {
        throw new Error(`不允许写入此类型的文件: .${ext}`)
      }

      await fs.writeFile(filePath, content, 'utf-8')
      ipcLog.info(`Wrote file: ${filePath}`)
    } catch (err) {
      ipcLog.error('fm:writeFile error:', err)
      throw err
    }
  })

  // Start watching a directory for changes
  ipcMain.on(IPC_CHANNELS.FM_WATCH_DIRECTORY, (_, dirPath: string) => {
    // Don't create duplicate watchers
    if (watchers.has(dirPath)) {
      ipcLog.info(`Already watching: ${dirPath}`)
      return
    }

    try {
      // Use Node.js native fs.watch with recursive option (macOS/Windows support)
      const watcher = watch(dirPath, { recursive: false }, (eventType, filename) => {
        if (!filename) return
        
        // Map Node.js event types to our event types
        const event: FMDirectoryChangeEvent = {
          type: eventType === 'rename' ? 'change' : 'change',
          path: path.join(dirPath, filename),
        }
        
        // Broadcast to all windows
        for (const win of BrowserWindow.getAllWindows()) {
          if (!win.isDestroyed() && !win.webContents.isDestroyed()) {
            win.webContents.send(IPC_CHANNELS.FM_DIRECTORY_CHANGED, event)
          }
        }
      })

      watcher.on('error', (err: Error) => {
        ipcLog.error(`Watcher error for ${dirPath}:`, err)
      })

      watchers.set(dirPath, watcher)
      ipcLog.info(`Started watching: ${dirPath}`)
    } catch (err) {
      ipcLog.error('fm:watchDirectory error:', err)
    }
  })

  // Stop watching a directory
  ipcMain.on(IPC_CHANNELS.FM_UNWATCH_DIRECTORY, (_, dirPath: string) => {
    const watcher = watchers.get(dirPath)
    if (watcher) {
      watcher.close()
      watchers.delete(dirPath)
      ipcLog.info(`Stopped watching: ${dirPath}`)
    }
  })
}

/**
 * Clean up all watchers (call on app quit)
 */
export function cleanupFileManagerWatchers(): void {
  for (const [dirPath, watcher] of watchers) {
    watcher.close()
    ipcLog.info(`Cleaned up watcher: ${dirPath}`)
  }
  watchers.clear()
}
