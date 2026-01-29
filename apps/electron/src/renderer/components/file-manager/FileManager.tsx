/**
 * FileManager - File browser component for browsing and managing files
 *
 * Features:
 * - Directory navigation with breadcrumb trail
 * - File/folder listing with icons and metadata
 * - Create, rename, delete, copy, move operations
 * - Real-time directory watching for auto-refresh
 * - Context menu for file operations
 * - Double-click to open files/navigate folders
 */

import * as React from 'react'
import { useState, useEffect, useCallback, useRef, useMemo } from 'react'
import { AnimatePresence, motion, type Variants } from 'motion/react'
import {
  File,
  Folder,
  FolderOpen,
  FileText,
  Image,
  FileCode,
  ChevronRight,
  Home,
  RefreshCw,
  FolderPlus,
  Trash2,
  Copy,
  Scissors,
  ClipboardPaste,
  Pencil,
  ArrowUp,
  HardDrive,
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ScrollArea } from '@/components/ui/scroll-area'
import {
  ContextMenu,
  ContextMenuTrigger,
  StyledContextMenuContent,
} from '@/components/ui/styled-context-menu'
import { ContextMenuProvider, useMenuComponents } from '@/components/ui/menu-context'
import { useT } from '@/context/LocaleContext'
import type { FMFileEntry, FMDirectoryChangeEvent } from '../../../shared/types'

/**
 * Stagger animation variants for items
 */
const containerVariants: Variants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: {
      staggerChildren: 0.02,
      delayChildren: 0.01,
    },
  },
}

const itemVariants: Variants = {
  hidden: { opacity: 0, y: 4 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.12, ease: 'easeOut' },
  },
}

export interface FileManagerProps {
  /** Initial directory path */
  initialPath?: string
  /** Root path that user cannot navigate above (workspace root) */
  rootPath?: string
  /** Called when a file is selected (single click) */
  onFileSelect?: (file: FMFileEntry) => void
  /** Called when a file is opened (double click) */
  onFileOpen?: (file: FMFileEntry) => void
  className?: string
}

/**
 * Format file size in human-readable format
 */
function formatFileSize(bytes?: number): string {
  if (bytes === undefined || bytes < 0) return ''
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(1)} GB`
}

/**
 * Format date in locale-friendly format
 */
function formatDate(timestamp?: number): string {
  if (!timestamp) return ''
  return new Date(timestamp).toLocaleDateString(undefined, {
    year: 'numeric',
    month: 'short',
    day: 'numeric',
    hour: '2-digit',
    minute: '2-digit',
  })
}

/**
 * Get icon for file based on name/extension
 */
function getFileIcon(file: FMFileEntry, isSelected?: boolean) {
  const iconClass = cn(
    "h-4 w-4 shrink-0",
    isSelected ? "text-primary-foreground" : "text-muted-foreground"
  )

  if (file.isDirectory) {
    return <Folder className={iconClass} />
  }

  const ext = file.name.split('.').pop()?.toLowerCase()

  // Document types
  if (['md', 'markdown', 'txt', 'doc', 'docx', 'pdf'].includes(ext || '')) {
    return <FileText className={iconClass} />
  }

  // Image types
  if (['png', 'jpg', 'jpeg', 'gif', 'webp', 'svg', 'ico', 'bmp'].includes(ext || '')) {
    return <Image className={iconClass} />
  }

  // Code types
  if (['ts', 'tsx', 'js', 'jsx', 'json', 'yaml', 'yml', 'py', 'rb', 'go', 'rs', 'java', 'c', 'cpp', 'h', 'css', 'scss', 'html', 'xml', 'sh', 'bash', 'zsh'].includes(ext || '')) {
    return <FileCode className={iconClass} />
  }

  return <File className={iconClass} />
}

/**
 * Context menu content for file operations
 */
function FileContextMenu({
  file,
  onRename,
  onDelete,
  onCopy,
  onCut,
  canPaste,
  onPaste,
}: {
  file: FMFileEntry
  onRename: () => void
  onDelete: () => void
  onCopy: () => void
  onCut: () => void
  canPaste: boolean
  onPaste: () => void
}) {
  const t = useT()
  const { MenuItem, Separator } = useMenuComponents()

  return (
    <>
      <MenuItem onClick={onRename}>
        <Pencil className="h-3.5 w-3.5" />
        <span className="flex-1">{t('重命名')}</span>
      </MenuItem>
      <MenuItem onClick={onCopy}>
        <Copy className="h-3.5 w-3.5" />
        <span className="flex-1">{t('复制')}</span>
      </MenuItem>
      <MenuItem onClick={onCut}>
        <Scissors className="h-3.5 w-3.5" />
        <span className="flex-1">{t('剪切')}</span>
      </MenuItem>
      {canPaste && (
        <MenuItem onClick={onPaste}>
          <ClipboardPaste className="h-3.5 w-3.5" />
          <span className="flex-1">{t('粘贴')}</span>
        </MenuItem>
      )}
      <Separator />
      <MenuItem onClick={onDelete}>
        <Trash2 className="h-3.5 w-3.5" />
        <span className="flex-1">{t('删除')}</span>
      </MenuItem>
    </>
  )
}

/**
 * Context menu for empty space (create folder, paste)
 */
function EmptyContextMenu({
  onCreateFolder,
  canPaste,
  onPaste,
}: {
  onCreateFolder: () => void
  canPaste: boolean
  onPaste: () => void
}) {
  const t = useT()
  const { MenuItem, Separator } = useMenuComponents()

  return (
    <>
      <MenuItem onClick={onCreateFolder}>
        <FolderPlus className="h-3.5 w-3.5" />
        <span className="flex-1">{t('新建文件夹')}</span>
      </MenuItem>
      {canPaste && (
        <>
          <Separator />
          <MenuItem onClick={onPaste}>
            <ClipboardPaste className="h-3.5 w-3.5" />
            <span className="flex-1">{t('粘贴')}</span>
          </MenuItem>
        </>
      )}
    </>
  )
}

/**
 * Single file/folder row item
 */
function FileItem({
  file,
  isSelected,
  isRenaming,
  renameValue,
  onRenameValueChange,
  onRenameSubmit,
  onRenameCancel,
  onClick,
  onDoubleClick,
  onRename,
  onDelete,
  onCopy,
  onCut,
  canPaste,
  onPaste,
}: {
  file: FMFileEntry
  isSelected: boolean
  isRenaming: boolean
  renameValue: string
  onRenameValueChange: (value: string) => void
  onRenameSubmit: () => void
  onRenameCancel: () => void
  onClick: () => void
  onDoubleClick: () => void
  onRename: () => void
  onDelete: () => void
  onCopy: () => void
  onCut: () => void
  canPaste: boolean
  onPaste: () => void
}) {
  const renameInputRef = useRef<HTMLInputElement>(null)

  // Focus input when renaming starts
  useEffect(() => {
    if (isRenaming && renameInputRef.current) {
      renameInputRef.current.focus()
      // Select filename without extension
      const lastDot = renameValue.lastIndexOf('.')
      if (lastDot > 0 && !file.isDirectory) {
        renameInputRef.current.setSelectionRange(0, lastDot)
      } else {
        renameInputRef.current.select()
      }
    }
  }, [isRenaming, renameValue, file.isDirectory])

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      e.preventDefault()
      onRenameSubmit()
    } else if (e.key === 'Escape') {
      e.preventDefault()
      onRenameCancel()
    }
  }

  return (
    <ContextMenu modal={true}>
      <ContextMenuTrigger asChild>
        <motion.div
          variants={itemVariants}
          className={cn(
            "flex items-center gap-3 px-3 py-2 rounded-md cursor-pointer select-none",
            "hover:bg-muted/50 transition-colors",
            isSelected && "bg-primary text-primary-foreground hover:bg-primary/90"
          )}
          onClick={onClick}
          onDoubleClick={onDoubleClick}
        >
          {getFileIcon(file, isSelected)}
          
          {isRenaming ? (
            <Input
              ref={renameInputRef}
              value={renameValue}
              onChange={(e) => onRenameValueChange(e.target.value)}
              onKeyDown={handleKeyDown}
              onBlur={onRenameSubmit}
              className="h-6 py-0 px-1 text-sm flex-1"
              onClick={(e) => e.stopPropagation()}
              onDoubleClick={(e) => e.stopPropagation()}
            />
          ) : (
            <span className="flex-1 truncate text-sm">{file.name}</span>
          )}
          
          {!isRenaming && !file.isDirectory && (
            <span className={cn(
              "text-xs shrink-0",
              isSelected ? "text-primary-foreground/70" : "text-muted-foreground"
            )}>
              {formatFileSize(file.size)}
            </span>
          )}
          
          {!isRenaming && (
            <span className={cn(
              "text-xs shrink-0 hidden sm:block",
              isSelected ? "text-primary-foreground/70" : "text-muted-foreground"
            )}>
              {formatDate(file.modifiedTime)}
            </span>
          )}
        </motion.div>
      </ContextMenuTrigger>
      <StyledContextMenuContent>
        <ContextMenuProvider>
          <FileContextMenu
            file={file}
            onRename={onRename}
            onDelete={onDelete}
            onCopy={onCopy}
            onCut={onCut}
            canPaste={canPaste}
            onPaste={onPaste}
          />
        </ContextMenuProvider>
      </StyledContextMenuContent>
    </ContextMenu>
  )
}

/**
 * Breadcrumb navigation component
 */
function Breadcrumbs({
  path,
  rootPath,
  onNavigate,
}: {
  path: string
  rootPath?: string
  onNavigate: (path: string) => void
}) {
  const t = useT()
  
  // Split path into segments, starting from rootPath if provided
  const segments = useMemo(() => {
    const parts = path.split('/').filter(Boolean)
    const result: { name: string; path: string }[] = []
    
    let currentPath = ''
    for (const part of parts) {
      currentPath += '/' + part
      // Only include segments at or below rootPath
      if (!rootPath || currentPath.length >= rootPath.length) {
        result.push({ name: part, path: currentPath })
      }
    }
    
    return result
  }, [path, rootPath])

  // Get the root name for display
  const rootName = rootPath ? rootPath.split('/').filter(Boolean).pop() || t('工作区') : t('根目录')

  return (
    <div className="flex items-center gap-1 text-sm overflow-x-auto">
      <button
        onClick={() => onNavigate(rootPath || '/')}
        className="p-1 hover:bg-muted rounded shrink-0 flex items-center gap-1"
        title={rootPath || t('根目录')}
      >
        <FolderOpen className="h-4 w-4 text-muted-foreground" />
        {rootPath && <span className="text-xs text-muted-foreground">{rootName}</span>}
      </button>
      
      {segments.filter(seg => seg.path !== rootPath).map((segment, index, arr) => (
        <React.Fragment key={segment.path}>
          <ChevronRight className="h-3 w-3 text-muted-foreground shrink-0" />
          <button
            onClick={() => onNavigate(segment.path)}
            className={cn(
              "px-1.5 py-0.5 rounded hover:bg-muted truncate max-w-[120px]",
              index === arr.length - 1 && "font-medium"
            )}
            title={segment.path}
          >
            {segment.name}
          </button>
        </React.Fragment>
      ))}
    </div>
  )
}

/**
 * FileManager component
 */
export function FileManager({
  initialPath,
  rootPath,
  onFileSelect,
  onFileOpen,
  className,
}: FileManagerProps) {
  const t = useT()
  // Use rootPath as default if initialPath not provided
  const effectiveInitialPath = initialPath || rootPath || '/'
  const [currentPath, setCurrentPath] = useState(effectiveInitialPath)
  const [files, setFiles] = useState<FMFileEntry[]>([])
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedFile, setSelectedFile] = useState<FMFileEntry | null>(null)
  
  // Rename state
  const [renamingFile, setRenamingFile] = useState<FMFileEntry | null>(null)
  const [renameValue, setRenameValue] = useState('')
  
  // Clipboard state
  const [clipboard, setClipboard] = useState<{ file: FMFileEntry; operation: 'copy' | 'cut' } | null>(null)
  
  // Creating new folder state
  const [isCreatingFolder, setIsCreatingFolder] = useState(false)
  const [newFolderName, setNewFolderName] = useState('')
  const newFolderInputRef = useRef<HTMLInputElement>(null)
  
  const mountedRef = useRef(true)

  // Load directory contents
  const loadDirectory = useCallback(async (path: string) => {
    setIsLoading(true)
    setError(null)
    try {
      const entries = await window.electronAPI.fm.listDirectory(path)
      if (mountedRef.current) {
        // Sort: folders first, then by name
        const sorted = [...entries].sort((a, b) => {
          if (a.isDirectory && !b.isDirectory) return -1
          if (!a.isDirectory && b.isDirectory) return 1
          return a.name.localeCompare(b.name)
        })
        setFiles(sorted)
        setSelectedFile(null)
      }
    } catch (err) {
      console.error('Failed to load directory:', err)
      if (mountedRef.current) {
        setError(err instanceof Error ? err.message : String(err))
        setFiles([])
      }
    } finally {
      if (mountedRef.current) {
        setIsLoading(false)
      }
    }
  }, [])

  // Initial load and path change
  useEffect(() => {
    mountedRef.current = true
    loadDirectory(currentPath)
    
    return () => {
      mountedRef.current = false
    }
  }, [currentPath, loadDirectory])

  // Watch directory for changes
  useEffect(() => {
    window.electronAPI.fm.watchDirectory(currentPath)
    
    const unsubscribe = window.electronAPI.fm.onDirectoryChanged((event: FMDirectoryChangeEvent) => {
      if (event.path === currentPath && mountedRef.current) {
        loadDirectory(currentPath)
      }
    })

    return () => {
      unsubscribe()
      window.electronAPI.fm.unwatchDirectory(currentPath)
    }
  }, [currentPath, loadDirectory])

  // Navigation handlers - respect rootPath boundary
  const navigateTo = useCallback((path: string) => {
    // Prevent navigation above rootPath
    if (rootPath && !path.startsWith(rootPath)) {
      return
    }
    setCurrentPath(path)
    setSelectedFile(null)
    setRenamingFile(null)
  }, [rootPath])

  const navigateUp = useCallback(() => {
    // Don't go above rootPath
    if (rootPath && currentPath === rootPath) {
      return
    }
    const parentPath = currentPath.split('/').slice(0, -1).join('/') || '/'
    // Ensure we don't go above rootPath
    if (rootPath && parentPath.length < rootPath.length) {
      navigateTo(rootPath)
    } else {
      navigateTo(parentPath)
    }
  }, [currentPath, rootPath, navigateTo])

  const navigateToRoot = useCallback(() => {
    navigateTo(rootPath || '/')
  }, [rootPath, navigateTo])

  // File operations
  const handleFileClick = useCallback((file: FMFileEntry) => {
    setSelectedFile(file)
    onFileSelect?.(file)
  }, [onFileSelect])

  const handleFileDoubleClick = useCallback((file: FMFileEntry) => {
    if (file.isDirectory) {
      navigateTo(file.path)
    } else {
      onFileOpen?.(file)
      window.electronAPI.openFile(file.path)
    }
  }, [navigateTo, onFileOpen])

  const handleRefresh = useCallback(() => {
    loadDirectory(currentPath)
  }, [currentPath, loadDirectory])

  // Rename operations
  const startRename = useCallback((file: FMFileEntry) => {
    setRenamingFile(file)
    setRenameValue(file.name)
  }, [])

  const submitRename = useCallback(async () => {
    if (!renamingFile || !renameValue.trim() || renameValue === renamingFile.name) {
      setRenamingFile(null)
      return
    }
    
    try {
      await window.electronAPI.fm.rename(renamingFile.path, renameValue.trim())
      setRenamingFile(null)
      loadDirectory(currentPath)
    } catch (err) {
      console.error('Failed to rename:', err)
      setError(err instanceof Error ? err.message : String(err))
    }
  }, [renamingFile, renameValue, currentPath, loadDirectory])

  const cancelRename = useCallback(() => {
    setRenamingFile(null)
    setRenameValue('')
  }, [])

  // Delete operation
  const handleDelete = useCallback(async (file: FMFileEntry) => {
    try {
      await window.electronAPI.fm.delete([file.path])
      loadDirectory(currentPath)
    } catch (err) {
      console.error('Failed to delete:', err)
      setError(err instanceof Error ? err.message : String(err))
    }
  }, [currentPath, loadDirectory])

  // Copy/Cut/Paste operations
  const handleCopy = useCallback((file: FMFileEntry) => {
    setClipboard({ file, operation: 'copy' })
  }, [])

  const handleCut = useCallback((file: FMFileEntry) => {
    setClipboard({ file, operation: 'cut' })
  }, [])

  const handlePaste = useCallback(async () => {
    if (!clipboard) return
    
    try {
      const destPath = currentPath + '/' + clipboard.file.name
      
      if (clipboard.operation === 'copy') {
        await window.electronAPI.fm.copy([clipboard.file.path], destPath)
      } else {
        await window.electronAPI.fm.move([clipboard.file.path], destPath)
        setClipboard(null) // Clear clipboard after cut
      }
      
      loadDirectory(currentPath)
    } catch (err) {
      console.error('Failed to paste:', err)
      setError(err instanceof Error ? err.message : String(err))
    }
  }, [clipboard, currentPath, loadDirectory])

  // Create folder
  const startCreateFolder = useCallback(() => {
    setIsCreatingFolder(true)
    setNewFolderName(t('新文件夹'))
    setTimeout(() => {
      newFolderInputRef.current?.focus()
      newFolderInputRef.current?.select()
    }, 0)
  }, [t])

  const submitCreateFolder = useCallback(async () => {
    if (!newFolderName.trim()) {
      setIsCreatingFolder(false)
      return
    }
    
    try {
      await window.electronAPI.fm.createFolder(currentPath + '/' + newFolderName.trim())
      setIsCreatingFolder(false)
      setNewFolderName('')
      loadDirectory(currentPath)
    } catch (err) {
      console.error('Failed to create folder:', err)
      setError(err instanceof Error ? err.message : String(err))
    }
  }, [newFolderName, currentPath, loadDirectory])

  const cancelCreateFolder = useCallback(() => {
    setIsCreatingFolder(false)
    setNewFolderName('')
  }, [])

  const handleNewFolderKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      e.preventDefault()
      submitCreateFolder()
    } else if (e.key === 'Escape') {
      e.preventDefault()
      cancelCreateFolder()
    }
  }

  // Keyboard navigation
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Backspace' && !renamingFile && !isCreatingFolder) {
        e.preventDefault()
        navigateUp()
      }
      if (e.key === 'Enter' && selectedFile && !renamingFile && !isCreatingFolder) {
        e.preventDefault()
        handleFileDoubleClick(selectedFile)
      }
      if (e.key === 'F2' && selectedFile && !renamingFile && !isCreatingFolder) {
        e.preventDefault()
        startRename(selectedFile)
      }
      if (e.key === 'Delete' && selectedFile && !renamingFile && !isCreatingFolder) {
        e.preventDefault()
        handleDelete(selectedFile)
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [selectedFile, renamingFile, isCreatingFolder, navigateUp, handleFileDoubleClick, startRename, handleDelete])

  return (
    <div className={cn("flex flex-col h-full", className)}>
      {/* Toolbar */}
      <div className="flex items-center gap-2 p-2 border-b shrink-0">
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={navigateUp}
          disabled={rootPath ? currentPath === rootPath : currentPath === '/'}
          title={t('上一级')}
        >
          <ArrowUp className="h-4 w-4" />
        </Button>
        
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={navigateToRoot}
          disabled={rootPath ? currentPath === rootPath : currentPath === '/'}
          title={t('工作区根目录')}
        >
          <Home className="h-4 w-4" />
        </Button>
        
        <div className="flex-1 min-w-0">
          <Breadcrumbs path={currentPath} rootPath={rootPath} onNavigate={navigateTo} />
        </div>
        
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={startCreateFolder}
          title={t('新建文件夹')}
        >
          <FolderPlus className="h-4 w-4" />
        </Button>
        
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={handleRefresh}
          disabled={isLoading}
          title={t('刷新')}
        >
          <RefreshCw className={cn("h-4 w-4", isLoading && "animate-spin")} />
        </Button>
      </div>

      {/* File list */}
      <ContextMenu modal={true}>
        <ContextMenuTrigger asChild>
          <ScrollArea className="flex-1">
            <div className="p-2">
              {error ? (
                <div className="text-center text-destructive py-8">
                  <p className="text-sm">{error}</p>
                  <Button
                    variant="link"
                    size="sm"
                    onClick={handleRefresh}
                    className="mt-2"
                  >
                    {t('重试')}
                  </Button>
                </div>
              ) : isLoading && files.length === 0 ? (
                <div className="text-center text-muted-foreground py-8">
                  <RefreshCw className="h-6 w-6 animate-spin mx-auto mb-2" />
                  <p className="text-sm">{t('加载中...')}</p>
                </div>
              ) : files.length === 0 ? (
                <div className="text-center text-muted-foreground py-8">
                  <Folder className="h-12 w-12 mx-auto mb-2 opacity-50" />
                  <p className="text-sm">{t('文件夹为空')}</p>
                </div>
              ) : (
                <motion.div
                  className="flex flex-col gap-0.5"
                  variants={containerVariants}
                  initial="hidden"
                  animate="visible"
                >
                  {/* New folder input */}
                  <AnimatePresence>
                    {isCreatingFolder && (
                      <motion.div
                        initial={{ opacity: 0, height: 0 }}
                        animate={{ opacity: 1, height: 'auto' }}
                        exit={{ opacity: 0, height: 0 }}
                        className="flex items-center gap-3 px-3 py-2 rounded-md bg-muted/50"
                      >
                        <Folder className="h-4 w-4 text-muted-foreground shrink-0" />
                        <Input
                          ref={newFolderInputRef}
                          value={newFolderName}
                          onChange={(e) => setNewFolderName(e.target.value)}
                          onKeyDown={handleNewFolderKeyDown}
                          onBlur={submitCreateFolder}
                          className="h-6 py-0 px-1 text-sm flex-1"
                        />
                      </motion.div>
                    )}
                  </AnimatePresence>
                  
                  {/* File items */}
                  {files.map((file) => (
                    <FileItem
                      key={file.path}
                      file={file}
                      isSelected={selectedFile?.path === file.path}
                      isRenaming={renamingFile?.path === file.path}
                      renameValue={renameValue}
                      onRenameValueChange={setRenameValue}
                      onRenameSubmit={submitRename}
                      onRenameCancel={cancelRename}
                      onClick={() => handleFileClick(file)}
                      onDoubleClick={() => handleFileDoubleClick(file)}
                      onRename={() => startRename(file)}
                      onDelete={() => handleDelete(file)}
                      onCopy={() => handleCopy(file)}
                      onCut={() => handleCut(file)}
                      canPaste={!!clipboard}
                      onPaste={handlePaste}
                    />
                  ))}
                </motion.div>
              )}
            </div>
          </ScrollArea>
        </ContextMenuTrigger>
        <StyledContextMenuContent>
          <ContextMenuProvider>
            <EmptyContextMenu
              onCreateFolder={startCreateFolder}
              canPaste={!!clipboard}
              onPaste={handlePaste}
            />
          </ContextMenuProvider>
        </StyledContextMenuContent>
      </ContextMenu>

      {/* Status bar */}
      <div className="flex items-center justify-between px-3 py-1.5 border-t text-xs text-muted-foreground shrink-0">
        <span>
          {files.length} {files.length === 1 ? t('项') : t('项')}
        </span>
        {clipboard && (
          <span>
            {clipboard.operation === 'copy' ? t('已复制') : t('已剪切')}: {clipboard.file.name}
          </span>
        )}
      </div>
    </div>
  )
}
