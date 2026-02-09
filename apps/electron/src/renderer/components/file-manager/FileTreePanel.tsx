/**
 * FileTreePanel - Tree view component for file system navigation
 *
 * Features:
 * - Recursive tree structure showing folders and files
 * - Single click to expand/collapse folders
 * - Single click on file to select and preview
 * - Real-time directory watching
 * - Keyboard navigation support
 */

import * as React from 'react'
import { useState, useEffect, useCallback, useMemo } from 'react'
import { AnimatePresence, motion } from 'motion/react'
import {
  ChevronRight,
  ChevronDown,
  File,
  Folder,
  FolderOpen,
  FileText,
  FileCode,
  Image,
  RefreshCw,
  FolderSearch,
  Sparkles,
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import {
  ContextMenu,
  ContextMenuTrigger,
  StyledContextMenuContent,
  StyledContextMenuItem,
  StyledContextMenuSeparator,
} from '@/components/ui/styled-context-menu'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import type { FMFileEntry } from '../../../shared/types'

export interface FileTreePanelProps {
  /** Root path to display (workspace root) */
  rootPath: string
  /** Currently selected file path */
  selectedPath?: string
  /** Called when a file is selected */
  onSelectFile?: (file: FMFileEntry) => void
  /** Called when a directory is selected (for context, not preview) */
  onSelectDirectory?: (path: string) => void
  className?: string
  /** 紧凑模式：隐藏 header，用于嵌入场景 */
  compact?: boolean
  /** Increment to trigger a refresh from outside */
  refreshTrigger?: number
}

interface TreeNode {
  entry: FMFileEntry
  children?: TreeNode[]
  isLoading?: boolean
  isExpanded?: boolean
}

/**
 * Get icon for file based on extension
 */
function getFileIcon(file: FMFileEntry, isExpanded?: boolean) {
  const iconClass = "h-4 w-4 shrink-0 text-muted-foreground"

  if (file.isDirectory) {
    return isExpanded 
      ? <FolderOpen className={iconClass} />
      : <Folder className={iconClass} />
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
 * Single tree item component
 */
function TreeItem({
  node,
  depth,
  selectedPath,
  expandedPaths,
  onToggleExpand,
  onSelectFile,
  onSelectDirectory,
}: {
  node: TreeNode
  depth: number
  selectedPath?: string
  expandedPaths: Set<string>
  onToggleExpand: (path: string) => void
  onSelectFile?: (file: FMFileEntry) => void
  onSelectDirectory?: (path: string) => void
}) {
  const t = useT()
  const isExpanded = expandedPaths.has(node.entry.path)
  const isSelected = selectedPath === node.entry.path
  const isDirectory = node.entry.isDirectory

  const handleClick = useCallback(() => {
    if (isDirectory) {
      onToggleExpand(node.entry.path)
      onSelectDirectory?.(node.entry.path)
    } else {
      onSelectFile?.(node.entry)
    }
  }, [isDirectory, node.entry, onToggleExpand, onSelectFile, onSelectDirectory])

  return (
    <div>
      <ContextMenu modal={true}>
        <ContextMenuTrigger asChild>
          <button
            onClick={handleClick}
            className={cn(
              "w-full flex items-center gap-1 py-1 px-2 text-sm text-left rounded-md transition-colors",
              "hover:bg-muted/50",
              isSelected && "bg-muted text-foreground"
            )}
            style={{ paddingLeft: `${depth * 16 + 8}px` }}
          >
            {/* Expand/collapse indicator for directories */}
            {isDirectory ? (
              <span className="w-4 h-4 flex items-center justify-center shrink-0">
                {node.isLoading ? (
                  <RefreshCw className="h-3 w-3 animate-spin text-muted-foreground" />
                ) : isExpanded ? (
                  <ChevronDown className="h-3 w-3 text-muted-foreground" />
                ) : (
                  <ChevronRight className="h-3 w-3 text-muted-foreground" />
                )}
              </span>
            ) : (
              <span className="w-4 h-4 shrink-0" />
            )}

            {/* File/folder icon */}
            {getFileIcon(node.entry, isExpanded)}

            {/* Name */}
            <span className="truncate flex-1">{node.entry.name}</span>
          </button>
        </ContextMenuTrigger>
        <StyledContextMenuContent>
          <StyledContextMenuItem onClick={() => {
            const ref = isDirectory ? `[folder:${node.entry.path}]` : `[file:${node.entry.path}]`
            window.dispatchEvent(new CustomEvent('craft:insert-text', { detail: { text: ref + ' ' } }))
          }}>
            <Sparkles />
            <span className="flex-1">{t('AI 编辑')}</span>
          </StyledContextMenuItem>
          <StyledContextMenuSeparator />
          <StyledContextMenuItem onClick={() => window.electronAPI.showInFolder(node.entry.path)}>
            <FolderSearch />
            <span className="flex-1">{t('在访达中显示')}</span>
          </StyledContextMenuItem>
        </StyledContextMenuContent>
      </ContextMenu>

      {/* Children (if expanded) */}
      <AnimatePresence initial={false}>
        {isDirectory && isExpanded && node.children && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.15, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            {node.children.map((child) => (
              <TreeItem
                key={child.entry.path}
                node={child}
                depth={depth + 1}
                selectedPath={selectedPath}
                expandedPaths={expandedPaths}
                onToggleExpand={onToggleExpand}
                onSelectFile={onSelectFile}
                onSelectDirectory={onSelectDirectory}
              />
            ))}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  )
}

/**
 * FileTreePanel component
 */
export function FileTreePanel({
  rootPath,
  selectedPath,
  onSelectFile,
  onSelectDirectory,
  className,
  compact,
  refreshTrigger,
}: FileTreePanelProps) {
  const t = useT()
  const [rootNode, setRootNode] = useState<TreeNode | null>(null)
  const [expandedPaths, setExpandedPaths] = useState<Set<string>>(new Set([rootPath]))
  const [loadingPaths, setLoadingPaths] = useState<Set<string>>(new Set())
  const [error, setError] = useState<string | null>(null)
  const [isRefreshing, setIsRefreshing] = useState(false)

  // Load directory contents
  const loadDirectory = useCallback(async (path: string): Promise<FMFileEntry[]> => {
    try {
      const entries = await window.electronAPI.fm.listDirectory(path)
      // Sort: directories first, then files, both alphabetically
      return entries.sort((a, b) => {
        if (a.isDirectory !== b.isDirectory) {
          return a.isDirectory ? -1 : 1
        }
        return a.name.localeCompare(b.name)
      })
    } catch (err) {
      console.error('Failed to load directory:', path, err)
      throw err
    }
  }, [])

  // Build tree node with children
  const buildTreeNode = useCallback(async (entry: FMFileEntry, loadChildren: boolean = false): Promise<TreeNode> => {
    const node: TreeNode = { entry }
    
    if (entry.isDirectory && loadChildren) {
      try {
        const children = await loadDirectory(entry.path)
        node.children = await Promise.all(
          children.map(child => buildTreeNode(child, false))
        )
      } catch {
        node.children = []
      }
    }
    
    return node
  }, [loadDirectory])

  // Load root directory on mount
  useEffect(() => {
    const loadRoot = async () => {
      try {
        setError(null)
        const entries = await loadDirectory(rootPath)
        const rootEntry: FMFileEntry = {
          path: rootPath,
          name: rootPath.split('/').pop() || t('工作区'),
          isDirectory: true,
        }
        
        const children = await Promise.all(
          entries.map(entry => buildTreeNode(entry, false))
        )
        
        setRootNode({
          entry: rootEntry,
          children,
          isExpanded: true,
        })
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      }
    }
    
    loadRoot()
  }, [rootPath, loadDirectory, buildTreeNode, t])

  // Toggle expand/collapse
  const handleToggleExpand = useCallback(async (path: string) => {
    const isExpanded = expandedPaths.has(path)
    
    if (isExpanded) {
      // Collapse
      setExpandedPaths(prev => {
        const next = new Set(prev)
        next.delete(path)
        return next
      })
    } else {
      // Expand - load children if needed
      setLoadingPaths(prev => new Set(prev).add(path))
      
      try {
        const entries = await loadDirectory(path)
        
        // Update the tree with loaded children
        setRootNode(prev => {
          if (!prev) return prev
          
          const updateNode = (node: TreeNode): TreeNode => {
            if (node.entry.path === path) {
              return {
                ...node,
                children: entries.map(entry => ({
                  entry,
                  children: entry.isDirectory ? undefined : undefined,
                })),
                isExpanded: true,
              }
            }
            
            if (node.children) {
              return {
                ...node,
                children: node.children.map(updateNode),
              }
            }
            
            return node
          }
          
          return updateNode(prev)
        })
        
        setExpandedPaths(prev => new Set(prev).add(path))
      } catch (err) {
        console.error('Failed to expand directory:', err)
      } finally {
        setLoadingPaths(prev => {
          const next = new Set(prev)
          next.delete(path)
          return next
        })
      }
    }
  }, [expandedPaths, loadDirectory])

  // Refresh tree
  const handleRefresh = useCallback(async () => {
    if (!rootPath || isRefreshing) return

    setIsRefreshing(true)
    try {
      setError(null)
      const entries = await loadDirectory(rootPath)
      const rootEntry: FMFileEntry = {
        path: rootPath,
        name: rootPath.split('/').pop() || t('工作区'),
        isDirectory: true,
      }

      const children = entries.map(entry => ({
        entry,
        children: entry.isDirectory ? undefined : undefined,
      }))

      setRootNode({
        entry: rootEntry,
        children,
        isExpanded: true,
      })

      // Keep root expanded
      setExpandedPaths(new Set([rootPath]))
      toast.success(t('刷新完成'))
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err))
      toast.error(t('刷新失败'))
    } finally {
      setIsRefreshing(false)
    }
  }, [rootPath, isRefreshing, loadDirectory, t])

  // External refresh trigger
  useEffect(() => {
    if (refreshTrigger && refreshTrigger > 0) {
      handleRefresh()
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [refreshTrigger])

  // Merge loading state into tree
  const treeWithLoadingState = useMemo(() => {
    if (!rootNode) return null
    
    const addLoadingState = (node: TreeNode): TreeNode => ({
      ...node,
      isLoading: loadingPaths.has(node.entry.path),
      children: node.children?.map(addLoadingState),
    })
    
    return addLoadingState(rootNode)
  }, [rootNode, loadingPaths])

  if (error) {
    return (
      <div className={cn("flex flex-col h-full", className)}>
        <div className="flex items-center justify-center h-full text-center p-4">
          <div>
            <p className="text-sm text-destructive mb-2">{error}</p>
            <Button variant="outline" size="sm" onClick={handleRefresh}>
              {t('重试')}
            </Button>
          </div>
        </div>
      </div>
    )
  }

  if (!treeWithLoadingState) {
    return (
      <div className={cn("flex flex-col h-full", className)}>
        <div className="flex items-center justify-center h-full">
          <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
        </div>
      </div>
    )
  }

  return (
    <div className={cn("flex flex-col h-full", className)}>
      {/* Header */}
      {!compact && (
      <div className="flex items-center justify-between px-3 py-2 border-b shrink-0 relative z-panel">
        <span className="text-sm font-medium truncate" title={rootPath}>{rootPath}</span>
        <Button
          variant="ghost"
          size="icon"
          className="h-6 w-6 titlebar-no-drag"
          onClick={handleRefresh}
          disabled={isRefreshing}
          title={t('刷新')}
        >
          <RefreshCw className={cn("h-3.5 w-3.5", isRefreshing && "animate-spin")} />
        </Button>
      </div>
      )}
      
      {/* Tree */}
      <ScrollArea className="flex-1">
        <div className="py-1">
          {treeWithLoadingState.children?.map((child) => (
            <TreeItem
              key={child.entry.path}
              node={child}
              depth={0}
              selectedPath={selectedPath}
              expandedPaths={expandedPaths}
              onToggleExpand={handleToggleExpand}
              onSelectFile={onSelectFile}
              onSelectDirectory={onSelectDirectory}
            />
          ))}
          
          {treeWithLoadingState.children?.length === 0 && (
            <div className="text-center text-muted-foreground py-8 text-sm">
              {t('文件夹为空')}
            </div>
          )}
        </div>
      </ScrollArea>
    </div>
  )
}
