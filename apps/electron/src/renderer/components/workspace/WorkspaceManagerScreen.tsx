/**
 * WorkspaceManagerScreen - Full-screen overlay for managing all workspaces
 *
 * Features:
 * - Card grid showing all workspaces
 * - Edit button to open workspace settings
 * - Open button to open workspace in new window
 */

import * as React from 'react'
import { useState, useEffect } from 'react'
import { motion, AnimatePresence } from 'motion/react'
import { X, ExternalLink, Settings, Folder, Plus } from 'lucide-react'
import { cn } from '@/lib/utils'
import { CrossfadeAvatar } from '@/components/ui/avatar'
import { ScrollArea } from '@/components/ui/scroll-area'
import { useT } from '@/context/LocaleContext'
import { WorkspaceCreationScreen } from './WorkspaceCreationScreen'
import type { Workspace } from '../../../shared/types'

interface WorkspaceManagerScreenProps {
  isOpen: boolean
  workspaces: Workspace[]
  activeWorkspaceId: string | null
  onClose: () => void
  onOpenWorkspace: (workspaceId: string) => void
  onEditWorkspace: (workspaceId: string) => void
  onWorkspaceCreated?: (workspace: Workspace) => void
  /** Initial marketplace app to use when creating workspace (skips app selection) */
  initialMarketplaceApp?: { id: string; name: string }
  /** Callback to clear initialMarketplaceApp after it's been used */
  onClearInitialApp?: () => void
}

interface WorkspaceCardProps {
  workspace: Workspace
  isActive: boolean
  iconDataUrl?: string
  onOpen: () => void
  onEdit: () => void
}

function WorkspaceCard({ workspace, isActive, iconDataUrl, onOpen, onEdit }: WorkspaceCardProps) {
  const t = useT()
  
  return (
    <div
      className={cn(
        'group relative flex flex-col rounded-xl bg-background shadow-minimal',
        'border border-transparent transition-all duration-200',
        'hover:shadow-md hover:border-border/50',
        isActive && 'ring-2 ring-accent/50'
      )}
    >
      {/* Card Header */}
      <div className="flex items-center gap-3 p-4 pb-2">
        <CrossfadeAvatar
          src={iconDataUrl}
          alt={workspace.name}
          className="h-10 w-10 rounded-lg ring-1 ring-border/50"
          fallbackClassName="bg-foreground/5 text-foreground/60 text-sm font-medium rounded-lg"
          fallback={workspace.name.charAt(0).toUpperCase()}
        />
        <div className="flex-1 min-w-0">
          <h3 className="font-medium text-sm truncate">{workspace.name}</h3>
          {isActive && (
            <span className="text-xs text-accent">{t('当前工作区')}</span>
          )}
        </div>
      </div>

      {/* Card Body - App Info */}
      <div className="px-4 py-2 flex-1">
        <div className="flex items-center gap-2 text-xs text-muted-foreground">
          <Folder className="h-3 w-3" />
          <span className="truncate">{workspace.rootPath}</span>
        </div>
        {workspace.appId && (
          <div className="mt-1 text-xs text-muted-foreground/70">
            {workspace.appId.replace('app.', '')}
          </div>
        )}
      </div>

      {/* Card Footer - Actions */}
      <div className="flex items-center gap-2 p-3 pt-2 border-t border-border/30">
        <button
          onClick={onEdit}
          className={cn(
            'flex-1 inline-flex items-center justify-center gap-1.5 h-8 px-3',
            'text-xs font-medium rounded-md',
            'bg-foreground/5 hover:bg-foreground/10',
            'transition-colors duration-150'
          )}
        >
          <Settings className="h-3.5 w-3.5" />
          {t('编辑')}
        </button>
        <button
          onClick={onOpen}
          className={cn(
            'flex-1 inline-flex items-center justify-center gap-1.5 h-8 px-3',
            'text-xs font-medium rounded-md',
            'bg-accent/10 text-accent hover:bg-accent/20',
            'transition-colors duration-150'
          )}
        >
          <ExternalLink className="h-3.5 w-3.5" />
          {t('打开')}
        </button>
      </div>
    </div>
  )
}

export function WorkspaceManagerScreen({
  isOpen,
  workspaces,
  activeWorkspaceId,
  onClose,
  onOpenWorkspace,
  onEditWorkspace,
  onWorkspaceCreated,
  initialMarketplaceApp,
  onClearInitialApp,
}: WorkspaceManagerScreenProps) {
  const t = useT()
  // Cache for workspace icons
  const [iconCache, setIconCache] = useState<Record<string, string>>({})
  // Show workspace creation screen
  const [showCreationScreen, setShowCreationScreen] = useState(false)
  
  // Auto-show creation screen if initialMarketplaceApp is provided
  React.useEffect(() => {
    if (isOpen && initialMarketplaceApp) {
      setShowCreationScreen(true)
    }
  }, [isOpen, initialMarketplaceApp])

  // Load workspace icons
  useEffect(() => {
    if (!isOpen) return
    const loadIcons = async () => {
      for (const workspace of workspaces) {
        // Skip if already cached or no icon
        if (iconCache[workspace.id]) continue
        if (!workspace.iconUrl?.startsWith('file://')) continue

        const urlWithoutQuery = workspace.iconUrl.split('?')[0]
        const iconFilename = urlWithoutQuery.split('/').pop()
        if (!iconFilename) continue

        try {
          const result = await window.electronAPI.readWorkspaceImage(workspace.id, iconFilename)
          if (result) {
            let dataUrl = result
            if (iconFilename.endsWith('.svg')) {
              dataUrl = `data:image/svg+xml;base64,${btoa(result)}`
            }
            setIconCache(prev => ({ ...prev, [workspace.id]: dataUrl }))
          }
        } catch (error) {
          console.error(`Failed to load icon for workspace ${workspace.id}:`, error)
        }
      }
    }
    loadIcons()
  }, [workspaces, isOpen])

  const handleOpen = (workspaceId: string) => {
    onOpenWorkspace(workspaceId)
    onClose()
  }

  const handleEdit = (workspaceId: string) => {
    onEditWorkspace(workspaceId)
    onClose()
  }

  const handleCreateWorkspace = () => {
    setShowCreationScreen(true)
  }

  const handleWorkspaceCreated = (workspace: Workspace) => {
    setShowCreationScreen(false)
    onClearInitialApp?.()
    onWorkspaceCreated?.(workspace)
    onClose()
  }
  
  const handleCreationClose = () => {
    setShowCreationScreen(false)
    onClearInitialApp?.()
  }

  if (!isOpen) return null

  return (
    <>
    {/* Workspace Creation Screen */}
    <AnimatePresence>
      {showCreationScreen && (
        <WorkspaceCreationScreen
          onWorkspaceCreated={handleWorkspaceCreated}
          onClose={handleCreationClose}
          initialMarketplaceApp={initialMarketplaceApp}
        />
      )}
    </AnimatePresence>
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      exit={{ opacity: 0 }}
      transition={{ duration: 0.2 }}
      className="fixed inset-0 z-50 flex items-center justify-center"
    >
      {/* Backdrop */}
      <div
        className="absolute inset-0 bg-background/80 backdrop-blur-sm"
        onClick={onClose}
      />

      {/* Content */}
      <motion.div
        initial={{ scale: 0.95, opacity: 0 }}
        animate={{ scale: 1, opacity: 1 }}
        exit={{ scale: 0.95, opacity: 0 }}
        transition={{ duration: 0.2, ease: [0.4, 0, 0.2, 1] }}
        className="relative w-full max-w-4xl max-h-[80vh] m-4 bg-muted/50 rounded-2xl shadow-xl overflow-hidden"
      >
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-border/30">
          <div>
            <h2 className="text-lg font-semibold">{t('管理工作区')}</h2>
            <p className="text-sm text-muted-foreground mt-0.5">
              {t('查看和管理所有工作区')}
            </p>
          </div>
          <div className="flex items-center gap-2">
            <button
              onClick={handleCreateWorkspace}
              className={cn(
                'inline-flex items-center gap-1.5 h-8 px-3',
                'text-sm font-medium rounded-lg',
                'bg-background text-foreground hover:bg-foreground/5',
                'border border-border/50 shadow-sm',
                'transition-colors duration-150'
              )}
            >
              <Plus className="h-4 w-4" />
              {t('新建')}
            </button>
            <button
              onClick={onClose}
              className="p-2 rounded-lg hover:bg-foreground/5 transition-colors"
            >
              <X className="h-5 w-5" />
            </button>
          </div>
        </div>

        {/* Body */}
        <ScrollArea className="h-[calc(80vh-80px)]">
          <div className="p-6">
            {workspaces.length === 0 ? (
              <div className="flex flex-col items-center justify-center py-12 text-muted-foreground">
                <Folder className="h-12 w-12 mb-4 opacity-50" />
                <p>{t('暂无工作区')}</p>
              </div>
            ) : (
              <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
                {workspaces.map((workspace) => (
                  <WorkspaceCard
                    key={workspace.id}
                    workspace={workspace}
                    isActive={workspace.id === activeWorkspaceId}
                    iconDataUrl={iconCache[workspace.id]}
                    onOpen={() => handleOpen(workspace.id)}
                    onEdit={() => handleEdit(workspace.id)}
                  />
                ))}
              </div>
            )}
          </div>
        </ScrollArea>
      </motion.div>
    </motion.div>
    </>
  )
}
