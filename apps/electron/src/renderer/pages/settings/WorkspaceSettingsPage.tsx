/**
 * WorkspaceSettingsPage
 *
 * Workspace-level settings for the active workspace.
 *
 * Settings:
 * - Identity (Name, Icon)
 * - Model
 * - Permissions (Default mode, Mode cycling)
 * - Advanced (Working directory, Local MCP servers)
 */

import * as React from 'react'
import { useState, useEffect, useCallback } from 'react'
import { motion, AnimatePresence } from 'motion/react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { HeaderMenu } from '@/components/ui/HeaderMenu'
import { useAppShellContext } from '@/context/AppShellContext'
import { useT } from '@/context/LocaleContext'
import { cn } from '@/lib/utils'
import { routes } from '@/lib/navigate'
import { Spinner } from '@creator-flow/ui'
import { RenameDialog } from '@/components/ui/rename-dialog'
import type { PermissionMode, ThinkingLevel, WorkspaceSettings } from '../../../shared/types'
import { PERMISSION_MODE_CONFIG } from '@creator-flow/shared/agent/mode-types'
import { DEFAULT_THINKING_LEVEL, THINKING_LEVELS } from '@creator-flow/shared/agent/thinking-levels'
import type { DetailsPageMeta } from '@/lib/navigation-registry'

import {
  SettingsSection,
  SettingsCard,
  SettingsRow,
  SettingsToggle,
  SettingsMenuSelectRow,
} from '@/components/settings'
import { useCloudModels } from '@/hooks/useCloudModels'
import { MODELS as FALLBACK_MODELS, type ModelDefinition } from '@creator-flow/shared/config/models'

export const meta: DetailsPageMeta = {
  navigator: 'settings',
  slug: 'workspace',
}

// ============================================
// Main Component
// ============================================

export default function WorkspaceSettingsPage() {
  const t = useT()
  // Get model, onModelChange, and active workspace from context
  const appShellContext = useAppShellContext()
  const onModelChange = appShellContext.onModelChange
  const activeWorkspaceId = appShellContext.activeWorkspaceId
  const onRefreshWorkspaces = appShellContext.onRefreshWorkspaces
  const customModel = appShellContext.customModel

  // Cloud models: dynamic list from backend, with local fallback
  const { models: cloudModels, loading: modelsLoading } = useCloudModels()
  const modelOptions = React.useMemo(() => {
    if (cloudModels.length > 0) {
      return cloudModels.map((m: ModelDefinition) => ({
        value: m.id,
        label: m.name,
        description: m.description || '',
      }))
    }
    return FALLBACK_MODELS.map(m => ({
      value: m.id,
      label: m.name,
      description: m.description,
    }))
  }, [cloudModels])

  // Workspace settings state
  const [wsName, setWsName] = useState('')
  const [wsNameEditing, setWsNameEditing] = useState('')
  const [renameDialogOpen, setRenameDialogOpen] = useState(false)
  const [wsIconUrl, setWsIconUrl] = useState<string | null>(null)
  const [isUploadingIcon, setIsUploadingIcon] = useState(false)
  const [wsModel, setWsModel] = useState('claude-sonnet-4-5-20250929')
  const [wsThinkingLevel, setWsThinkingLevel] = useState<ThinkingLevel>(DEFAULT_THINKING_LEVEL)
  const [permissionMode, setPermissionMode] = useState<PermissionMode>('ask')
  const [workingDirectory, setWorkingDirectory] = useState('')
  const [localMcpEnabled, setLocalMcpEnabled] = useState(true)
  const [isLoadingWorkspace, setIsLoadingWorkspace] = useState(true)

  // Mode cycling state
  const [enabledModes, setEnabledModes] = useState<PermissionMode[]>(['safe', 'ask', 'allow-all'])
  const [modeCyclingError, setModeCyclingError] = useState<string | null>(null)

  // Load workspace settings when active workspace changes
  useEffect(() => {
    const loadWorkspaceSettings = async () => {
      if (!window.electronAPI || !activeWorkspaceId) {
        setIsLoadingWorkspace(false)
        return
      }

      setIsLoadingWorkspace(true)
      try {
        const settings = await window.electronAPI.getWorkspaceSettings(activeWorkspaceId)
        if (settings) {
          setWsName(settings.name || '')
          setWsNameEditing(settings.name || '')
          setWsModel(settings.model || 'claude-sonnet-4-5-20250929')
          setWsThinkingLevel(settings.thinkingLevel || DEFAULT_THINKING_LEVEL)
          setPermissionMode(settings.permissionMode || 'ask')
          setWorkingDirectory(settings.workingDirectory || '')
          setLocalMcpEnabled(settings.localMcpEnabled ?? true)
          // Load cyclable permission modes from workspace settings
          if (settings.cyclablePermissionModes && settings.cyclablePermissionModes.length >= 2) {
            setEnabledModes(settings.cyclablePermissionModes)
          }
        }

        // Try to load workspace icon (check common extensions)
        const ICON_EXTENSIONS = ['png', 'jpg', 'jpeg', 'svg', 'webp', 'gif']
        let iconFound = false
        for (const ext of ICON_EXTENSIONS) {
          try {
            const iconData = await window.electronAPI.readWorkspaceImage(activeWorkspaceId, `./icon.${ext}`)
            // IPC returns null for missing files - continue to next extension
            if (!iconData) {
              continue
            }
            // For SVG, wrap in data URL
            if (ext === 'svg' && !iconData.startsWith('data:')) {
              setWsIconUrl(`data:image/svg+xml;base64,${btoa(iconData)}`)
            } else {
              setWsIconUrl(iconData)
            }
            iconFound = true
            break
          } catch {
            // Icon not found with this extension, try next
          }
        }
        if (!iconFound) {
          setWsIconUrl(null)
        }
      } catch (error) {
        console.error('Failed to load workspace settings:', error)
      } finally {
        setIsLoadingWorkspace(false)
      }
    }

    loadWorkspaceSettings()
  }, [activeWorkspaceId])

  // Save workspace setting
  const updateWorkspaceSetting = useCallback(
    async <K extends keyof WorkspaceSettings>(key: K, value: WorkspaceSettings[K]) => {
      if (!window.electronAPI || !activeWorkspaceId) return

      try {
        await window.electronAPI.updateWorkspaceSetting(activeWorkspaceId, key, value)
      } catch (error) {
        console.error(`Failed to save ${key}:`, error)
      }
    },
    [activeWorkspaceId]
  )

  // Workspace icon upload handler
  const handleIconUpload = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file || !activeWorkspaceId || !window.electronAPI) return

    // Validate file type
    const validTypes = ['image/png', 'image/jpeg', 'image/svg+xml', 'image/webp', 'image/gif']
    if (!validTypes.includes(file.type)) {
      console.error('Invalid file type:', file.type)
      return
    }

    setIsUploadingIcon(true)
    try {
      // Read file as base64
      const buffer = await file.arrayBuffer()
      const base64 = btoa(
        new Uint8Array(buffer).reduce((data, byte) => data + String.fromCharCode(byte), '')
      )

      // Determine extension from mime type
      const extMap: Record<string, string> = {
        'image/png': 'png',
        'image/jpeg': 'jpg',
        'image/svg+xml': 'svg',
        'image/webp': 'webp',
        'image/gif': 'gif',
      }
      const ext = extMap[file.type] || 'png'

      // Upload to workspace
      await window.electronAPI.writeWorkspaceImage(activeWorkspaceId, `./icon.${ext}`, base64, file.type)

      // Reload the icon locally for settings display
      const iconData = await window.electronAPI.readWorkspaceImage(activeWorkspaceId, `./icon.${ext}`)
      if (iconData) {
        if (ext === 'svg' && !iconData.startsWith('data:')) {
          setWsIconUrl(`data:image/svg+xml;base64,${btoa(iconData)}`)
        } else {
          setWsIconUrl(iconData)
        }
      }

      // Refresh workspaces to update sidebar icon
      onRefreshWorkspaces?.()
    } catch (error) {
      console.error('Failed to upload icon:', error)
    } finally {
      setIsUploadingIcon(false)
      // Reset the input so the same file can be selected again
      e.target.value = ''
    }
  }, [activeWorkspaceId, onRefreshWorkspaces])

  // Workspace settings handlers
  const handleModelChange = useCallback(
    async (newModel: string) => {
      setWsModel(newModel)
      await updateWorkspaceSetting('model', newModel)
      // Also update the global model context so it takes effect immediately
      onModelChange?.(newModel)
    },
    [updateWorkspaceSetting, onModelChange]
  )

  const handleThinkingLevelChange = useCallback(
    async (newLevel: ThinkingLevel) => {
      setWsThinkingLevel(newLevel)
      await updateWorkspaceSetting('thinkingLevel', newLevel)
    },
    [updateWorkspaceSetting]
  )

  const handlePermissionModeChange = useCallback(
    async (newMode: PermissionMode) => {
      setPermissionMode(newMode)
      await updateWorkspaceSetting('permissionMode', newMode)
    },
    [updateWorkspaceSetting]
  )

  const handleChangeWorkingDirectory = useCallback(async () => {
    if (!window.electronAPI) return

    try {
      const selectedPath = await window.electronAPI.openFolderDialog()
      if (selectedPath) {
        setWorkingDirectory(selectedPath)
        await updateWorkspaceSetting('workingDirectory', selectedPath)
      }
    } catch (error) {
      console.error('Failed to change working directory:', error)
    }
  }, [updateWorkspaceSetting])

  const handleClearWorkingDirectory = useCallback(async () => {
    if (!window.electronAPI) return

    try {
      setWorkingDirectory('')
      await updateWorkspaceSetting('workingDirectory', undefined)
    } catch (error) {
      console.error('Failed to clear working directory:', error)
    }
  }, [updateWorkspaceSetting])

  const handleLocalMcpEnabledChange = useCallback(
    async (enabled: boolean) => {
      setLocalMcpEnabled(enabled)
      await updateWorkspaceSetting('localMcpEnabled', enabled)
    },
    [updateWorkspaceSetting]
  )

  const handleModeToggle = useCallback(
    async (mode: PermissionMode, checked: boolean) => {
      if (!window.electronAPI) return

      // Calculate what the new modes would be
      const newModes = checked
        ? [...enabledModes, mode]
        : enabledModes.filter((m) => m !== mode)

      // Validate: at least 2 modes required
      if (newModes.length < 2) {
        setModeCyclingError('At least 2 modes required')
        // Auto-dismiss after 2 seconds
        setTimeout(() => {
          setModeCyclingError(null)
        }, 2000)
        return
      }

      // Update state and persist
      setEnabledModes(newModes)
      setModeCyclingError(null)
      try {
        await updateWorkspaceSetting('cyclablePermissionModes', newModes)
      } catch (error) {
        console.error('Failed to save mode cycling settings:', error)
      }
    },
    [enabledModes, updateWorkspaceSetting]
  )

  // Show empty state if no workspace is active
  if (!activeWorkspaceId) {
    return (
      <div className="h-full flex flex-col">
        <PanelHeader title={t('工作区设置')} actions={<HeaderMenu route={routes.view.settings('workspace')} helpFeature="workspaces" />} />
        <div className="flex-1 flex items-center justify-center">
          <p className="text-sm text-muted-foreground">{t('未选择工作区')}</p>
        </div>
      </div>
    )
  }

  // Show loading state
  if (isLoadingWorkspace) {
    return (
      <div className="h-full flex flex-col">
        <PanelHeader title={t('工作区设置')} actions={<HeaderMenu route={routes.view.settings('workspace')} helpFeature="workspaces" />} />
        <div className="flex-1 flex items-center justify-center">
          <Spinner className="text-muted-foreground" />
        </div>
      </div>
    )
  }

  return (
    <div className="h-full flex flex-col">
      <PanelHeader title={t('工作区设置')} actions={<HeaderMenu route={routes.view.settings('workspace')} helpFeature="workspaces" />} />
      <div className="flex-1 min-h-0 mask-fade-y">
        <ScrollArea className="h-full">
          <div className="px-5 py-7 max-w-3xl mx-auto">
          <div className="space-y-8">
            {/* Workspace Info */}
            <SettingsSection title={t('工作区信息')}>
              <SettingsCard>
                <SettingsRow
                  label={t('名称')}
                  description={wsName || t('未命名')}
                  action={
                    <button
                      type="button"
                      onClick={() => {
                        setWsNameEditing(wsName)
                        setRenameDialogOpen(true)
                      }}
                      className="inline-flex items-center h-8 px-3 text-sm rounded-lg bg-background shadow-minimal hover:bg-foreground/[0.02] transition-colors"
                    >
                      {t('编辑')}
                    </button>
                  }
                />
                <SettingsRow
                  label={t('图标')}
                  action={
                    <label className="cursor-pointer">
                      <input
                        type="file"
                        accept="image/png,image/jpeg,image/svg+xml,image/webp,image/gif"
                        onChange={handleIconUpload}
                        className="sr-only"
                        disabled={isUploadingIcon}
                      />
                      <span className="inline-flex items-center h-8 px-3 text-sm rounded-lg bg-background shadow-minimal hover:bg-foreground/[0.02] transition-colors">
                        {isUploadingIcon ? t('上传中...') : t('更改')}
                      </span>
                    </label>
                  }
                >
                  <div
                    className={cn(
                      'w-6 h-6 rounded-full overflow-hidden bg-foreground/5 flex items-center justify-center',
                      'ring-1 ring-border/50'
                    )}
                  >
                    {isUploadingIcon ? (
                      <Spinner className="text-muted-foreground text-[8px]" />
                    ) : wsIconUrl ? (
                      <img src={wsIconUrl} alt="" className="w-full h-full object-cover" />
                    ) : (
                      <span className="text-xs font-medium text-muted-foreground">
                        {wsName?.charAt(0)?.toUpperCase() || 'W'}
                      </span>
                    )}
                  </div>
                </SettingsRow>
              </SettingsCard>

              <RenameDialog
                open={renameDialogOpen}
                onOpenChange={setRenameDialogOpen}
                title={t('重命名工作区')}
                value={wsNameEditing}
                onValueChange={setWsNameEditing}
                onSubmit={() => {
                  const newName = wsNameEditing.trim()
                  if (newName && newName !== wsName) {
                    setWsName(newName)
                    updateWorkspaceSetting('name', newName)
                    onRefreshWorkspaces?.()
                  }
                  setRenameDialogOpen(false)
                }}
                placeholder={t('输入工作区名称...')}
              />
            </SettingsSection>

            {/* Model */}
            <SettingsSection title={t('模型')}>
              <SettingsCard>
                {/* When a custom API connection is active, model is fixed — show info instead of selector */}
                {customModel ? (
                  <SettingsRow
                    label={t('默认模型')}
                    description={t('通过 API 连接设置')}
                  >
                    <span className="text-sm text-muted-foreground">{customModel}</span>
                  </SettingsRow>
                ) : (
                  <SettingsMenuSelectRow
                    label={t('默认模型')}
                    description={t('新聊天的 AI 模型')}
                    value={wsModel}
                    onValueChange={handleModelChange}
                    options={modelOptions}
                  />
                )}
                <SettingsMenuSelectRow
                  label={t('思考级别')}
                  description={t('新聊天的推理深度')}
                  value={wsThinkingLevel}
                  onValueChange={(v) => handleThinkingLevelChange(v as ThinkingLevel)}
                  options={THINKING_LEVELS.map(({ id, name, description }) => ({
                    value: id,
                    label: name,
                    description,
                  }))}
                />
              </SettingsCard>
            </SettingsSection>

            {/* Permissions */}
            <SettingsSection title={t('权限')}>
              <SettingsCard>
                <SettingsMenuSelectRow
                  label={t('默认模式')}
                  description={t('控制 AI 可以做什么')}
                  value={permissionMode}
                  onValueChange={(v) => handlePermissionModeChange(v as PermissionMode)}
                  options={[
                    { value: 'safe', label: PERMISSION_MODE_CONFIG['safe'].shortName, description: t('只读，不允许更改') },
                    { value: 'ask', label: PERMISSION_MODE_CONFIG['ask'].shortName, description: t('编辑前提示') },
                    { value: 'allow-all', label: PERMISSION_MODE_CONFIG['allow-all'].shortName, description: t('完全自主执行') },
                  ]}
                />
              </SettingsCard>
            </SettingsSection>

            {/* Mode Cycling */}
            <SettingsSection
              title={t('模式循环')}
              description={t('选择使用 Shift+Tab 循环的模式')}
            >
              <SettingsCard>
                {(['safe', 'ask', 'allow-all'] as const).map((m) => {
                  const config = PERMISSION_MODE_CONFIG[m]
                  const isEnabled = enabledModes.includes(m)
                  return (
                    <SettingsToggle
                      key={m}
                      label={config.displayName}
                      description={config.description}
                      checked={isEnabled}
                      onCheckedChange={(checked) => handleModeToggle(m, checked)}
                    />
                  )
                })}
              </SettingsCard>
              <AnimatePresence>
                {modeCyclingError && (
                  <motion.p
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    transition={{ duration: 0.2, ease: [0.4, 0, 0.2, 1] }}
                    className="text-xs text-destructive mt-1 overflow-hidden"
                  >
                    {modeCyclingError}
                  </motion.p>
                )}
              </AnimatePresence>
            </SettingsSection>

            {/* Advanced */}
            <SettingsSection title={t('高级')}>
              <SettingsCard>
                <SettingsRow
                  label={t('默认工作目录')}
                  description={workingDirectory || t('未设置（使用会话文件夹）')}
                  action={
                    <div className="flex items-center gap-2">
                      {workingDirectory && (
                        <button
                          type="button"
                          onClick={handleClearWorkingDirectory}
                          className="inline-flex items-center h-8 px-3 text-sm rounded-lg bg-background shadow-minimal hover:bg-foreground/[0.02] transition-colors text-foreground/60 hover:text-foreground"
                        >
                          {t('清除')}
                        </button>
                      )}
                      <button
                        type="button"
                        onClick={handleChangeWorkingDirectory}
                        className="inline-flex items-center h-8 px-3 text-sm rounded-lg bg-background shadow-minimal hover:bg-foreground/[0.02] transition-colors"
                      >
                        {t('更改...')}
                      </button>
                    </div>
                  }
                />
                <SettingsToggle
                  label={t('本地 MCP 服务器')}
                  description={t('启用 stdio 子进程服务器')}
                  checked={localMcpEnabled}
                  onCheckedChange={handleLocalMcpEnabledChange}
                />
              </SettingsCard>
            </SettingsSection>

          </div>
        </div>
        </ScrollArea>
      </div>
    </div>
  )
}
