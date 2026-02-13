/**
 * SourceInfoPage
 *
 * Displays source details including connection info, authentication status,
 * documentation (guide.md), and metadata. View-only.
 */

import * as React from 'react'
import { useEffect, useState, useMemo, useCallback } from 'react'
import { AlertCircle, Loader2, Plug } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { EditPopover, EditButton, getEditConfig } from '@/components/ui/EditPopover'
import { SourceAvatar } from '@/components/ui/source-avatar'
import { SourceMenu } from '@/components/app-shell/SourceMenu'
import { cn } from '@/lib/utils'
import { routes, navigate } from '@/lib/navigate'
import { useNavigation } from '@/contexts/NavigationContext'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import {
  Info_Page,
  Info_Section,
  Info_Table,
  Info_Alert,
  Info_Markdown,
  PermissionsDataTable,
  ToolsDataTable,
  type PermissionRow,
  type ToolRow,
} from '@/components/info'
import type { LoadedSource, McpToolWithPermission } from '../../shared/types'
import type { PermissionsConfigFile } from '@sprouty-ai/shared/agent/modes'

interface SourceInfoPageProps {
  sourceSlug: string
  workspaceId: string
  /** Optional callback when source is deleted */
  onDelete?: () => void
}

/**
 * Format timestamp to relative time
 */
function formatRelativeTime(timestamp: number | undefined, t: (key: string) => string): string {
  if (!timestamp) return t('从未')

  const now = Date.now()
  const diff = now - timestamp
  const minutes = Math.floor(diff / 60000)
  const hours = Math.floor(diff / 3600000)
  const days = Math.floor(diff / 86400000)

  if (minutes < 1) return t('刚刚')
  if (minutes < 60) return t('{count}分钟前').replace('{count}', String(minutes))
  if (hours < 24) return t('{count}小时前').replace('{count}', String(hours))
  return t('{count}天前').replace('{count}', String(days))
}

/**
 * Get source URL for display
 */
function getSourceUrl(source: LoadedSource): string | null {
  const { type, mcp, api, local } = source.config

  if (type === 'mcp' && mcp?.url) return mcp.url
  if (type === 'api' && api?.baseUrl) return api.baseUrl
  if (type === 'local' && local?.path) return local.path

  return null
}

/**
 * Convert permissions config to PermissionRow[] for API/local sources
 */
function buildApiPermissionsData(config: PermissionsConfigFile): PermissionRow[] {
  const rows: PermissionRow[] = []

  // Blocked Tools
  config.blockedTools?.forEach((item) => {
    const pattern = typeof item === 'string' ? item : item.pattern
    const comment = typeof item === 'string' ? null : item.comment
    rows.push({ access: 'blocked', type: 'tool', pattern, comment })
  })

  // Allowed Bash Patterns
  config.allowedBashPatterns?.forEach((item) => {
    const pattern = typeof item === 'string' ? item : item.pattern
    const comment = typeof item === 'string' ? null : item.comment
    rows.push({ access: 'allowed', type: 'bash', pattern, comment })
  })

  // Allowed API Endpoints
  config.allowedApiEndpoints?.forEach((item) => {
    const pattern = `${item.method} ${item.path}`
    const comment = typeof item === 'object' && 'comment' in item ? item.comment : null
    rows.push({ access: 'allowed', type: 'api', pattern, comment })
  })

  return rows
}

/**
 * Convert permissions config to PermissionRow[] for MCP sources
 */
function buildMcpPermissionsData(config: PermissionsConfigFile): PermissionRow[] {
  const rows: PermissionRow[] = []

  // Blocked Tools
  config.blockedTools?.forEach((item) => {
    const pattern = typeof item === 'string' ? item : item.pattern
    const comment = typeof item === 'string' ? null : item.comment
    rows.push({ access: 'blocked', type: 'mcp', pattern, comment })
  })

  // Allowed MCP Patterns
  config.allowedMcpPatterns?.forEach((item) => {
    const pattern = typeof item === 'string' ? item : item.pattern
    const comment = typeof item === 'string' ? null : item.comment
    rows.push({ access: 'allowed', type: 'mcp', pattern, comment })
  })

  return rows
}

/**
 * Convert MCP tools to ToolRow[]
 */
function buildToolsData(tools: McpToolWithPermission[]): ToolRow[] {
  return tools.map((tool) => ({
    name: tool.name,
    description: tool.description || '',
    permission: tool.allowed ? 'allowed' : 'requires-permission',
  }))
}

/**
 * Get contextual description for Connection section based on source type
 */
function getConnectionDescription(source: LoadedSource, t: (key: string) => string): string {
  const { type, mcp } = source.config

  if (type === 'mcp') {
    if (mcp?.transport === 'stdio') {
      return t('启动此 MCP 服务器的本地命令')
    }
    return t('服务器 URL 和连接状态')
  }
  if (type === 'api') {
    return t('API 请求的基础 URL')
  }
  if (type === 'local') {
    return t('此数据源的文件系统路径')
  }
  return t('连接详情')
}

/**
 * Get contextual description for Permissions section based on source type
 */
function getPermissionsDescription(source: LoadedSource, t: (key: string) => string): string {
  const { type } = source.config

  if (type === 'mcp') {
    return t('探索模式下允许的工具模式')
  }
  if (type === 'api') {
    return t('探索模式下允许的 API 端点')
  }
  return t('探索模式的访问规则')
}

export default function SourceInfoPage({ sourceSlug, workspaceId, onDelete }: SourceInfoPageProps) {
  const t = useT()
  const { navigateToSource } = useNavigation()
  const [source, setSource] = useState<LoadedSource | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [permissionsConfig, setPermissionsConfig] = useState<PermissionsConfigFile | null>(null)
  const [mcpTools, setMcpTools] = useState<McpToolWithPermission[] | null>(null)
  const [mcpToolsLoading, setMcpToolsLoading] = useState(false)
  const [mcpToolsError, setMcpToolsError] = useState<string | null>(null)
  const [localMcpEnabled, setLocalMcpEnabled] = useState(true)
  const [testingConnection, setTestingConnection] = useState(false)


  // Load source data
  useEffect(() => {
    let isMounted = true
    setLoading(true)
    setError(null)

    const loadSource = async () => {
      try {
        const sources = await window.electronAPI.getSources(workspaceId)

        if (!isMounted) return

        const found = sources.find((s) => s.config.slug === sourceSlug)
        if (found) {
          setSource(found)

          const config = await window.electronAPI.getSourcePermissionsConfig(workspaceId, sourceSlug)
          if (isMounted) {
            setPermissionsConfig(config)
          }
        } else {
          setError('Source not found')
        }
      } catch (err) {
        if (!isMounted) return
        setError(err instanceof Error ? err.message : 'Failed to load source')
      } finally {
        if (isMounted) setLoading(false)
      }
    }

    loadSource()

    return () => {
      isMounted = false
    }
  }, [workspaceId, sourceSlug])

  // Reusable function to load MCP tools
  const loadMcpTools = useCallback(async () => {
    setMcpToolsLoading(true)
    setMcpToolsError(null)
    try {
      const result = await window.electronAPI.getMcpTools(workspaceId, sourceSlug)
      if (result.success && result.tools) {
        setMcpTools(result.tools)
      } else {
        setMcpToolsError(result.error || 'Failed to load tools')
      }
    } catch (err) {
      setMcpToolsError(err instanceof Error ? err.message : 'Failed to load tools')
    } finally {
      setMcpToolsLoading(false)
    }
  }, [workspaceId, sourceSlug])

  // Reusable function to reload source data
  const reloadSource = useCallback(async () => {
    try {
      const sources = await window.electronAPI.getSources(workspaceId)
      const found = sources.find((s) => s.config.slug === sourceSlug)
      if (found) {
        setSource(found)
      }
    } catch (err) {
      console.error('[SourceInfoPage] Failed to reload source:', err)
    }
  }, [workspaceId, sourceSlug])

  // Load MCP tools when source is loaded and is MCP type
  useEffect(() => {
    if (!source || source.config.type !== 'mcp') {
      setMcpTools(null)
      setMcpToolsError(null)
      return
    }

    loadMcpTools()
  }, [source, workspaceId, sourceSlug, loadMcpTools])

  // Load workspace settings (for localMcpEnabled)
  useEffect(() => {
    if (!workspaceId) return
    window.electronAPI.getWorkspaceSettings(workspaceId).then((settings) => {
      if (settings) {
        setLocalMcpEnabled(settings.localMcpEnabled ?? true)
      }
    }).catch((err) => {
      console.error('[SourceInfoPage] Failed to load workspace settings:', err)
    })
  }, [workspaceId])

  // Listen for source folder changes
  useEffect(() => {
    if (!window.electronAPI?.onSourcesChanged) return

    const cleanup = window.electronAPI.onSourcesChanged((sources) => {
      const updated = sources.find((s) => s.config.slug === sourceSlug)

      if (updated) {
        setSource(updated)

        const loadPermissionsConfig = async () => {
          try {
            const config = await window.electronAPI.getSourcePermissionsConfig(workspaceId, sourceSlug)
            setPermissionsConfig(config)
          } catch (err) {
            console.error('[SourceInfoPage] Failed to reload permissions config:', err)
          }
        }
        loadPermissionsConfig()
      }
    })

    return cleanup
  }, [sourceSlug, workspaceId])

  // Compute source URL
  const sourceUrl = useMemo(() => source ? getSourceUrl(source) : null, [source])

  // Build data for PermissionsDataTable
  const apiPermissionsData = useMemo(() => {
    if (!permissionsConfig || source?.config.type === 'mcp') return []
    return buildApiPermissionsData(permissionsConfig)
  }, [permissionsConfig, source])

  const mcpPermissionsData = useMemo(() => {
    if (!permissionsConfig || source?.config.type !== 'mcp') return []
    return buildMcpPermissionsData(permissionsConfig)
  }, [permissionsConfig, source])

  // Build data for ToolsDataTable
  const toolsData = useMemo(() => {
    if (!mcpTools) return []
    return buildToolsData(mcpTools)
  }, [mcpTools])

  // Handle opening URL (website or folder)
  const handleOpenUrl = useCallback(async () => {
    if (!source || !sourceUrl) return
    if (window.electronAPI) {
      if (sourceUrl.startsWith('http://') || sourceUrl.startsWith('https://')) {
        await window.electronAPI.openUrl(sourceUrl)
      } else {
        await window.electronAPI.showInFolder(sourceUrl)
      }
    }
  }, [source, sourceUrl])

  // Handle opening source folder
  const handleOpenSourceFolder = useCallback(async () => {
    if (!source) return
    if (window.electronAPI) {
      await window.electronAPI.showInFolder(source.folderPath)
    }
  }, [source])

  // Handle editing guide.md - opens in system default text editor
  const handleEditGuide = useCallback(async () => {
    if (!source) return

    const guidePath = `${source.folderPath}/guide.md`
    await window.electronAPI.openFile(guidePath)
  }, [source])

  // Handle editing config.json - opens in system default text editor
  const handleEditConfig = useCallback(async () => {
    if (!source) return

    const configPath = `${source.folderPath}/config.json`
    await window.electronAPI.openFile(configPath)
  }, [source])

  // Handle editing permissions.json - opens in system default text editor
  const handleEditPermissions = useCallback(async () => {
    if (!source) return

    const permissionsPath = `${source.folderPath}/permissions.json`
    await window.electronAPI.openFile(permissionsPath)
  }, [source])

  // Handle deleting source (navigates to source list, preserving current filter)
  const handleDelete = useCallback(async () => {
    if (!source) return
    try {
      await window.electronAPI.deleteSource(workspaceId, sourceSlug)
      toast.success(t('已删除数据源') + `: ${source.config.name}`)
      navigateToSource() // Navigate to source list, preserving filter
      onDelete?.()
    } catch (err) {
      toast.error(t('删除数据源失败'), {
        description: err instanceof Error ? err.message : t('未知错误'),
      })
    }
  }, [source, workspaceId, sourceSlug, onDelete, navigateToSource, t])

  // Handle opening in new window
  const handleOpenInNewWindow = useCallback(() => {
    window.electronAPI.openUrl(`sproutyai://sources/source/${sourceSlug}?window=focused`)
  }, [sourceSlug])

  // Handle testing MCP connection
  const handleTestConnection = useCallback(async () => {
    if (!source || source.config.type !== 'mcp') return
    setTestingConnection(true)
    try {
      const result = await window.electronAPI.testSourceConnection(workspaceId, sourceSlug)
      if (result.success) {
        toast.success(t('连接成功'))
        // Reload source data to reflect updated connectionStatus
        await reloadSource()
        // Reload MCP tools
        await loadMcpTools()
      } else {
        toast.error(t('连接失败'), {
          description: result.error || t('未知错误'),
        })
        // Reload source to show updated status
        await reloadSource()
      }
    } catch (err) {
      toast.error(t('连接测试出错'), {
        description: err instanceof Error ? err.message : t('未知错误'),
      })
    } finally {
      setTestingConnection(false)
    }
  }, [source, workspaceId, sourceSlug, reloadSource, loadMcpTools, t])

  // Get source name for header
  const sourceName = source?.config.name || sourceSlug

  return (
    <Info_Page
      loading={loading}
      error={error ?? undefined}
      empty={!source && !loading && !error ? t('未找到数据源') : undefined}
    >
      <Info_Page.Header
        title={sourceName}
        titleMenu={
          <SourceMenu
            sourceSlug={sourceSlug}
            sourceName={sourceName}
            onOpenInNewWindow={handleOpenInNewWindow}
            onShowInFinder={handleOpenSourceFolder}
            onDelete={handleDelete}
          />
        }
      />

      {source && (
        <Info_Page.Content>
          {/* Hero: Avatar, title, and tagline */}
          <Info_Page.Hero
            avatar={<SourceAvatar source={source} fluid />}
            title={source.config.name}
            tagline={source.config.tagline}
          />

          {/* Disabled Warning */}
          {source.config.mcp?.transport === 'stdio' && !localMcpEnabled && (
            <Info_Alert variant="warning" icon={<AlertCircle className="h-4 w-4" />}>
              <Info_Alert.Title>{t('数据源已禁用')}</Info_Alert.Title>
              <Info_Alert.Description>
                {t('本地 MCP 服务器已在设置 > 高级中禁用。启用它们以使用此数据源。')}
              </Info_Alert.Description>
            </Info_Alert>
          )}

          {/* Connection */}
          <Info_Section
            title={t('连接')}
            description={getConnectionDescription(source, t)}
            actions={
              <div className="flex items-center gap-1.5">
                {source.config.type === 'mcp' && (
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={handleTestConnection}
                    disabled={testingConnection}
                  >
                    {testingConnection ? (
                      <Loader2 className="h-3.5 w-3.5 animate-spin" />
                    ) : (
                      <Plug className="h-3.5 w-3.5" />
                    )}
                    {t('测试连接')}
                  </Button>
                )}
                {/* EditPopover for AI-assisted config.json editing with "Edit File" as secondary action */}
                <EditPopover
                  trigger={<EditButton />}
                  {...getEditConfig('source-config', source.folderPath)}
                  secondaryAction={{
                    label: t('编辑文件'),
                    onClick: handleEditConfig,
                  }}
                />
              </div>
            }
          >
            <Info_Table
              footer={source.config.connectionError && (
                <div className="px-4 py-2 border-t border-border/30 bg-destructive/5">
                  <div className="flex items-start gap-2 text-sm text-destructive">
                    <AlertCircle className="h-4 w-4 shrink-0 mt-0.5" />
                    <span>{source.config.connectionError}</span>
                  </div>
                </div>
              )}
            >
              <Info_Table.Row label={t('类型')} value={source.config.type.toUpperCase()} />
              {sourceUrl && (
                <Info_Table.Row label="URL">
                  <button
                    onClick={handleOpenUrl}
                    className="truncate hover:underline text-foreground focus:outline-none focus-visible:underline text-left block w-full"
                  >
                    {sourceUrl}
                  </button>
                </Info_Table.Row>
              )}
              <Info_Table.Row label={t('上次测试')} value={formatRelativeTime(source.config.lastTestedAt, t)} />
            </Info_Table>
          </Info_Section>

          {/* Permissions - for API and local sources */}
          {source.config.type !== 'mcp' && permissionsConfig && apiPermissionsData.length > 0 && (
            <Info_Section
              title={t('权限')}
              description={getPermissionsDescription(source, t)}
              actions={
                // EditPopover for AI-assisted permissions.json editing
                <EditPopover
                  trigger={<EditButton />}
                  {...getEditConfig('source-permissions', source.folderPath)}
                  secondaryAction={{
                    label: t('编辑文件'),
                    onClick: handleEditPermissions,
                  }}
                />
              }
            >
              <PermissionsDataTable data={apiPermissionsData} fullscreen fullscreenTitle={t('权限')} />
            </Info_Section>
          )}

          {/* Tools - for MCP sources */}
          {source.config.type === 'mcp' && (
            <Info_Section
              title={t('工具')}
              description={t('此服务器提供的操作')}
              actions={
                // EditPopover for AI-assisted tool permissions editing
                <EditPopover
                  trigger={<EditButton />}
                  {...getEditConfig('source-tool-permissions', source.folderPath)}
                  secondaryAction={{
                    label: t('编辑文件'),
                    onClick: handleEditPermissions,
                  }}
                />
              }
            >
              <ToolsDataTable
                data={toolsData}
                loading={mcpToolsLoading}
                error={mcpToolsError ?? undefined}
              />
            </Info_Section>
          )}

          {/* Permissions - for MCP sources */}
          {source.config.type === 'mcp' && permissionsConfig && mcpPermissionsData.length > 0 && (
            <Info_Section
              title={t('权限')}
              description={getPermissionsDescription(source, t)}
              actions={
                // EditPopover for AI-assisted permissions.json editing
                <EditPopover
                  trigger={<EditButton />}
                  {...getEditConfig('source-permissions', source.folderPath)}
                  secondaryAction={{
                    label: t('编辑文件'),
                    onClick: handleEditPermissions,
                  }}
                />
              }
            >
              <PermissionsDataTable data={mcpPermissionsData} hideTypeColumn fullscreen fullscreenTitle={t('权限')} />
            </Info_Section>
          )}

          {/* Documentation */}
          {source.guide?.raw && (
            <Info_Section
              title={t('文档')}
              description={t('智能体的上下文和指南')}
              actions={
                // EditPopover for AI-assisted guide.md editing with "Edit File" as secondary action
                <EditPopover
                  trigger={<EditButton />}
                  {...getEditConfig('source-guide', source.folderPath)}
                  secondaryAction={{
                    label: t('编辑文件'),
                    onClick: handleEditGuide,
                  }}
                />
              }
            >
              <Info_Markdown maxHeight={540} fullscreen>
                {source.guide.raw}
              </Info_Markdown>
            </Info_Section>
          )}
        </Info_Page.Content>
      )}
    </Info_Page>
  )
}
