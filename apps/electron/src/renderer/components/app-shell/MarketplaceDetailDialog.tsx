/**
 * MarketplaceDetailDialog
 *
 * Dialog component for displaying marketplace skill or app details.
 * Supports installation, version history, and related info.
 */

import * as React from 'react'
import { useState, useEffect } from 'react'
import { toast } from 'sonner'
import { 
  Zap, 
  Package, 
  Download, 
  RefreshCw, 
  CheckCircle2, 
  User, 
  Calendar, 
  Tag,
  ExternalLink,
  History,
  X,
} from 'lucide-react'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
} from '@/components/ui/dialog'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { Progress } from '@/components/ui/progress'
import { cn } from '@/lib/utils'
import { useT } from '@/context/LocaleContext'
import { useAppShellContext } from '@/context/AppShellContext'
import type {
  MarketplaceSkill, 
  MarketplaceApp,
  MarketplaceSkillVersion,
  MarketplaceAppVersion,
  InstallProgress,
  InstalledSkillInfo,
} from '@creator-flow/shared/marketplace'

export interface MarketplaceDetailDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  type: 'skill' | 'app' | null
  itemId: string | null
  workspaceId?: string
  /** Callback when user clicks "Use" for an app */
  onUseApp?: (appId: string, appName: string) => void
}

export function MarketplaceDetailDialog({
  open,
  onOpenChange,
  type,
  itemId,
  workspaceId,
  onUseApp,
}: MarketplaceDetailDialogProps) {
  const t = useT()
  const { onOpenFile } = useAppShellContext()
  const [skill, setSkill] = useState<MarketplaceSkill | null>(null)
  const [app, setApp] = useState<MarketplaceApp | null>(null)
  const [skillVersions, setSkillVersions] = useState<MarketplaceSkillVersion[]>([])
  const [appVersions, setAppVersions] = useState<MarketplaceAppVersion[]>([])
  const [installedInfo, setInstalledInfo] = useState<InstalledSkillInfo | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [isInstalling, setIsInstalling] = useState(false)
  const [installProgress, setInstallProgress] = useState<InstallProgress | null>(null)
  const [error, setError] = useState<string | null>(null)

  // Load item details when dialog opens with valid item
  useEffect(() => {
    const loadItemDetails = async () => {
      if (!open || !itemId || !type) return
      
      setIsLoading(true)
      setError(null)
      try {
        if (type === 'skill') {
          const result = await window.electronAPI.marketplaceGetSkill(itemId)
          console.log('[MarketplaceDetailDialog] API result:', result)
          setSkill(result.skill)
          setSkillVersions(result.versions || [])
          // Check if installed
          if (workspaceId) {
            const installed = await window.electronAPI.marketplaceGetInstalled(workspaceId)
            const found = installed.find(s => s.skillId === itemId)
            setInstalledInfo(found || null)
          }
        } else {
          const result = await window.electronAPI.marketplaceGetApp(itemId)
          setApp(result.app)
          setAppVersions(result.versions)
        }
      } catch (err) {
        console.error('Failed to load item details:', err)
        setError(err instanceof Error ? err.message : 'Failed to load details')
      } finally {
        setIsLoading(false)
      }
    }

    loadItemDetails()
  }, [open, itemId, type, workspaceId])

  // Reset state when dialog closes or opens
  useEffect(() => {
    if (!open) {
      // Reset all state when closing
      setSkill(null)
      setApp(null)
      setSkillVersions([])
      setAppVersions([])
      setInstalledInfo(null)
      setInstallProgress(null)
      setError(null)
      setIsLoading(true) // Reset to loading state for next open
    }
  }, [open])

  // Listen for install progress events
  useEffect(() => {
    if (!workspaceId || !open) return
    const cleanup = window.electronAPI.onMarketplaceInstallProgress((progress: InstallProgress) => {
      setInstallProgress(progress)
      if (progress.stage === 'complete' || progress.stage === 'error') {
        setIsInstalling(false)
        if (progress.stage === 'complete') {
          // Show success toast and update installed status
          toast.success(t('安装成功'), {
            description: `${skill?.name || itemId} ${t('已安装到工作区')}`,
          })
          // Update installedInfo to reflect installed state
          setInstalledInfo({
            skillId: itemId!,
            version: progress.skillId ? 'latest' : 'latest',
            path: '',
            hasUpdate: false,
            isModified: false,
          })
        } else if (progress.stage === 'error') {
          toast.error(t('安装失败'), {
            description: progress.error || t('未知错误'),
          })
        }
      }
    })
    return cleanup
  }, [workspaceId, open, itemId, skill?.name, t])

  // Install handler
  const handleInstall = async () => {
    console.log('[MarketplaceDetailDialog] handleInstall called:', { workspaceId, itemId, type })
    if (!workspaceId || !itemId) {
      console.warn('[MarketplaceDetailDialog] Missing workspaceId or itemId')
      return
    }
    setIsInstalling(true)
    setInstallProgress(null)
    try {
      const latestVersion = skillVersions.find(v => v.is_latest)?.version || skillVersions[0]?.version
      const installResult = await window.electronAPI.marketplaceInstallSkill(workspaceId, itemId, latestVersion)
      if (!installResult.success) {
        setIsInstalling(false)
        toast.error(t('安装失败'), {
          description: installResult.error || t('未知错误'),
        })
      }
    } catch (err) {
      console.error('[MarketplaceDetailDialog] Install failed:', err)
      setIsInstalling(false)
      toast.error(t('安装失败'), {
        description: err instanceof Error ? err.message : t('安装失败'),
      })
    }
  }

  // Update handler
  const handleUpdate = async () => {
    if (!workspaceId || !itemId) return
    setIsInstalling(true)
    setInstallProgress(null)
    try {
      const latestVersion = skillVersions.find(v => v.is_latest)?.version
      await window.electronAPI.marketplaceUpdateSkill(workspaceId, itemId, latestVersion)
    } catch (err) {
      console.error('Update failed:', err)
      setIsInstalling(false)
    }
  }

  const item = type === 'skill' ? skill : app
  const versions = type === 'skill' ? skillVersions : appVersions
  const latestVersion = versions?.find(v => v.is_latest)
  const isInstalled = installedInfo !== null
  const hasUpdate = isInstalled && installedInfo?.hasUpdate

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-w-3xl max-h-[85vh] p-0 gap-0 overflow-hidden" showCloseButton={false} aria-describedby={undefined}>
        {/* Loading state */}
        {isLoading && (
          <>
            <DialogTitle className="sr-only">{t('加载中')}</DialogTitle>
            <div className="flex items-center justify-center h-80">
              <RefreshCw className="h-6 w-6 animate-spin text-muted-foreground" />
            </div>
          </>
        )}

        {/* Error state */}
        {!isLoading && error && (
          <>
            <DialogTitle className="sr-only">{t('加载失败')}</DialogTitle>
            <div className="flex flex-col items-center justify-center h-80 gap-4">
              <Package className="h-16 w-16 text-muted-foreground" />
              <div className="text-center">
                <div className="text-lg font-medium">{t('加载失败')}</div>
                <div className="text-sm text-muted-foreground mt-1">{error}</div>
              </div>
            </div>
          </>
        )}

        {/* Not found */}
        {!isLoading && !error && !item && (
          <>
            <DialogTitle className="sr-only">{t('未找到')}</DialogTitle>
            <div className="flex flex-col items-center justify-center h-80 gap-4">
              <Package className="h-16 w-16 text-muted-foreground" />
              <div className="text-center">
                <div className="text-lg font-medium">{t('未找到')}</div>
                <div className="text-sm text-muted-foreground mt-1">{t('该项目不存在或已被删除')}</div>
              </div>
            </div>
          </>
        )}

        {/* Content */}
        {!isLoading && !error && item && (
          <>
            {/* Header - Icon, Title, Badges */}
            <DialogHeader className="p-6 pb-4">
              <div className="flex items-start gap-4">
                {/* Icon */}
                <div className="shrink-0">
                  {item.icon_url ? (
                    <img src={item.icon_url} alt="" className="w-16 h-16 rounded-xl" />
                  ) : (
                    <div className={cn(
                      "w-16 h-16 rounded-xl flex items-center justify-center",
                      type === 'skill' ? "bg-accent/10" : "bg-purple-500/10"
                    )}>
                      {type === 'skill' ? (
                        <Zap className="h-8 w-8 text-accent" />
                      ) : (
                        <Package className="h-8 w-8 text-purple-500" />
                      )}
                    </div>
                  )}
                </div>

                {/* Title & meta */}
                <div className="flex-1 min-w-0">
                  <div className="flex items-center gap-2 flex-wrap">
                    <DialogTitle className="text-xl font-semibold">{item.name}</DialogTitle>
                    {item.is_official && (
                      <Badge variant="secondary" className="bg-blue-500/10 text-blue-500 border-0">
                        {t('官方')}
                      </Badge>
                    )}
                    {isInstalled && (
                      <Badge variant="secondary" className="bg-green-500/10 text-green-500 border-0">
                        <CheckCircle2 className="h-3 w-3 mr-1" />
                        {t('已安装')}
                      </Badge>
                    )}
                  </div>
                  <div className="flex items-center gap-4 mt-2 text-sm text-muted-foreground">
                    {item.author_name && (
                      <span className="flex items-center gap-1">
                        <User className="h-3.5 w-3.5" />
                        {item.author_name}
                      </span>
                    )}
                    <span className="flex items-center gap-1">
                      <Download className="h-3.5 w-3.5" />
                      {item.download_count} {t('次下载')}
                    </span>
                    {item.category && (
                      <span className="flex items-center gap-1">
                        <Tag className="h-3.5 w-3.5" />
                        {item.category}
                      </span>
                    )}
                  </div>
                </div>
              </div>
            </DialogHeader>

            <Separator />

            {/* Scrollable content */}
            <ScrollArea className="flex-1 max-h-[55vh]">
              <div className="p-6 space-y-6">
                {/* Description */}
                {item.description && (
                  <div>
                    <h3 className="text-sm font-medium mb-3">{t('介绍')}</h3>
                    <div className="text-sm text-muted-foreground whitespace-pre-wrap leading-relaxed">
                      {item.description}
                    </div>
                  </div>
                )}

                {/* Version info */}
                {latestVersion && (
                  <div>
                    <h3 className="text-sm font-medium mb-3 flex items-center gap-2">
                      <History className="h-4 w-4" />
                      {t('最新版本')}
                    </h3>
                    <div className="p-4 rounded-lg bg-foreground/[0.02] border border-foreground/5">
                      <div className="flex items-center justify-between mb-2">
                        <span className="font-medium">v{latestVersion.version}</span>
                        <span className="text-xs text-muted-foreground">
                          <Calendar className="h-3 w-3 inline mr-1" />
                          {new Date(latestVersion.published_at).toLocaleDateString()}
                        </span>
                      </div>
                      {latestVersion.changelog && (
                        <div className="text-sm text-muted-foreground whitespace-pre-wrap">
                          {latestVersion.changelog}
                        </div>
                      )}
                    </div>
                  </div>
                )}

                {/* Version history */}
                {versions.length > 1 && (
                  <div>
                    <h3 className="text-sm font-medium mb-3">{t('版本历史')}</h3>
                    <div className="space-y-2">
                      {versions.slice(1, 5).map((version) => (
                        <div
                          key={'version' in version ? version.version : version.version}
                          className="p-3 rounded-lg bg-foreground/[0.02] border border-foreground/5"
                        >
                          <div className="flex items-center justify-between">
                            <span className="text-sm">v{version.version}</span>
                            <span className="text-xs text-muted-foreground">
                              {new Date(version.published_at).toLocaleDateString()}
                            </span>
                          </div>
                          {version.changelog && (
                            <div className="text-xs text-muted-foreground mt-1 line-clamp-2">
                              {version.changelog}
                            </div>
                          )}
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {/* Tags */}
                {type === 'skill' && skill?.tags && (
                  <div>
                    <h3 className="text-sm font-medium mb-3">{t('标签')}</h3>
                    <div className="flex flex-wrap gap-2">
                      {skill.tags.split(',').map((tag) => (
                        <Badge key={tag} variant="secondary" className="text-xs">
                          {tag.trim()}
                        </Badge>
                      ))}
                    </div>
                  </div>
                )}

                {/* App dependencies */}
                {type === 'app' && app?.skill_dependencies && (
                  <div>
                    <h3 className="text-sm font-medium mb-3">{t('依赖技能')}</h3>
                    <div className="p-4 rounded-lg bg-foreground/[0.02] border border-foreground/5">
                      <div className="text-sm text-muted-foreground">
                        {app.skill_dependencies.split(',').map((dep) => (
                          <div key={dep} className="flex items-center gap-2 py-1">
                            <Zap className="h-3.5 w-3.5 text-accent" />
                            <span>{dep.trim()}</span>
                          </div>
                        ))}
                      </div>
                    </div>
                  </div>
                )}

                {/* Installed info */}
                {isInstalled && installedInfo && (
                  <div>
                    <h3 className="text-sm font-medium mb-3">{t('安装信息')}</h3>
                    <div className="p-4 rounded-lg bg-foreground/[0.02] border border-foreground/5 space-y-2">
                      <div className="flex items-center justify-between text-sm">
                        <span className="text-muted-foreground">{t('已安装版本')}</span>
                        <span>v{installedInfo.version}</span>
                      </div>
                      <div className="flex items-center justify-between text-sm">
                        <span className="text-muted-foreground">{t('本地路径')}</span>
                        <button
                          className="text-xs text-accent hover:underline flex items-center gap-1"
                          onClick={() => onOpenFile?.(installedInfo.path)}
                        >
                          {t('打开文件夹')}
                          <ExternalLink className="h-3 w-3" />
                        </button>
                      </div>
                      {installedInfo.isModified && (
                        <div className="flex items-center gap-2 text-sm text-amber-500">
                          <span>{t('本地已修改')}</span>
                        </div>
                      )}
                    </div>
                  </div>
                )}
              </div>
            </ScrollArea>

            {/* Footer - Install button */}
            <div className="border-t border-border/40 p-4 flex items-center justify-between bg-background">
              <div className="text-sm text-muted-foreground">
                {item.pricing_type === 'free' ? (
                  <span className="text-green-500 font-medium">{t('免费')}</span>
                ) : (
                  <span className="text-amber-500 font-medium">¥{item.price}</span>
                )}
              </div>
              <div className="flex items-center gap-3">
                {isInstalling ? (
                  <div className="flex items-center gap-3">
                    <div className="text-sm text-muted-foreground">
                      {installProgress?.message || t('准备中...')}
                    </div>
                    <div className="w-32">
                      <Progress value={installProgress?.percent || 0} className="h-2" />
                    </div>
                  </div>
                ) : (
                  <>
                    <Button variant="outline" onClick={() => onOpenChange(false)}>
                      {t('关闭')}
                    </Button>
                    {type === 'app' ? (
                      // Apps: "Use" button to create workspace
                      <Button onClick={() => {
                        if (app && onUseApp) {
                          onOpenChange(false)
                          onUseApp(app.app_id, app.name)
                        }
                      }}>
                        <Package className="h-4 w-4 mr-2" />
                        {t('使用')}
                      </Button>
                    ) : hasUpdate ? (
                      <Button onClick={handleUpdate}>
                        <RefreshCw className="h-4 w-4 mr-2" />
                        {t('更新')}
                      </Button>
                    ) : !isInstalled && (
                      <Button onClick={handleInstall}>
                        <Download className="h-4 w-4 mr-2" />
                        {t('安装')}
                      </Button>
                    )}
                  </>
                )}
              </div>
            </div>
          </>
        )}
      </DialogContent>
    </Dialog>
  )
}
