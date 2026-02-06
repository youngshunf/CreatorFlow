/**
 * MarketplaceInfoPage
 *
 * Detail page for displaying marketplace skill or app information.
 * Supports installation, version history, and related info.
 */

import * as React from 'react'
import { useState, useEffect, useCallback } from 'react'
import { 
  Zap, 
  Package, 
  Download, 
  RefreshCw, 
  CheckCircle2, 
  User, 
  Calendar, 
  Tag,
  ArrowLeft,
  ExternalLink,
  History,
} from 'lucide-react'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { Progress } from '@/components/ui/progress'
import { cn } from '@/lib/utils'
import { useT } from '@/context/LocaleContext'
import { navigate, routes } from '@/lib/navigate'
import type { 
  MarketplaceSkill, 
  MarketplaceApp,
  MarketplaceSkillVersion,
  MarketplaceAppVersion,
  InstallProgress,
  InstalledSkillInfo,
} from '@creator-flow/shared/marketplace'

export interface MarketplaceInfoPageProps {
  type: 'skill' | 'app'
  itemId: string
  workspaceId?: string
  className?: string
}

export function MarketplaceInfoPage({
  type,
  itemId,
  workspaceId,
  className,
}: MarketplaceInfoPageProps) {
  const t = useT()
  const [skill, setSkill] = useState<MarketplaceSkill | null>(null)
  const [app, setApp] = useState<MarketplaceApp | null>(null)
  const [skillVersions, setSkillVersions] = useState<MarketplaceSkillVersion[]>([])
  const [appVersions, setAppVersions] = useState<MarketplaceAppVersion[]>([])
  const [installedInfo, setInstalledInfo] = useState<InstalledSkillInfo | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [isInstalling, setIsInstalling] = useState(false)
  const [installProgress, setInstallProgress] = useState<InstallProgress | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [imgError, setImgError] = useState(false)

  // Load item details
  const loadItemDetails = useCallback(async () => {
    if (!itemId) return
    setIsLoading(true)
    setError(null)
    try {
      if (type === 'skill') {
        const result = await window.electronAPI.marketplaceGetSkill(itemId)
        setSkill(result.skill)
        setSkillVersions(result.versions)
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
  }, [itemId, type, workspaceId])

  useEffect(() => {
    loadItemDetails()
  }, [loadItemDetails])

  // Listen for install progress events
  useEffect(() => {
    if (!workspaceId) return
    const cleanup = window.electronAPI.onMarketplaceInstallProgress((progress: InstallProgress) => {
      setInstallProgress(progress)
      if (progress.stage === 'complete' || progress.stage === 'error') {
        setIsInstalling(false)
        if (progress.stage === 'complete') {
          loadItemDetails() // Refresh to update installed status
        }
      }
    })
    return cleanup
  }, [workspaceId, loadItemDetails])

  // Reset imgError when item changes
  useEffect(() => {
    setImgError(false)
  }, [skill?.icon_url, app?.icon_url])

  // Install handler
  const handleInstall = async () => {
    if (!workspaceId || !itemId) return
    setIsInstalling(true)
    setInstallProgress(null)
    try {
      const latestVersion = skillVersions.find(v => v.is_latest)?.version
      await window.electronAPI.marketplaceInstallSkill(workspaceId, itemId, latestVersion)
    } catch (err) {
      console.error('Install failed:', err)
      setIsInstalling(false)
      setInstallProgress({
        stage: 'error',
        percent: 0,
        message: t('安装失败'),
        error: err instanceof Error ? err.message : 'Install failed',
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

  // Loading state
  if (isLoading) {
    return (
      <div className={cn('flex-1 flex items-center justify-center', className)}>
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  // Error state
  if (error) {
    return (
      <div className={cn('flex-1 flex flex-col items-center justify-center gap-4', className)}>
        <Package className="h-12 w-12 text-muted-foreground" />
        <div className="text-center">
          <div className="font-medium">{t('加载失败')}</div>
          <div className="text-sm text-muted-foreground">{error}</div>
        </div>
        <Button variant="outline" onClick={() => navigate(routes.view.marketplace())}>
          <ArrowLeft className="h-4 w-4 mr-2" />
          {t('返回广场')}
        </Button>
      </div>
    )
  }

  // Not found
  if ((type === 'skill' && !skill) || (type === 'app' && !app)) {
    return (
      <div className={cn('flex-1 flex flex-col items-center justify-center gap-4', className)}>
        <Package className="h-12 w-12 text-muted-foreground" />
        <div className="text-center">
          <div className="font-medium">{t('未找到')}</div>
          <div className="text-sm text-muted-foreground">{t('该项目不存在或已被删除')}</div>
        </div>
        <Button variant="outline" onClick={() => navigate(routes.view.marketplace())}>
          <ArrowLeft className="h-4 w-4 mr-2" />
          {t('返回广场')}
        </Button>
      </div>
    )
  }

  const item = type === 'skill' ? skill! : app!
  const versions = type === 'skill' ? skillVersions : appVersions
  const latestVersion = versions.find(v => v.is_latest)
  const isInstalled = installedInfo !== null
  const hasUpdate = isInstalled && installedInfo.hasUpdate

  return (
    <ScrollArea className={cn('flex-1 h-full', className)}>
      <div className="max-w-3xl mx-auto p-6">
        {/* Header */}
        <div className="flex items-start gap-4 mb-6">
          {/* Icon */}
          <div className="shrink-0">
            {item.icon_url && !imgError ? (
              <img 
                src={item.icon_url} 
                alt="" 
                className="w-16 h-16 rounded-xl object-cover" 
                onError={() => setImgError(true)}
              />
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
              <h1 className="text-xl font-semibold">{item.name}</h1>
              {latestVersion && (
                <Badge variant="secondary" className="bg-foreground/5 text-muted-foreground border-0">
                  v{latestVersion.version}
                </Badge>
              )}
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
            <div className="text-sm text-muted-foreground mt-1">
              {item.description || t('暂无描述')}
            </div>
            <div className="flex items-center gap-4 mt-3 text-sm text-muted-foreground">
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

          {/* Install/Update button */}
          <div className="shrink-0">
            {isInstalling ? (
              <div className="w-32">
                <div className="text-xs text-muted-foreground mb-1.5">
                  {installProgress?.message || t('准备中...')}
                </div>
                <Progress value={installProgress?.percent || 0} className="h-1.5" />
              </div>
            ) : hasUpdate ? (
              <Button onClick={handleUpdate}>
                <RefreshCw className="h-4 w-4 mr-2" />
                {t('更新')}
              </Button>
            ) : isInstalled ? (
              <Button variant="outline" onClick={handleInstall}>
                <RefreshCw className="h-4 w-4 mr-2" />
                {t('重新安装')}
              </Button>
            ) : (
              <Button onClick={handleInstall}>
                <Download className="h-4 w-4 mr-2" />
                {item.pricing_type === 'free' ? t('免费安装') : `¥${item.price}`}
              </Button>
            )}
          </div>
        </div>

        <Separator className="my-6" />

        {/* Version info */}
        {latestVersion && (
          <div className="mb-6">
            <h2 className="text-sm font-medium mb-3 flex items-center gap-2">
              <History className="h-4 w-4" />
              {t('最新版本')}
            </h2>
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
          <div className="mb-6">
            <h2 className="text-sm font-medium mb-3">{t('版本历史')}</h2>
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
          <div className="mb-6">
            <h2 className="text-sm font-medium mb-3">{t('标签')}</h2>
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
          <div className="mb-6">
            <h2 className="text-sm font-medium mb-3">{t('依赖技能')}</h2>
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
          <div className="mb-6">
            <h2 className="text-sm font-medium mb-3">{t('安装信息')}</h2>
            <div className="p-4 rounded-lg bg-foreground/[0.02] border border-foreground/5 space-y-2">
              <div className="flex items-center justify-between text-sm">
                <span className="text-muted-foreground">{t('已安装版本')}</span>
                <span>v{installedInfo.version}</span>
              </div>
              <div className="flex items-center justify-between text-sm">
                <span className="text-muted-foreground">{t('本地路径')}</span>
                <button
                  className="text-xs text-accent hover:underline flex items-center gap-1"
                  onClick={() => window.electronAPI.openFile(installedInfo.path)}
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
  )
}
