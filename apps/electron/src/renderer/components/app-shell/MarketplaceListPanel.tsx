/**
 * MarketplaceListPanel
 *
 * Panel component for displaying cloud marketplace skills/apps in the sidebar.
 * Supports filtering by skills, apps, or category.
 */

import * as React from 'react'
import { useState, useEffect, useCallback } from 'react'
import { Search, Zap, Package, Download, RefreshCw, Tag, MoreHorizontal, CheckCircle2 } from 'lucide-react'
import { SkillAvatar } from '@/components/ui/skill-avatar'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Empty, EmptyHeader, EmptyMedia, EmptyTitle, EmptyDescription } from '@/components/ui/empty'
import { Separator } from '@/components/ui/separator'
import { Input } from '@/components/ui/input'
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
  DropdownMenu,
  DropdownMenuTrigger,
  StyledDropdownMenuContent,
  StyledDropdownMenuItem,
} from '@/components/ui/styled-dropdown'
import { DropdownMenuProvider } from '@/components/ui/menu-context'
import { cn } from '@/lib/utils'
import { useT } from '@/context/LocaleContext'
import type { MarketplaceFilter } from '../../../shared/types'
import type { 
  MarketplaceSkill, 
  MarketplaceApp, 
  MarketplaceCategory,
  InstalledSkillInfo,
} from '@creator-flow/shared/marketplace'

export interface MarketplaceListPanelProps {
  filter: MarketplaceFilter
  onFilterChange: (filter: MarketplaceFilter) => void
  onSkillClick: (skillId: string) => void
  onAppClick: (appId: string) => void
  selectedSkillId?: string | null
  selectedAppId?: string | null
  workspaceId?: string
  className?: string
}

export function MarketplaceListPanel({
  filter,
  onFilterChange,
  onSkillClick,
  onAppClick,
  selectedSkillId,
  selectedAppId,
  workspaceId,
  className,
}: MarketplaceListPanelProps) {
  const t = useT()
  const [skills, setSkills] = useState<MarketplaceSkill[]>([])
  const [apps, setApps] = useState<MarketplaceApp[]>([])
  const [categories, setCategories] = useState<MarketplaceCategory[]>([])
  const [installedSkills, setInstalledSkills] = useState<InstalledSkillInfo[]>([])
  const [searchQuery, setSearchQuery] = useState('')
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  // Load marketplace data
  const loadMarketplaceData = useCallback(async () => {
    if (!workspaceId) return
    setIsLoading(true)
    setError(null)
    try {
      const [skillsResult, appsResult, categoriesResult, installedResult] = await Promise.all([
        window.electronAPI.marketplaceListSkills(),
        window.electronAPI.marketplaceListApps(),
        window.electronAPI.marketplaceListCategories(),
        window.electronAPI.marketplaceGetInstalled(workspaceId),
      ])
      setSkills(skillsResult.items)
      setApps(appsResult.items)
      setCategories(categoriesResult)
      setInstalledSkills(installedResult)
    } catch (err) {
      console.error('Failed to load marketplace data:', err)
      setError(err instanceof Error ? err.message : 'Failed to load marketplace data')
    } finally {
      setIsLoading(false)
    }
  }, [workspaceId])

  useEffect(() => {
    loadMarketplaceData()
  }, [loadMarketplaceData])

  // Search handler
  const handleSearch = useCallback(async (query: string) => {
    if (!query.trim()) {
      loadMarketplaceData()
      return
    }
    setIsLoading(true)
    try {
      const result = await window.electronAPI.marketplaceSearch(query)
      setSkills(result.skills)
      setApps(result.apps)
    } catch (err) {
      console.error('Search failed:', err)
    } finally {
      setIsLoading(false)
    }
  }, [loadMarketplaceData])

  // Filter items based on current filter
  const filteredSkills = React.useMemo(() => {
    if (filter.kind === 'apps') return []
    if (filter.kind === 'category' && filter.categoryId) {
      return skills.filter(s => s.category === filter.categoryId)
    }
    return skills
  }, [skills, filter])

  const filteredApps = React.useMemo(() => {
    if (filter.kind === 'skills') return []
    return apps
  }, [apps, filter])

  // Check if a skill is installed
  const isSkillInstalled = useCallback((skillId: string) => {
    return installedSkills.some(s => s.skillId === skillId)
  }, [installedSkills])

  // Tab change handler
  const handleTabChange = (value: string) => {
    if (value === 'all') {
      onFilterChange({ kind: 'all' })
    } else if (value === 'skills') {
      onFilterChange({ kind: 'skills' })
    } else if (value === 'apps') {
      onFilterChange({ kind: 'apps' })
    }
  }

  // Loading state
  if (isLoading && skills.length === 0 && apps.length === 0) {
    return (
      <div className={cn('flex-1 flex items-center justify-center', className)}>
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  // Error state
  if (error) {
    return (
      <Empty className={cn('flex-1', className)}>
        <EmptyHeader>
          <EmptyMedia variant="icon">
            <Package />
          </EmptyMedia>
          <EmptyTitle>{t('加载失败')}</EmptyTitle>
          <EmptyDescription>{error}</EmptyDescription>
        </EmptyHeader>
      </Empty>
    )
  }

  // Empty state
  if (skills.length === 0 && apps.length === 0) {
    return (
      <Empty className={cn('flex-1', className)}>
        <EmptyHeader>
          <EmptyMedia variant="icon">
            <Package />
          </EmptyMedia>
          <EmptyTitle>{t('广场暂无内容')}</EmptyTitle>
          <EmptyDescription>
            {t('探索并安装技能和应用以增强您的智能体能力。')}
          </EmptyDescription>
        </EmptyHeader>
      </Empty>
    )
  }

  return (
    <div className={cn('flex flex-col flex-1 h-full', className)}>
      {/* Search bar */}
      <div className="px-3 py-2">
        <div className="relative">
          <Search className="absolute left-2.5 top-1/2 -translate-y-1/2 h-3.5 w-3.5 text-muted-foreground" />
          <Input
            placeholder={t('搜索技能和应用...')}
            value={searchQuery}
            onChange={(e) => {
              setSearchQuery(e.target.value)
              handleSearch(e.target.value)
            }}
            className="pl-8 h-8 text-sm"
          />
        </div>
      </div>

      {/* Filter tabs */}
      <div className="px-3 pb-2">
        <Tabs value={filter.kind} onValueChange={handleTabChange}>
          <TabsList className="h-7 p-0.5 w-full">
            <TabsTrigger value="all" className="flex-1 h-6 text-xs">
              {t('全部')}
            </TabsTrigger>
            <TabsTrigger value="skills" className="flex-1 h-6 text-xs">
              <Zap className="h-3 w-3 mr-1" />
              {t('技能')}
            </TabsTrigger>
            <TabsTrigger value="apps" className="flex-1 h-6 text-xs">
              <Package className="h-3 w-3 mr-1" />
              {t('应用')}
            </TabsTrigger>
          </TabsList>
        </Tabs>
      </div>

      {/* Category chips */}
      {categories.length > 0 && filter.kind !== 'apps' && (
        <div className="px-3 pb-2 flex flex-wrap gap-1.5">
          <button
            className={cn(
              "inline-flex items-center h-6 px-2 text-xs rounded-full transition-colors",
              filter.kind === 'all' || (filter.kind !== 'category')
                ? "bg-accent text-accent-foreground"
                : "bg-foreground/5 hover:bg-foreground/10"
            )}
            onClick={() => onFilterChange({ kind: 'all' })}
          >
            {t('全部')}
          </button>
          {categories.slice(0, 5).map((cat) => (
            <button
              key={cat.id}
              className={cn(
                "inline-flex items-center h-6 px-2 text-xs rounded-full transition-colors",
                filter.kind === 'category' && filter.categoryId === cat.slug
                  ? "bg-accent text-accent-foreground"
                  : "bg-foreground/5 hover:bg-foreground/10"
              )}
              onClick={() => onFilterChange({ kind: 'category', categoryId: cat.slug })}
            >
              <Tag className="h-3 w-3 mr-1" />
              {cat.name}
            </button>
          ))}
        </div>
      )}

      <Separator className="mb-2" />

      {/* Items list */}
      <ScrollArea className="flex-1">
        <div className="pb-2">
          {/* Skills section */}
          {filteredSkills.length > 0 && (
            <div className="pt-1">
              {filter.kind === 'all' && (
                <div className="px-4 py-1.5 text-xs font-medium text-muted-foreground flex items-center gap-1.5">
                  <Zap className="h-3 w-3" />
                  {t('技能')}
                  <span className="text-foreground/50">({filteredSkills.length})</span>
                </div>
              )}
              {filteredSkills.map((skill, index) => (
                <MarketplaceSkillItem
                  key={skill.skill_id}
                  skill={skill}
                  isSelected={selectedSkillId === skill.skill_id}
                  isInstalled={isSkillInstalled(skill.skill_id)}
                  isFirst={index === 0 || filter.kind === 'all'}
                  onClick={() => onSkillClick(skill.skill_id)}
                  workspaceId={workspaceId}
                />
              ))}
            </div>
          )}

          {/* Apps section */}
          {filteredApps.length > 0 && (
            <div className={cn("pt-1", filteredSkills.length > 0 && "mt-2")}>
              {filter.kind === 'all' && (
                <div className="px-4 py-1.5 text-xs font-medium text-muted-foreground flex items-center gap-1.5">
                  <Package className="h-3 w-3" />
                  {t('应用')}
                  <span className="text-foreground/50">({filteredApps.length})</span>
                </div>
              )}
              {filteredApps.map((app, index) => (
                <MarketplaceAppItem
                  key={app.app_id}
                  app={app}
                  isSelected={selectedAppId === app.app_id}
                  isFirst={index === 0 || filter.kind === 'all'}
                  onClick={() => onAppClick(app.app_id)}
                />
              ))}
            </div>
          )}
        </div>
      </ScrollArea>
    </div>
  )
}

interface MarketplaceSkillItemProps {
  skill: MarketplaceSkill
  isSelected: boolean
  isInstalled: boolean
  isFirst: boolean
  workspaceId?: string
  onClick: () => void
}

function MarketplaceSkillItem({ 
  skill, 
  isSelected, 
  isInstalled,
  isFirst, 
  workspaceId,
  onClick, 
}: MarketplaceSkillItemProps) {
  const t = useT()
  const [menuOpen, setMenuOpen] = useState(false)
  const [imgError, setImgError] = useState(false)

  return (
    <div className="marketplace-item" data-selected={isSelected || undefined}>
      {/* Separator - only show if not first */}
      {!isFirst && (
        <div className="marketplace-separator pl-12 pr-4">
          <Separator />
        </div>
      )}
      <div className="marketplace-content relative group select-none pl-2 mr-2">
        {/* Avatar - positioned absolutely */}
        <div className="absolute left-[18px] top-3.5 z-10 flex items-center justify-center">
          {skill.icon_url && !imgError ? (
            <img 
              src={skill.icon_url} 
              alt="" 
              className="w-5 h-5 rounded object-cover" 
              onError={() => setImgError(true)}
            />
          ) : (
            <div className="w-5 h-5 rounded bg-accent/10 flex items-center justify-center">
              <Zap className="h-3 w-3 text-accent" />
            </div>
          )}
        </div>
        {/* Main content button */}
        <button
          className={cn(
            "flex w-full items-start gap-2 pl-2 pr-4 py-3 text-left text-sm transition-all outline-none rounded-[8px]",
            isSelected
              ? "bg-foreground/5 hover:bg-foreground/7"
              : "hover:bg-foreground/2"
          )}
          onClick={onClick}
        >
          {/* Spacer for avatar */}
          <div className="w-5 h-5 shrink-0" />
          {/* Content column */}
          <div className="flex flex-col gap-1 min-w-0 flex-1">
            {/* Title - skill name */}
            <div className="flex items-center gap-2 w-full pr-6 min-w-0">
              <span className="font-medium font-sans line-clamp-1 min-w-0">
                {skill.name}
              </span>
              {skill.latest_version && (
                <span className="text-[10px] text-muted-foreground shrink-0">v{skill.latest_version}</span>
              )}
              {isInstalled && (
                <CheckCircle2 className="h-3.5 w-3.5 text-green-500 shrink-0" />
              )}
              {skill.is_official && (
                <span className="px-1.5 py-0.5 text-[10px] bg-blue-500/10 text-blue-500 rounded shrink-0">
                  {t('官方')}
                </span>
              )}
            </div>
            {/* Subtitle - description */}
            <div className="flex items-center gap-1.5 text-xs text-foreground/70 w-full -mb-[2px] pr-6 min-w-0">
              <span className="truncate">
                {skill.description || t('暂无描述')}
              </span>
            </div>
            {/* Meta info */}
            <div className="flex items-center gap-3 text-xs text-foreground/50 mt-0.5">
              {skill.author_name && (
                <span>{skill.author_name}</span>
              )}
              <span className="flex items-center gap-1">
                <Download className="h-3 w-3" />
                {skill.download_count}
              </span>
              {skill.pricing_type === 'free' ? (
                <span className="text-green-500">{t('免费')}</span>
              ) : (
                <span className="text-amber-500">¥{skill.price}</span>
              )}
            </div>
          </div>
        </button>
        {/* Action buttons */}
        <div
          className={cn(
            "absolute right-2 top-2 transition-opacity z-10",
            menuOpen ? "opacity-100" : "opacity-0 group-hover:opacity-100"
          )}
        >
          <div className="flex items-center rounded-[8px] overflow-hidden border border-transparent hover:border-border/50">
            <DropdownMenu modal={true} onOpenChange={setMenuOpen}>
              <DropdownMenuTrigger asChild>
                <div className="p-1.5 hover:bg-foreground/10 data-[state=open]:bg-foreground/10 cursor-pointer">
                  <MoreHorizontal className="h-4 w-4 text-muted-foreground" />
                </div>
              </DropdownMenuTrigger>
              <StyledDropdownMenuContent align="end">
                <DropdownMenuProvider>
                  <StyledDropdownMenuItem
                    onClick={async () => {
                      if (workspaceId) {
                        await window.electronAPI.marketplaceInstallSkill(workspaceId, skill.skill_id)
                      }
                    }}
                  >
                    <Download className="h-4 w-4" />
                    {isInstalled ? t('重新安装') : t('安装')}
                  </StyledDropdownMenuItem>
                </DropdownMenuProvider>
              </StyledDropdownMenuContent>
            </DropdownMenu>
          </div>
        </div>
      </div>
    </div>
  )
}

interface MarketplaceAppItemProps {
  app: MarketplaceApp
  isSelected: boolean
  isFirst: boolean
  onClick: () => void
}

function MarketplaceAppItem({ app, isSelected, isFirst, onClick }: MarketplaceAppItemProps) {
  const t = useT()
  const [imgError, setImgError] = useState(false)

  return (
    <div className="marketplace-item" data-selected={isSelected || undefined}>
      {/* Separator - only show if not first */}
      {!isFirst && (
        <div className="marketplace-separator pl-12 pr-4">
          <Separator />
        </div>
      )}
      <div className="marketplace-content relative group select-none pl-2 mr-2">
        {/* Avatar - positioned absolutely */}
        <div className="absolute left-[18px] top-3.5 z-10 flex items-center justify-center">
          {app.icon_url && !imgError ? (
            <img 
              src={app.icon_url} 
              alt="" 
              className="w-5 h-5 rounded object-cover" 
              onError={() => setImgError(true)}
            />
          ) : (
            <div className="w-5 h-5 rounded bg-purple-500/10 flex items-center justify-center">
              <Package className="h-3 w-3 text-purple-500" />
            </div>
          )}
        </div>
        {/* Main content button */}
        <button
          className={cn(
            "flex w-full items-start gap-2 pl-2 pr-4 py-3 text-left text-sm transition-all outline-none rounded-[8px]",
            isSelected
              ? "bg-foreground/5 hover:bg-foreground/7"
              : "hover:bg-foreground/2"
          )}
          onClick={onClick}
        >
          {/* Spacer for avatar */}
          <div className="w-5 h-5 shrink-0" />
          {/* Content column */}
          <div className="flex flex-col gap-1 min-w-0 flex-1">
            {/* Title - app name */}
            <div className="flex items-center gap-2 w-full pr-6 min-w-0">
              <span className="font-medium font-sans line-clamp-1 min-w-0">
                {app.name}
              </span>
              {app.latest_version && (
                <span className="text-[10px] text-muted-foreground shrink-0">v{app.latest_version}</span>
              )}
              {app.is_official && (
                <span className="px-1.5 py-0.5 text-[10px] bg-blue-500/10 text-blue-500 rounded shrink-0">
                  {t('官方')}
                </span>
              )}
            </div>
            {/* Subtitle - description */}
            <div className="flex items-center gap-1.5 text-xs text-foreground/70 w-full -mb-[2px] pr-6 min-w-0">
              <span className="truncate">
                {app.description || t('暂无描述')}
              </span>
            </div>
            {/* Meta info */}
            <div className="flex items-center gap-3 text-xs text-foreground/50 mt-0.5">
              {app.author_name && (
                <span>{app.author_name}</span>
              )}
              <span className="flex items-center gap-1">
                <Download className="h-3 w-3" />
                {app.download_count}
              </span>
              {app.pricing_type === 'free' ? (
                <span className="text-green-500">{t('免费')}</span>
              ) : app.pricing_type === 'subscription' ? (
                <span className="text-amber-500">{t('订阅')}</span>
              ) : (
                <span className="text-amber-500">¥{app.price}</span>
              )}
            </div>
          </div>
        </button>
      </div>
    </div>
  )
}
