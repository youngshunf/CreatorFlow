/**
 * MarketplacePage
 *
 * Full-page marketplace component with card grid layout.
 * Skills and apps are displayed as cards, clicking opens detail dialog.
 */

import * as React from 'react'
import { useState, useEffect, useCallback, useMemo, useRef } from 'react'
import { 
  Search, 
  Zap, 
  Package, 
  Download, 
  RefreshCw, 
  Tag, 
  CheckCircle2,
  X,
} from 'lucide-react'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Empty, EmptyHeader, EmptyMedia, EmptyTitle, EmptyDescription } from '@/components/ui/empty'
import { Input } from '@/components/ui/input'
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { cn } from '@/lib/utils'
import { useT } from '@/context/LocaleContext'
import { MarketplaceDetailDialog } from '@/components/app-shell/MarketplaceDetailDialog'
import type { MarketplaceFilter } from '../../shared/types'
import type { 
  MarketplaceSkill, 
  MarketplaceApp, 
  MarketplaceCategory,
  InstalledSkillInfo,
} from '@creator-flow/shared/marketplace'

export interface MarketplacePageProps {
  filter: MarketplaceFilter
  onFilterChange: (filter: MarketplaceFilter) => void
  workspaceId?: string
  className?: string
}

export function MarketplacePage({
  filter,
  onFilterChange,
  workspaceId,
  className,
}: MarketplacePageProps) {
  const t = useT()
  const [skills, setSkills] = useState<MarketplaceSkill[]>([])
  const [apps, setApps] = useState<MarketplaceApp[]>([])
  const [categories, setCategories] = useState<MarketplaceCategory[]>([])
  const [installedSkills, setInstalledSkills] = useState<InstalledSkillInfo[]>([])
  const [searchQuery, setSearchQuery] = useState('')
  const [isSearching, setIsSearching] = useState(false)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const searchTimerRef = useRef<NodeJS.Timeout | null>(null)
  
  // Keep original data for when search is cleared
  const [originalSkills, setOriginalSkills] = useState<MarketplaceSkill[]>([])
  const [originalApps, setOriginalApps] = useState<MarketplaceApp[]>([])
  
  // Dialog state
  const [dialogOpen, setDialogOpen] = useState(false)
  const [selectedItem, setSelectedItem] = useState<{ type: 'skill' | 'app'; id: string } | null>(null)

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
      setOriginalSkills(skillsResult.items)
      setOriginalApps(appsResult.items)
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

  // Search handler with debounce
  const handleSearchInput = useCallback((query: string) => {
    setSearchQuery(query)
    
    // Clear previous timer
    if (searchTimerRef.current) {
      clearTimeout(searchTimerRef.current)
    }
    
    // If query is empty, restore original data
    if (!query.trim()) {
      setSkills(originalSkills)
      setApps(originalApps)
      setIsSearching(false)
      return
    }
    
    // Only search if query has at least 2 characters (1 Chinese character = 1-3 bytes but 1 char length)
    if (query.trim().length < 2) {
      return
    }
    
    // Debounce search
    searchTimerRef.current = setTimeout(async () => {
      setIsSearching(true)
      try {
        const result = await window.electronAPI.marketplaceSearch(query)
        setSkills(result.skills)
        setApps(result.apps)
      } catch (err) {
        console.error('Search failed:', err)
      } finally {
        setIsSearching(false)
      }
    }, 300)
  }, [originalSkills, originalApps])
  
  // Cleanup timer on unmount
  useEffect(() => {
    return () => {
      if (searchTimerRef.current) {
        clearTimeout(searchTimerRef.current)
      }
    }
  }, [])

  // Filter items based on current filter
  const filteredSkills = useMemo(() => {
    if (filter.kind === 'apps') return []
    if (filter.kind === 'category' && filter.categoryId) {
      return skills.filter(s => s.category === filter.categoryId)
    }
    return skills
  }, [skills, filter])

  const filteredApps = useMemo(() => {
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

  // Handle item click - open dialog
  const handleSkillClick = (skillId: string) => {
    setSelectedItem({ type: 'skill', id: skillId })
    setDialogOpen(true)
  }

  const handleAppClick = (appId: string) => {
    setSelectedItem({ type: 'app', id: appId })
    setDialogOpen(true)
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

  // Empty state - only show if no original data (not during search)
  if (originalSkills.length === 0 && originalApps.length === 0 && !searchQuery) {
    return (
      <Empty className={cn('flex-1', className)}>
        <EmptyHeader>
          <EmptyMedia variant="icon">
            <Package />
          </EmptyMedia>
          <EmptyTitle>{t('市场暂无内容')}</EmptyTitle>
          <EmptyDescription>
            {t('探索并安装技能和应用以增强您的智能体能力。')}
          </EmptyDescription>
        </EmptyHeader>
      </Empty>
    )
  }

  return (
    <div className={cn('flex flex-col flex-1 h-full', className)}>
      {/* Header with search and filters */}
      <div className="px-6 pt-6 pb-4 border-b border-border/40">
        {/* Title - centered */}
        <h1 className="text-base font-semibold text-foreground mb-4 text-center">{t('市场')}</h1>
        
        {/* Search bar - centered */}
        <div className="flex justify-center mb-4">
          <div className="relative w-full max-w-md">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
            <Input
              placeholder={t('搜索技能和应用...')}
              value={searchQuery}
              onChange={(e) => handleSearchInput(e.target.value)}
              className="pl-9 pr-9 h-9"
            />
            {searchQuery && (
              <button
                type="button"
                onClick={() => handleSearchInput('')}
                className="absolute right-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground hover:text-foreground transition-colors cursor-pointer"
              >
                <X className="h-4 w-4" />
              </button>
            )}
          </div>
        </div>

        {/* Filter tabs - centered */}
        <div className="flex justify-center mb-3">
          <Tabs value={filter.kind} onValueChange={handleTabChange}>
            <TabsList className="h-9 p-1">
              <TabsTrigger value="all" className="h-7 px-6 text-sm">
                {t('全部')}
              </TabsTrigger>
              <TabsTrigger value="skills" className="h-7 px-6 text-sm">
                <Zap className="h-3.5 w-3.5 mr-1.5" />
                {t('技能')}
              </TabsTrigger>
              <TabsTrigger value="apps" className="h-7 px-6 text-sm">
                <Package className="h-3.5 w-3.5 mr-1.5" />
                {t('应用')}
              </TabsTrigger>
            </TabsList>
          </Tabs>
        </div>

        {/* Category chips - left aligned */}
        {categories.length > 0 && filter.kind !== 'apps' && (
          <div className="flex flex-wrap gap-1.5">
            <button
              className={cn(
                "inline-flex items-center h-7 px-2.5 text-xs rounded-full transition-colors cursor-pointer",
                filter.kind === 'all' || (filter.kind !== 'category')
                  ? "bg-accent text-accent-foreground"
                  : "bg-foreground/5 hover:bg-foreground/10"
              )}
              onClick={() => onFilterChange({ kind: 'all' })}
            >
              {t('全部')}
            </button>
            {categories.map((cat) => (
              <button
                key={cat.id}
                className={cn(
                  "inline-flex items-center h-7 px-2.5 text-xs rounded-full transition-colors cursor-pointer",
                  filter.kind === 'category' && filter.categoryId === cat.slug
                    ? "bg-accent text-accent-foreground"
                    : "bg-foreground/5 hover:bg-foreground/10"
                )}
                onClick={() => onFilterChange({ kind: 'category', categoryId: cat.slug })}
              >
                {cat.icon && <span className="mr-1">{cat.icon}</span>}
                {cat.name}
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Card grid */}
      <ScrollArea className="flex-1">
        <div className="p-6">
          {/* Searching indicator */}
          {isSearching && (
            <div className="flex items-center justify-center py-12">
              <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
            </div>
          )}

          {/* No search results */}
          {!isSearching && searchQuery.trim().length >= 2 && filteredSkills.length === 0 && filteredApps.length === 0 && (
            <div className="flex flex-col items-center justify-center py-12">
              <Search className="h-12 w-12 text-muted-foreground/50 mb-4" />
              <p className="text-sm text-muted-foreground">{t('无搜索结果')}</p>
              <p className="text-xs text-muted-foreground/70 mt-1">{t('请尝试其他关键词')}</p>
            </div>
          )}

          {/* Skills section */}
          {!isSearching && filteredSkills.length > 0 && (
            <div className="mb-8">
              {filter.kind === 'all' && (
                <div className="flex items-center gap-2 mb-4">
                  <Zap className="h-4 w-4 text-accent" />
                  <h2 className="text-sm font-medium">{t('技能')}</h2>
                  <span className="text-xs text-muted-foreground">({filteredSkills.length})</span>
                </div>
              )}
              <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4">
                {filteredSkills.map((skill) => (
                  <SkillCard
                    key={skill.skill_id}
                    skill={skill}
                    isInstalled={isSkillInstalled(skill.skill_id)}
                    onClick={() => handleSkillClick(skill.skill_id)}
                  />
                ))}
              </div>
            </div>
          )}

          {/* Apps section */}
          {!isSearching && filteredApps.length > 0 && (
            <div>
              {filter.kind === 'all' && (
                <div className="flex items-center gap-2 mb-4">
                  <Package className="h-4 w-4 text-purple-500" />
                  <h2 className="text-sm font-medium">{t('应用')}</h2>
                  <span className="text-xs text-muted-foreground">({filteredApps.length})</span>
                </div>
              )}
              <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4">
                {filteredApps.map((app) => (
                  <AppCard
                    key={app.app_id}
                    app={app}
                    onClick={() => handleAppClick(app.app_id)}
                  />
                ))}
              </div>
            </div>
          )}
        </div>
      </ScrollArea>

      {/* Detail dialog */}
      <MarketplaceDetailDialog
        open={dialogOpen}
        onOpenChange={setDialogOpen}
        type={selectedItem?.type || null}
        itemId={selectedItem?.id || null}
        workspaceId={workspaceId}
      />
    </div>
  )
}

interface SkillCardProps {
  skill: MarketplaceSkill
  isInstalled: boolean
  onClick: () => void
}

function SkillCard({ skill, isInstalled, onClick }: SkillCardProps) {
  const t = useT()

  return (
    <button
      type="button"
      onClick={onClick}
      className="group flex flex-col items-start rounded-xl border border-border/60 bg-background/40 p-4 text-left shadow-[0_0_0_1px_rgba(15,23,42,0.02)] hover:border-border hover:bg-foreground/[0.02] transition-colors cursor-pointer"
    >
      {/* Header: Icon + Title */}
      <div className="flex items-start gap-3 w-full mb-3">
        <div className="shrink-0">
          {skill.icon_url ? (
            <img src={skill.icon_url} alt="" className="w-10 h-10 rounded-lg" />
          ) : (
            <div className="w-10 h-10 rounded-lg bg-accent/10 flex items-center justify-center">
              <Zap className="h-5 w-5 text-accent" />
            </div>
          )}
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2">
            <span className="font-medium text-sm line-clamp-1">{skill.name}</span>
            {isInstalled && (
              <CheckCircle2 className="h-3.5 w-3.5 text-green-500 shrink-0" />
            )}
          </div>
          {skill.is_official && (
            <span className="inline-block mt-1 px-1.5 py-0.5 text-[10px] bg-blue-500/10 text-blue-500 rounded">
              {t('官方')}
            </span>
          )}
        </div>
      </div>

      {/* Description */}
      <p className="text-xs text-muted-foreground line-clamp-2 mb-3 min-h-[2.5em]">
        {skill.description || t('暂无描述')}
      </p>

      {/* Footer: Meta info */}
      <div className="flex items-center gap-3 text-xs text-muted-foreground mt-auto w-full">
        {skill.author_name && (
          <span className="truncate max-w-[100px]">{skill.author_name}</span>
        )}
        <span className="flex items-center gap-1">
          <Download className="h-3 w-3" />
          {skill.download_count}
        </span>
        <span className="ml-auto">
          {skill.pricing_type === 'free' ? (
            <span className="text-green-500">{t('免费')}</span>
          ) : (
            <span className="text-amber-500">¥{skill.price}</span>
          )}
        </span>
      </div>
    </button>
  )
}

interface AppCardProps {
  app: MarketplaceApp
  onClick: () => void
}

function AppCard({ app, onClick }: AppCardProps) {
  const t = useT()

  return (
    <button
      type="button"
      onClick={onClick}
      className="group flex flex-col items-start rounded-xl border border-border/60 bg-background/40 p-4 text-left shadow-[0_0_0_1px_rgba(15,23,42,0.02)] hover:border-border hover:bg-foreground/[0.02] transition-colors cursor-pointer"
    >
      {/* Header: Icon + Title */}
      <div className="flex items-start gap-3 w-full mb-3">
        <div className="shrink-0">
          {app.icon_url ? (
            <img src={app.icon_url} alt="" className="w-10 h-10 rounded-lg" />
          ) : (
            <div className="w-10 h-10 rounded-lg bg-purple-500/10 flex items-center justify-center">
              <Package className="h-5 w-5 text-purple-500" />
            </div>
          )}
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2">
            <span className="font-medium text-sm line-clamp-1">{app.name}</span>
          </div>
          {app.is_official && (
            <span className="inline-block mt-1 px-1.5 py-0.5 text-[10px] bg-blue-500/10 text-blue-500 rounded">
              {t('官方')}
            </span>
          )}
        </div>
      </div>

      {/* Description */}
      <p className="text-xs text-muted-foreground line-clamp-2 mb-3 min-h-[2.5em]">
        {app.description || t('暂无描述')}
      </p>

      {/* Footer: Meta info */}
      <div className="flex items-center gap-3 text-xs text-muted-foreground mt-auto w-full">
        {app.author_name && (
          <span className="truncate max-w-[100px]">{app.author_name}</span>
        )}
        <span className="flex items-center gap-1">
          <Download className="h-3 w-3" />
          {app.download_count}
        </span>
        <span className="ml-auto">
          {app.pricing_type === 'free' ? (
            <span className="text-green-500">{t('免费')}</span>
          ) : app.pricing_type === 'subscription' ? (
            <span className="text-amber-500">{t('订阅')}</span>
          ) : (
            <span className="text-amber-500">¥{app.price}</span>
          )}
        </span>
      </div>
    </button>
  )
}
