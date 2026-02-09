/**
 * SourcesSettingsPage - 设置页面中的数据源管理
 *
 * 使用卡片网格布局展示数据源，点击卡片在新窗口打开详情
 */

import * as React from 'react'
import { Plus, DatabaseZap, Folder, Trash2 } from 'lucide-react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Empty, EmptyMedia, EmptyTitle, EmptyDescription } from '@/components/ui/empty'
import { SourceAvatar } from '@/components/ui/source-avatar'
import { EditPopover, getEditConfig } from '@/components/ui/EditPopover'
import { useAppShellContext, useActiveWorkspace } from '@/context/AppShellContext'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import type { LoadedSource } from '../../shared/types'

export function SourcesSettingsPage() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const { enabledSources = [] } = useAppShellContext()

  const handleSourceClick = React.useCallback((source: LoadedSource) => {
    // 使用 deep link 在新窗口打开详情
    window.electronAPI.openUrl(`sproutyai://sources/source/${source.config.slug}?window=focused`)
  }, [])

  const handleDeleteSource = React.useCallback(async (sourceSlug: string) => {
    if (!activeWorkspace) return
    try {
      await window.electronAPI.deleteSource(activeWorkspace.id, sourceSlug)
      toast.success(t('已删除数据源'))
    } catch (error) {
      console.error('[SourcesSettingsPage] Failed to delete source:', error)
      toast.error(t('删除数据源失败'))
    }
  }, [activeWorkspace, t])

  return (
    <div className="flex flex-col h-full">
      <PanelHeader
        title={t('数据源')}
        actions={
          activeWorkspace && (
            <EditPopover
              side="bottom"
              align="start"
              width={320}
              trigger={
                <button className="inline-flex items-center gap-1.5 h-8 px-3 text-xs font-medium rounded-md bg-accent text-white hover:bg-accent/90 transition-colors">
                  <Plus className="h-3.5 w-3.5" />
                  {t('添加数据源')}
                </button>
              }
              {...getEditConfig('add-source', activeWorkspace.rootPath)}
            />
          )
        }
      />

      <ScrollArea className="flex-1">
        <div className="p-6">
          {enabledSources.length === 0 ? (
            <Empty>
              <EmptyMedia>
                <DatabaseZap className="h-12 w-12 text-muted-foreground/50" />
              </EmptyMedia>
              <EmptyTitle>{t('未配置数据源')}</EmptyTitle>
              <EmptyDescription>
                {t('添加数据源以连接外部服务和数据')}
              </EmptyDescription>
            </Empty>
          ) : (
            <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 lg:grid-cols-3">
              {enabledSources.map((source) => (
                <SourceCard
                  key={source.config.slug}
                  source={source}
                  onClick={() => handleSourceClick(source)}
                  onDelete={() => handleDeleteSource(source.config.slug)}
                />
              ))}
            </div>
          )}
        </div>
      </ScrollArea>
    </div>
  )
}

interface SourceCardProps {
  source: LoadedSource
  onClick: () => void
  onDelete: () => void
}

function SourceCard({ source, onClick, onDelete }: SourceCardProps) {
  const t = useT()
  const { config } = source

  // 获取类型标签
  const getTypeLabel = () => {
    switch (config.type) {
      case 'mcp': return { label: 'MCP', color: 'bg-blue-500/10 text-blue-500' }
      case 'api': return { label: 'API', color: 'bg-purple-500/10 text-purple-500' }
      case 'local': return { label: t('本地'), color: 'bg-green-500/10 text-green-500' }
      case 'gmail': return { label: 'Gmail', color: 'bg-red-500/10 text-red-500' }
      default: return { label: config.type, color: 'bg-gray-500/10 text-gray-500' }
    }
  }

  const typeInfo = getTypeLabel()

  return (
    <button
      onClick={onClick}
      className="group flex flex-col items-start rounded-xl border border-border/60
        bg-background/40 p-4 text-left shadow-[0_0_0_1px_rgba(15,23,42,0.02)]
        hover:border-border hover:bg-foreground/[0.02] transition-colors cursor-pointer"
    >
      {/* 头部: 图标 + 标题 + 类型徽章 */}
      <div className="flex items-start gap-3 w-full mb-3">
        <div className="shrink-0">
          <SourceAvatar
            source={source}
            size="lg"
            className="w-10 h-10 rounded-lg"
          />
        </div>
        <div className="flex-1 min-w-0">
          <span className="font-medium text-sm line-clamp-1">{config.name}</span>
          <div className="flex items-center gap-1.5 mt-0.5">
            <span className={`px-1.5 py-0.5 text-[10px] rounded font-medium ${typeInfo.color}`}>
              {typeInfo.label}
            </span>
          </div>
        </div>
      </div>

      {/* 描述 */}
      <p className="text-xs text-muted-foreground line-clamp-2 mb-3 min-h-[2.5em]">
        {config.tagline || config.provider || t('无描述')}
      </p>

      {/* 底部元信息 */}
      <div className="flex items-center gap-3 text-xs text-muted-foreground mt-auto w-full">
        <span className="flex items-center gap-1 truncate">
          <Folder className="h-3 w-3 shrink-0" />
          <span className="truncate">{config.slug}</span>
        </span>
        <button
          onClick={(e) => {
            e.stopPropagation()
            onDelete()
          }}
          className="ml-auto text-destructive hover:text-destructive/80 shrink-0"
        >
          <Trash2 className="h-3 w-3" />
        </button>
      </div>
    </button>
  )
}
