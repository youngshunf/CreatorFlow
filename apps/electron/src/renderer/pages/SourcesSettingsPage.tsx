/**
 * SourcesSettingsPage - 设置页面中的数据源管理
 *
 * 在设置导航器中显示数据源列表，复用 SourcesListPanel 组件
 */

import * as React from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { SourcesListPanel } from '@/components/app-shell/SourcesListPanel'
import { useAppShellContext, useActiveWorkspace } from '@/context/AppShellContext'
import { useNavigation, routes } from '@/contexts/NavigationContext'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import type { LoadedSource } from '../../shared/types'

export function SourcesSettingsPage() {
  const t = useT()
  const { navigate } = useNavigation()
  const activeWorkspace = useActiveWorkspace()
  const {
    enabledSources = [],
    activeWorkspaceId,
  } = useAppShellContext()

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

  const handleSourceClick = React.useCallback((source: LoadedSource) => {
    navigate(routes.view.sources({ sourceSlug: source.config.slug }))
  }, [navigate])

  return (
    <div className="flex flex-col h-full">
      <PanelHeader
        title={t('数据源')}
      />
      <div className="flex-1 overflow-hidden">
        <SourcesListPanel
          sources={enabledSources}
          sourceFilter={null}
          workspaceRootPath={activeWorkspace?.rootPath}
          onDeleteSource={handleDeleteSource}
          onSourceClick={handleSourceClick}
        />
      </div>
    </div>
  )
}
