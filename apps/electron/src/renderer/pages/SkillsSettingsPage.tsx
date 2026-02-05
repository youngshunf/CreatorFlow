/**
 * SkillsSettingsPage - 设置页面中的技能管理
 *
 * 在设置导航器中显示技能列表，复用 SkillsListPanel 组件
 */

import * as React from 'react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { SkillsListPanel } from '@/components/app-shell/SkillsListPanel'
import { useAppShellContext, useActiveWorkspace } from '@/context/AppShellContext'
import { useNavigation, routes } from '@/contexts/NavigationContext'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import type { LoadedSkill } from '../../shared/types'

export function SkillsSettingsPage() {
  const t = useT()
  const { navigate } = useNavigation()
  const activeWorkspace = useActiveWorkspace()
  const {
    skills = [],
    activeWorkspaceId,
  } = useAppShellContext()

  const handleDeleteSkill = React.useCallback(async (skillSlug: string) => {
    if (!activeWorkspace) return
    try {
      await window.electronAPI.deleteSkill(activeWorkspace.id, skillSlug)
      toast.success(t('已删除技能'))
    } catch (error) {
      console.error('[SkillsSettingsPage] Failed to delete skill:', error)
      toast.error(t('删除技能失败'))
    }
  }, [activeWorkspace, t])

  const handleSkillClick = React.useCallback((skill: LoadedSkill) => {
    navigate(routes.view.skills(skill.slug))
  }, [navigate])

  return (
    <div className="flex flex-col h-full">
      <PanelHeader
        title={t('技能')}
      />
      <div className="flex-1 overflow-hidden">
        <SkillsListPanel
          skills={skills}
          onDeleteSkill={handleDeleteSkill}
          onSkillClick={handleSkillClick}
          workspaceId={activeWorkspaceId || ''}
          workspaceRootPath={activeWorkspace?.rootPath}
        />
      </div>
    </div>
  )
}
