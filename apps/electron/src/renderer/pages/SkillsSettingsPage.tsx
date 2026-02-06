/**
 * SkillsSettingsPage - 设置页面中的技能管理
 *
 * 使用卡片网格布局展示技能，点击卡片在新窗口打开详情
 */

import * as React from 'react'
import { Plus, Zap, Folder, Trash2 } from 'lucide-react'
import { PanelHeader } from '@/components/app-shell/PanelHeader'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Empty, EmptyMedia, EmptyTitle, EmptyDescription } from '@/components/ui/empty'
import { SkillAvatar } from '@/components/ui/skill-avatar'
import { EditPopover, getEditConfig } from '@/components/ui/EditPopover'
import { useAppShellContext, useActiveWorkspace } from '@/context/AppShellContext'
import { useT } from '@/context/LocaleContext'
import { toast } from 'sonner'
import type { LoadedSkill } from '../../shared/types'

export function SkillsSettingsPage() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const { skills = [], activeWorkspaceId } = useAppShellContext()

  const handleSkillClick = React.useCallback((skill: LoadedSkill) => {
    // 使用 deep link 在新窗口打开详情
    window.electronAPI.openUrl(`craftagents://skills/skill/${skill.slug}?window=focused`)
  }, [])

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

  return (
    <div className="flex flex-col h-full">
      <PanelHeader
        title={t('技能')}
        actions={
          activeWorkspace && (
            <EditPopover
              side="bottom"
              align="start"
              width={320}
              trigger={
                <button className="inline-flex items-center gap-1.5 h-8 px-3 text-xs font-medium rounded-md bg-accent text-white hover:bg-accent/90 transition-colors">
                  <Plus className="h-3.5 w-3.5" />
                  {t('添加技能')}
                </button>
              }
              {...getEditConfig('add-skill', activeWorkspace.rootPath)}
            />
          )
        }
      />

      <ScrollArea className="flex-1">
        <div className="p-6">
          {skills.length === 0 ? (
            <Empty>
              <EmptyMedia>
                <Zap className="h-12 w-12 text-muted-foreground/50" />
              </EmptyMedia>
              <EmptyTitle>{t('未安装技能')}</EmptyTitle>
              <EmptyDescription>
                {t('添加技能以扩展 AI 助手的能力')}
              </EmptyDescription>
            </Empty>
          ) : (
            <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 lg:grid-cols-3">
              {skills.map((skill) => (
                <SkillCard
                  key={skill.slug}
                  skill={skill}
                  workspaceId={activeWorkspaceId}
                  onClick={() => handleSkillClick(skill)}
                  onDelete={() => handleDeleteSkill(skill.slug)}
                />
              ))}
            </div>
          )}
        </div>
      </ScrollArea>
    </div>
  )
}

interface SkillCardProps {
  skill: LoadedSkill
  workspaceId?: string
  onClick: () => void
  onDelete: () => void
}

function SkillCard({ skill, workspaceId, onClick, onDelete }: SkillCardProps) {
  const t = useT()
  const { metadata } = skill

  return (
    <button
      onClick={onClick}
      className="group flex flex-col items-start rounded-xl border border-border/60
        bg-background/40 p-4 text-left shadow-[0_0_0_1px_rgba(15,23,42,0.02)]
        hover:border-border hover:bg-foreground/[0.02] transition-colors cursor-pointer"
    >
      {/* 头部: 图标 + 标题 */}
      <div className="flex items-start gap-3 w-full mb-3">
        <div className="shrink-0">
          <SkillAvatar
            skill={skill}
            workspaceId={workspaceId}
            size="lg"
            className="w-10 h-10 rounded-lg"
          />
        </div>
        <div className="flex-1 min-w-0">
          <span className="font-medium text-sm line-clamp-1">{metadata.name}</span>
          <div className="flex items-center gap-1.5 mt-0.5">
            <span className="px-1.5 py-0.5 text-[10px] bg-amber-500/10 text-amber-600 dark:text-amber-400 rounded font-medium">
              {t('技能')}
            </span>
            {metadata.version && (
              <span className="text-[10px] text-muted-foreground">
                v{metadata.version}
              </span>
            )}
          </div>
        </div>
      </div>

      {/* 描述 */}
      <p className="text-xs text-muted-foreground line-clamp-2 mb-3 min-h-[2.5em]">
        {metadata.description || t('无描述')}
      </p>

      {/* 底部元信息 */}
      <div className="flex items-center gap-3 text-xs text-muted-foreground mt-auto w-full">
        <span className="flex items-center gap-1 truncate">
          <Folder className="h-3 w-3 shrink-0" />
          <span className="truncate">{skill.slug}</span>
        </span>
        {metadata.author && (
          <span className="truncate max-w-[100px]">{metadata.author}</span>
        )}
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
