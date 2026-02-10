/**
 * MainContentPanel - Right panel component for displaying content
 *
 * Renders content based on the unified NavigationState:
 * - Chats navigator: ChatPage for selected session, or empty state
 * - Sources navigator: SourceInfoPage for selected source, or empty state
 * - Settings navigator: Settings, Preferences, or Shortcuts page
 *
 * The NavigationState is the single source of truth for what to display.
 *
 * In focused mode (single window), wraps content with StoplightProvider
 * so PanelHeader components automatically compensate for macOS traffic lights.
 */

import * as React from 'react'
import { Panel } from './Panel'
import { cn } from '@/lib/utils'
import { useAppShellContext, useActiveWorkspace } from '@/context/AppShellContext'
import { StoplightProvider } from '@/context/StoplightContext'
import { useT } from '@/context/LocaleContext'
import {
  useNavigationState,
  useNavigation,
  isChatsNavigation,
  isSourcesNavigation,
  isSettingsNavigation,
  isSkillsNavigation,
  isFilesNavigation,
  routes,
} from '@/contexts/NavigationContext'
import { isMarketplaceNavigation, isAppViewNavigation, isVideoNavigation } from '../../../shared/types'
import { APP_VIEW_REGISTRY } from '../../pages/creator-media/registry'
import { UserProfilePage, UserProfileEditPage, AppSettingsPage, AppearanceSettingsPage, InputSettingsPage, WorkspaceSettingsPage, PermissionsSettingsPage, LabelsSettingsPage, PreferencesPage, ShortcutsPage, SourceInfoPage, ChatPage, SubscriptionSettingsPage, SourcesSettingsPage, SkillsSettingsPage } from '@/pages'
import SkillInfoPage from '@/pages/SkillInfoPage'
import { MarketplacePage } from '@/pages/MarketplacePage'
import { FileManager } from '@/components/file-manager'
import { VideoEditor } from '@/components/video'
import { SkillAvatar } from '@/components/ui/skill-avatar'
import { toast } from 'sonner'
import type { LoadedSkill } from '../../../shared/types'

export interface MainContentPanelProps {
  /** Whether the app is in focused mode (single chat, no sidebar) */
  isFocusedMode?: boolean
  /** Optional className for the container */
  className?: string
  /** Callback when user clicks "Use" on a marketplace app */
  onUseMarketplaceApp?: (appId: string, appName: string) => void
}

export function MainContentPanel({
  isFocusedMode = false,
  className,
  onUseMarketplaceApp,
}: MainContentPanelProps) {
  const t = useT()
  const navState = useNavigationState()
  const { navigate } = useNavigation()
  const { activeWorkspaceId, skills = [], onCreateSession, onRenameSession, onSendMessage } = useAppShellContext()
  const activeWorkspace = useActiveWorkspace()

  // Wrap content with StoplightProvider so PanelHeaders auto-compensate in focused mode
  const wrapWithStoplight = (content: React.ReactNode) => (
    <StoplightProvider value={isFocusedMode}>
      {content}
    </StoplightProvider>
  )

  const handleSkillStartChat = React.useCallback(
    async (skill: LoadedSkill) => {
      if (!activeWorkspaceId) return

      try {
        const session = await onCreateSession(activeWorkspaceId)

        const name = skill.metadata.name?.trim()
        if (name) {
          onRenameSession(session.id, name)
        }

        // Build initial prompt that activates the skill via mention syntax
        const initialMessage = activeWorkspaceId
          ? `[skill:${activeWorkspaceId}:${skill.slug}] ${t('请根据该技能的说明开始处理当前工作。')}`
          : `${skill.metadata.name}\n\n${t('请根据该技能的说明开始处理当前工作。')}`

        // Auto-send the initial message with the selected skill as context
        onSendMessage(session.id, initialMessage, undefined, [skill.slug])

        // Navigate to the new chat session
        navigate(routes.view.allChats(session.id))
      } catch (error) {
        console.error('[MainContentPanel] Failed to start chat with skill:', skill.slug, error)
        toast.error(t('创建会话失败'), {
          description: error instanceof Error ? error.message : String(error),
        })
      }
    },
    [activeWorkspaceId, onCreateSession, onRenameSession, onSendMessage, navigate, t]
  )

  // Settings navigator - always has content (subpage determines which page)
  if (isSettingsNavigation(navState)) {
    switch (navState.subpage) {
      case 'user-profile':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <UserProfilePage />
          </Panel>
        )
      case 'user-profile-edit':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <UserProfileEditPage />
          </Panel>
        )
      case 'appearance':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <AppearanceSettingsPage />
          </Panel>
        )
      case 'input':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <InputSettingsPage />
          </Panel>
        )
      case 'workspace':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <WorkspaceSettingsPage />
          </Panel>
        )
      case 'sources':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <SourcesSettingsPage />
          </Panel>
        )
      case 'skills':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <SkillsSettingsPage />
          </Panel>
        )
      case 'permissions':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <PermissionsSettingsPage />
          </Panel>
        )
      case 'labels':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <LabelsSettingsPage />
          </Panel>
        )
      case 'shortcuts':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <ShortcutsPage />
          </Panel>
        )
      case 'preferences':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <PreferencesPage />
          </Panel>
        )
      case 'subscription':
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <SubscriptionSettingsPage />
          </Panel>
        )
      case 'app':
      default:
        return wrapWithStoplight(
          <Panel variant="grow" className={className}>
            <AppSettingsPage />
          </Panel>
        )
    }
  }

  // Sources navigator - show source info or empty state
  if (isSourcesNavigation(navState)) {
    if (navState.details) {
      return wrapWithStoplight(
        <Panel variant="grow" className={className}>
          <SourceInfoPage
            sourceSlug={navState.details.sourceSlug}
            workspaceId={activeWorkspaceId || ''}
          />
        </Panel>
      )
    }
    // No source selected - empty state
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <div className="flex items-center justify-center h-full text-muted-foreground">
          <p className="text-sm">{t('暂无配置的数据源')}</p>
        </div>
      </Panel>
    )
  }

  // Skills navigator - show skill info or empty state
  if (isSkillsNavigation(navState)) {
    if (navState.details) {
      return wrapWithStoplight(
        <Panel variant="grow" className={className}>
          <SkillInfoPage
            skillSlug={navState.details.skillSlug}
            workspaceId={activeWorkspaceId || ''}
          />
        </Panel>
      )
    }
    // No skill selected - empty state
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <div className="flex items-center justify-center h-full text-muted-foreground">
          <p className="text-sm">{t('暂无配置的技能')}</p>
        </div>
      </Panel>
    )
  }

  // Files navigator - show file manager
  if (isFilesNavigation(navState)) {
    // Use workspace rootPath as boundary - users cannot navigate above it
    const workspaceRoot = activeWorkspace?.rootPath
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <FileManager
          initialPath={navState.path || workspaceRoot}
          rootPath={workspaceRoot}
          className="h-full"
        />
      </Panel>
    )
  }

  // Marketplace navigator - show full marketplace page with card grid
  if (isMarketplaceNavigation(navState)) {
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <MarketplacePage
          filter={navState.filter}
          onFilterChange={(filter) => {
            // Update filter via navigation
            if (filter.kind === 'all') {
              navigate(routes.view.marketplace())
            } else if (filter.kind === 'skills') {
              navigate(routes.view.marketplaceSkills())
            } else if (filter.kind === 'apps') {
              navigate(routes.view.marketplaceApps())
            } else if (filter.kind === 'category' && filter.categoryId) {
              navigate(routes.view.marketplace({ filter: 'category', categoryId: filter.categoryId }))
            }
          }}
          workspaceId={activeWorkspaceId || undefined}
          onUseMarketplaceApp={onUseMarketplaceApp}
        />
      </Panel>
    )
  }

  // APP 视图 — 从注册表动态加载视图组件
  if (isAppViewNavigation(navState)) {
    const appViews = APP_VIEW_REGISTRY[navState.appId]
    const ViewComponent = appViews?.[navState.viewId]
    if (ViewComponent) {
      return wrapWithStoplight(
        <Panel variant="grow" className={className}>
          <React.Suspense fallback={<div className="flex items-center justify-center h-full"><div className="text-sm text-muted-foreground">{t('加载中...')}</div></div>}>
            <ViewComponent />
          </React.Suspense>
        </Panel>
      )
    }
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <div className="flex items-center justify-center h-full text-muted-foreground">
          <p className="text-sm">{t('未知视图')}</p>
        </div>
      </Panel>
    )
  }

  // Video navigator - show VideoEditor
  if (isVideoNavigation(navState)) {
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        <VideoEditor
          workspaceId={activeWorkspaceId || ''}
          className="h-full"
        />
      </Panel>
    )
  }

  // Chats navigator - show chat or empty state
  if (isChatsNavigation(navState)) {
    if (navState.details) {
      return wrapWithStoplight(
        <Panel variant="grow" className={className}>
          <ChatPage sessionId={navState.details.sessionId} />
        </Panel>
      )
    }

    const hasSkills = skills.length > 0

    // No session selected - show skills grid home if skills exist, otherwise fallback empty state
    return wrapWithStoplight(
      <Panel variant="grow" className={className}>
        {hasSkills ? (
          <div className="flex flex-col h-full">
            <div className="px-8 pt-6 pb-4 border-b border-border/40">
              <h1 className="text-base font-semibold text-foreground">
                {t('选择一个技能开始对话')}
              </h1>
              <p className="mt-1 text-xs text-muted-foreground">
                {t('根据当前工作区的技能快速启动一个新的对话。')}
              </p>
            </div>
            <div className="flex-1 overflow-auto px-6 py-6">
              <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 xl:grid-cols-3">
                {skills.map((skill) => (
                  <button
                    key={skill.slug}
                    type="button"
                    onClick={() => handleSkillStartChat(skill)}
                    className="group flex flex-col items-start gap-3 rounded-[12px] border border-border/60 bg-background/40 px-4 py-3 text-left shadow-[0_0_0_1px_rgba(15,23,42,0.02)] hover:border-border hover:bg-foreground/[0.02] transition-colors"
                  >
                    <div className="flex items-center gap-3 w-full">
                      <div className="shrink-0">
                        <SkillAvatar
                          skill={skill}
                          size="md"
                          workspaceId={activeWorkspaceId || undefined}
                        />
                      </div>
                      <div className="flex flex-col min-w-0">
                        <div className="text-sm font-medium text-foreground line-clamp-2">
                          {skill.metadata.name}
                        </div>
                        {skill.metadata.description && (
                          <div className="mt-0.5 text-xs text-muted-foreground line-clamp-2">
                            {skill.metadata.description}
                          </div>
                        )}
                      </div>
                    </div>
                    <div className="mt-1 text-[11px] text-primary opacity-0 group-hover:opacity-100 transition-opacity">
                      {t('使用该技能新建一个对话')}
                    </div>
                  </button>
                ))}
              </div>
            </div>
          </div>
        ) : (
          <div className="flex items-center justify-center h-full text-muted-foreground">
            <p className="text-sm">
              {navState.filter.kind === 'flagged'
                ? t('暂无标记的对话')
                : t('暂无对话')}
            </p>
          </div>
        )}
      </Panel>
    )
  }

  // Fallback (should not happen with proper NavigationState)
  return wrapWithStoplight(
    <Panel variant="grow" className={className}>
      <div className="flex items-center justify-center h-full text-muted-foreground">
        <p className="text-sm">{t('选择一个对话开始')}</p>
      </div>
    </Panel>
  )
}
