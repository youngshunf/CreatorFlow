import { useState, useMemo } from 'react'
import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { StatCards } from './components/StatCard'
import { ContentTable } from './components/ContentTable'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { CreateProjectDialog } from './components/CreateProjectDialog'
import { ProjectSettingsDialog } from './components/ProjectSettingsDialog'
import { ProfileEditDialog } from './components/ProfileEditDialog'
import { CreateContentDialog } from './components/CreateContentDialog'
import { VersionHistoryDialog } from './components/VersionHistoryDialog'
import { HotTopicsPanel } from './components/HotTopicsPanel'
import { TopicRecommendPanel } from './components/TopicRecommendPanel'
import type { Content } from '@sprouty-ai/shared/db/types'
import { POSTING_FREQUENCY_LIST } from '@sprouty-ai/shared/db/types'

/**
 * 创作面板 — 项目仪表盘视图
 * 展示项目统计、内容列表、账号画像摘要
 */
export default function ProjectDashboard() {
  const t = useT()
  const {
    projects, activeProject, contents, profile, stats, loading,
    switchProject, createProject, updateProject, deleteProject,
    upsertProfile, createContent, deleteContent, updateContentStatus,
    listContentVersions, rollbackContentVersion,
  } = useCreatorMedia()

  const [showCreateProject, setShowCreateProject] = useState(false)
  const [showProjectSettings, setShowProjectSettings] = useState(false)
  const [showProfileEdit, setShowProfileEdit] = useState(false)
  const [showCreateContent, setShowCreateContent] = useState(false)
  const [versionHistoryContent, setVersionHistoryContent] = useState<Content | null>(null)

  // 视频相关统计
  const videoStats = useMemo(() => {
    const videoContents = contents.filter(c => c.content_type === 'video' || c.content_type === 'short-video')
    const creating = videoContents.filter(c => c.status === 'creating').length
    let completed = 0
    for (const c of videoContents) {
      if (c.metadata) {
        try {
          const meta = JSON.parse(c.metadata)
          if (meta.video_render_status === 'completed') completed++
        } catch { /* ignore */ }
      }
    }
    return { creating, completed }
  }, [contents])

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  // 无项目时显示引导创建界面
  if (projects.length === 0) {
    return (
      <div className="flex flex-col h-full">
        <div className="flex-1 flex items-center justify-center">
          <div className="text-center space-y-4 max-w-sm">
            <div className="mx-auto w-12 h-12 rounded-full bg-muted/60 flex items-center justify-center">
              <svg className="w-6 h-6 text-muted-foreground" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
              </svg>
            </div>
            <div>
              <h2 className="text-base font-semibold text-foreground">{t('开始你的创作之旅')}</h2>
              <p className="mt-1 text-sm text-muted-foreground">{t('创建第一个项目，开始管理你的自媒体内容')}</p>
            </div>
            <button
              type="button"
              onClick={() => setShowCreateProject(true)}
              className="inline-flex items-center gap-1.5 rounded-md bg-foreground px-4 py-2 text-sm font-medium text-background hover:bg-foreground/90 transition-colors"
            >
              <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
              </svg>
              {t('新建项目')}
            </button>
          </div>
        </div>

        <CreateProjectDialog
          open={showCreateProject}
          onOpenChange={setShowCreateProject}
          onCreateProject={createProject}
          onUpsertProfile={async (_projectId, data) => upsertProfile(data)}
        />
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* 头部 — relative z-panel 确保在 titlebar 拖拽区域之上 */}
      <div className="relative z-panel flex items-center justify-between px-6 py-4 border-b border-border/40">
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('创作面板')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('项目概览与内容管理')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <button
            type="button"
            onClick={() => setShowCreateProject(true)}
            className="inline-flex items-center gap-1 rounded-md border border-border/60 bg-background px-2.5 py-1 text-xs text-foreground hover:bg-muted/40 transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
            </svg>
            {t('新建项目')}
          </button>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
          {activeProject && (
            <button
              type="button"
              onClick={() => setShowProjectSettings(true)}
              className="inline-flex items-center justify-center rounded-md w-7 h-7 text-muted-foreground hover:text-foreground hover:bg-muted/40 transition-colors"
              title={t('项目设置')}
            >
              <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9.594 3.94c.09-.542.56-.94 1.11-.94h2.593c.55 0 1.02.398 1.11.94l.213 1.281c.063.374.313.686.645.87.074.04.147.083.22.127.325.196.72.257 1.075.124l1.217-.456a1.125 1.125 0 0 1 1.37.49l1.296 2.247a1.125 1.125 0 0 1-.26 1.431l-1.003.827c-.293.241-.438.613-.43.992a7.723 7.723 0 0 1 0 .255c-.008.378.137.75.43.991l1.004.827c.424.35.534.955.26 1.43l-1.298 2.247a1.125 1.125 0 0 1-1.369.491l-1.217-.456c-.355-.133-.75-.072-1.076.124a6.47 6.47 0 0 1-.22.128c-.331.183-.581.495-.644.869l-.213 1.281c-.09.543-.56.94-1.11.94h-2.594c-.55 0-1.019-.398-1.11-.94l-.213-1.281c-.062-.374-.312-.686-.644-.87a6.52 6.52 0 0 1-.22-.127c-.325-.196-.72-.257-1.076-.124l-1.217.456a1.125 1.125 0 0 1-1.369-.49l-1.297-2.247a1.125 1.125 0 0 1 .26-1.431l1.004-.827c.292-.24.437-.613.43-.991a6.932 6.932 0 0 1 0-.255c.007-.38-.138-.751-.43-.992l-1.004-.827a1.125 1.125 0 0 1-.26-1.43l1.297-2.247a1.125 1.125 0 0 1 1.37-.491l1.216.456c.356.133.751.072 1.076-.124.072-.044.146-.086.22-.128.332-.183.582-.495.644-.869l.214-1.28Z" />
                <path strokeLinecap="round" strokeLinejoin="round" d="M15 12a3 3 0 1 1-6 0 3 3 0 0 1 6 0Z" />
              </svg>
            </button>
          )}
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-6">
        {/* 统计卡片 */}
        <StatCards stats={stats} videoCreating={videoStats.creating} videoCompleted={videoStats.completed} />

        {/* 最近内容 + 账号画像：大屏左右，小屏上下 */}
        <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">
          {/* 最近内容 */}
          <div className="min-w-0">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-sm font-medium text-foreground">{t('最近内容')}</h2>
              <button
                type="button"
                onClick={() => setShowCreateContent(true)}
                className="inline-flex items-center gap-1 rounded-md text-xs text-muted-foreground hover:text-foreground transition-colors"
              >
                <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
                </svg>
                {t('新建内容')}
              </button>
            </div>
            <ContentTable
              contents={contents}
              maxItems={10}
              onStatusChange={updateContentStatus}
              onDelete={deleteContent}
              onVersionHistory={setVersionHistoryContent}
            />
          </div>

          {/* 账号画像 */}
          <div className="min-w-0">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-sm font-medium text-foreground">{t('账号画像')}</h2>
              {activeProject && (
                <button
                  type="button"
                  onClick={() => setShowProfileEdit(true)}
                  className="inline-flex items-center gap-1 rounded-md text-xs text-muted-foreground hover:text-foreground transition-colors"
                >
                  <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                    <path strokeLinecap="round" strokeLinejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L6.832 19.82a4.5 4.5 0 0 1-1.897 1.13l-2.685.8.8-2.685a4.5 4.5 0 0 1 1.13-1.897L16.863 4.487Z" />
                  </svg>
                  {t('编辑画像')}
                </button>
              )}
            </div>
            {profile ? (
              <div className="rounded-lg border border-border/60 bg-background/40 px-4 py-3 space-y-3">
                {/* 核心字段 */}
                <div className="grid grid-cols-2 gap-x-4 gap-y-2 text-xs">
                  {profile.niche && (
                    <div>
                      <span className="text-muted-foreground">{t('领域')}</span>
                      <p className="text-foreground mt-0.5">{profile.niche}</p>
                    </div>
                  )}
                  {profile.sub_niche && (
                    <div>
                      <span className="text-muted-foreground">{t('细分领域')}</span>
                      <p className="text-foreground mt-0.5">{profile.sub_niche}</p>
                    </div>
                  )}
                  {profile.persona && (
                    <div>
                      <span className="text-muted-foreground">{t('人设')}</span>
                      <p className="text-foreground mt-0.5">{profile.persona}</p>
                    </div>
                  )}
                  {profile.target_audience && (
                    <div>
                      <span className="text-muted-foreground">{t('目标受众')}</span>
                      <p className="text-foreground mt-0.5">{profile.target_audience}</p>
                    </div>
                  )}
                  {profile.tone && (
                    <div>
                      <span className="text-muted-foreground">{t('调性')}</span>
                      <p className="text-foreground mt-0.5">{profile.tone}</p>
                    </div>
                  )}
                  {profile.posting_frequency && (
                    <div>
                      <span className="text-muted-foreground">{t('发布频率')}</span>
                      <p className="text-foreground mt-0.5">{t(POSTING_FREQUENCY_LIST.find(f => f.id === profile.posting_frequency)?.label ?? profile.posting_frequency)}</p>
                    </div>
                  )}
                  {profile.best_posting_time && (
                    <div>
                      <span className="text-muted-foreground">{t('最佳发布时间')}</span>
                      <p className="text-foreground mt-0.5">{profile.best_posting_time}</p>
                    </div>
                  )}
                </div>

                {/* 简介 */}
                {profile.bio && (
                  <div className="text-xs">
                    <span className="text-muted-foreground">{t('简介')}</span>
                    <p className="text-foreground mt-0.5 leading-relaxed">{profile.bio}</p>
                  </div>
                )}

                {/* 关键词 */}
                {profile.keywords && (
                  <div className="text-xs">
                    <span className="text-muted-foreground">{t('关键词')}</span>
                    <div className="flex flex-wrap gap-1.5 mt-1">
                      {profile.keywords.split(',').map((kw: string, i: number) => (
                        <span key={i} className="inline-flex rounded-full bg-muted/60 px-2 py-0.5 text-[10px] text-muted-foreground">
                          {kw.trim()}
                        </span>
                      ))}
                    </div>
                  </div>
                )}

                {/* 内容支柱 */}
                {profile.content_pillars && (
                  <div className="text-xs">
                    <span className="text-muted-foreground">{t('内容支柱')}</span>
                    <div className="flex flex-wrap gap-1.5 mt-1">
                      {profile.content_pillars.split(',').map((p: string, i: number) => {
                        let weight: string | null = null
                        if (profile.pillar_weights) {
                          try {
                            const weights = JSON.parse(profile.pillar_weights)
                            const w = weights[p.trim()]
                            if (w != null) weight = `${Math.round(w * 100)}%`
                          } catch { /* ignore */ }
                        }
                        return (
                          <span key={i} className="inline-flex items-center rounded-full bg-accent/10 px-2 py-0.5 text-[10px] text-accent-foreground">
                            {p.trim()}
                            {weight && <span className="ml-1 opacity-60">{weight}</span>}
                          </span>
                        )
                      })}
                    </div>
                  </div>
                )}

                {/* 风格参考 */}
                {profile.style_references && (
                  <div className="text-xs">
                    <span className="text-muted-foreground">{t('风格参考')}</span>
                    <p className="text-foreground mt-0.5">{profile.style_references}</p>
                  </div>
                )}

                {/* 禁忌话题 */}
                {profile.taboo_topics && (
                  <div className="text-xs">
                    <span className="text-muted-foreground">{t('禁忌话题')}</span>
                    <p className="text-foreground mt-0.5">{profile.taboo_topics}</p>
                  </div>
                )}
              </div>
            ) : (
              <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-6 text-center">
                <p className="text-sm text-muted-foreground">{t('尚未设置账号画像')}</p>
                {activeProject && (
                  <button
                    type="button"
                    onClick={() => setShowProfileEdit(true)}
                    className="mt-2 inline-flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground transition-colors"
                  >
                    <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
                    </svg>
                    {t('设置画像')}
                  </button>
                )}
              </div>
            )}
          </div>
        </div>

        {/* 热榜 + 选题推荐：大屏左右，小屏上下 */}
        {activeProject && (
          <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">
            <HotTopicsPanel />
            <TopicRecommendPanel projectId={activeProject.id} />
          </div>
        )}
      </div>

      {/* 对话框 */}
      <CreateProjectDialog
        open={showCreateProject}
        onOpenChange={setShowCreateProject}
        onCreateProject={createProject}
        onUpsertProfile={async (_projectId, data) => upsertProfile(data)}
      />

      {activeProject && (
        <>
          <ProjectSettingsDialog
            open={showProjectSettings}
            onOpenChange={setShowProjectSettings}
            project={activeProject}
            profile={profile}
            onUpdate={updateProject}
            onDelete={deleteProject}
          />

          <ProfileEditDialog
            open={showProfileEdit}
            onOpenChange={setShowProfileEdit}
            projectId={activeProject.id}
            profile={profile}
            onSave={upsertProfile}
          />

          <CreateContentDialog
            open={showCreateContent}
            onOpenChange={setShowCreateContent}
            onCreateContent={createContent}
          />

          <VersionHistoryDialog
            open={versionHistoryContent !== null}
            onOpenChange={(open) => { if (!open) setVersionHistoryContent(null) }}
            contentId={versionHistoryContent?.id || ''}
            contentTitle={versionHistoryContent?.title || null}
            onRollback={() => updateContentStatus(versionHistoryContent!.id, versionHistoryContent!.status)}
            listContentVersions={listContentVersions}
            rollbackContentVersion={rollbackContentVersion}
          />
        </>
      )}
    </div>
  )
}
