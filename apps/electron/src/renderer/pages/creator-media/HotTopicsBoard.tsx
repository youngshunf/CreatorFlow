import { useState } from 'react'
import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { PLATFORM_LIST } from '@sprouty-ai/shared/db/types'

/** 平台标签（热点源包含非标准平台） */
const PLATFORMS = [
  { id: 'all', label: '全部' },
  ...PLATFORM_LIST.map(p => ({ id: p.id, label: p.shortLabel })),
  { id: 'baidu', label: '百度' },
]

/**
 * 热点看板 — 实时热点流 + 智能选题推荐
 * 当前为引导式 UI，topic_cache IPC 尚未集成
 */
export default function HotTopicsBoard() {
  const t = useT()
  const {
    projects, activeProject, profile, loading, switchProject,
  } = useCreatorMedia()

  const [activePlatform, setActivePlatform] = useState('all')

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* 头部 */}
      <div className="relative z-panel flex items-center justify-between px-6 py-4 border-b border-border/40">
        <div>
          <h1 className="text-base font-semibold text-foreground">{t('热点看板')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('实时热点流与智能选题')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-6">
        {/* 推荐选题区 */}
        <div>
          <h2 className="text-sm font-medium text-foreground mb-3">{t('推荐选题')}</h2>
          {profile?.keywords ? (
            <div className="rounded-lg border border-border/60 bg-background/40 px-4 py-4">
              <p className="text-xs text-muted-foreground mb-2">
                {t('基于你的账号关键词匹配热点选题')}
              </p>
              <div className="flex flex-wrap gap-1.5 mb-3">
                {profile.keywords.split(',').map((kw: string, i: number) => (
                  <span key={i} className="inline-flex rounded-full bg-blue-100 dark:bg-blue-900/30 px-2 py-0.5 text-[10px] text-blue-700 dark:text-blue-400">
                    {kw.trim()}
                  </span>
                ))}
              </div>
              <div className="flex items-center justify-center py-6 border border-dashed border-border/60 rounded-md">
                <div className="text-center space-y-2">
                  <svg className="mx-auto w-8 h-8 text-muted-foreground/50" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                    <path strokeLinecap="round" strokeLinejoin="round" d="M2.25 18 9 11.25l4.306 4.306a11.95 11.95 0 0 1 5.814-5.518l2.74-1.22m0 0-5.94-2.281m5.94 2.28-2.28 5.941" />
                  </svg>
                  <p className="text-sm text-muted-foreground">{t('在对话中使用"帮我找选题"获取热点推荐')}</p>
                  <p className="text-xs text-muted-foreground/70">{t('将调用 topic-generator 技能生成选题')}</p>
                </div>
              </div>
            </div>
          ) : (
            <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-6 text-center">
              <p className="text-sm text-muted-foreground">{t('请先设置账号画像和关键词')}</p>
              <p className="mt-1 text-xs text-muted-foreground/70">{t('在创作面板中编辑账号画像')}</p>
            </div>
          )}
        </div>

        {/* 实时热榜区 */}
        <div>
          <div className="flex items-center justify-between mb-3">
            <h2 className="text-sm font-medium text-foreground">{t('实时热榜')}</h2>
          </div>

          {/* 平台标签切换 */}
          <div className="flex items-center gap-1 mb-4">
            {PLATFORMS.map((p) => (
              <button
                key={p.id}
                type="button"
                onClick={() => setActivePlatform(p.id)}
                className={`rounded-full px-3 py-1 text-xs font-medium transition-colors ${
                  activePlatform === p.id
                    ? 'bg-foreground text-background'
                    : 'bg-muted/50 text-muted-foreground hover:bg-muted'
                }`}
              >
                {t(p.label)}
              </button>
            ))}
          </div>

          {/* 空状态 — 引导用户获取热点 */}
          <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-8">
            <div className="text-center space-y-3">
              <svg className="mx-auto w-10 h-10 text-muted-foreground/40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M15.362 5.214A8.252 8.252 0 0 1 12 21 8.25 8.25 0 0 1 6.038 7.047 8.287 8.287 0 0 0 9 9.601a8.983 8.983 0 0 1 3.361-6.867 8.21 8.21 0 0 0 3 2.48Z" />
                <path strokeLinecap="round" strokeLinejoin="round" d="M12 18a3.75 3.75 0 0 0 .495-7.468 5.99 5.99 0 0 0-1.925 3.547 5.975 5.975 0 0 1-2.133-1.001A3.75 3.75 0 0 0 12 18Z" />
              </svg>
              <div>
                <p className="text-sm font-medium text-foreground">{t('热点数据尚未加载')}</p>
                <p className="mt-1 text-xs text-muted-foreground">
                  {t('在对话中说"有什么热点"或"帮我找选题"，AI 将通过 NewsNow 和 TrendRadar 获取实时热点')}
                </p>
              </div>
              <div className="flex items-center justify-center gap-2 pt-1">
                <span className="inline-flex items-center rounded-full bg-muted/60 px-2.5 py-0.5 text-[10px] text-muted-foreground">
                  NewsNow MCP
                </span>
                <span className="inline-flex items-center rounded-full bg-muted/60 px-2.5 py-0.5 text-[10px] text-muted-foreground">
                  TrendRadar MCP
                </span>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
