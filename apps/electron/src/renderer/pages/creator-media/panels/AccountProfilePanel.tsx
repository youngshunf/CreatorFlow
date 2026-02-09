import { useT } from '@/context/LocaleContext'
import { useCreatorMedia } from '../hooks/useCreatorMedia'

/**
 * 账号画像面板 — 在侧边区域显示当前账号画像摘要
 * 当前作为可折叠侧边区域使用，后续迭代集成到右侧边栏
 */
export function AccountProfilePanel() {
  const t = useT()
  const { activeProject, profile } = useCreatorMedia()

  if (!activeProject) {
    return (
      <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-4 text-center">
        <p className="text-xs text-muted-foreground">{t('请先选择项目')}</p>
      </div>
    )
  }

  if (!profile) {
    return (
      <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-4 text-center">
        <p className="text-xs text-muted-foreground">{t('尚未设置账号画像')}</p>
      </div>
    )
  }

  return (
    <div className="rounded-lg border border-border/60 bg-background/40 px-4 py-3 space-y-3">
      {/* 项目名称 + 平台 */}
      <div className="flex items-center gap-2">
        <span className="text-sm font-medium text-foreground">{activeProject.name}</span>
        <span className="inline-flex rounded-full bg-muted/60 px-2 py-0.5 text-[10px] text-muted-foreground">
          {activeProject.platform}
        </span>
      </div>

      {/* 画像信息 */}
      <div className="space-y-1.5 text-xs">
        {profile.niche && (
          <div className="flex gap-2">
            <span className="text-muted-foreground flex-shrink-0">{t('领域')}:</span>
            <span className="text-foreground">{profile.niche}</span>
          </div>
        )}
        {profile.persona && (
          <div className="flex gap-2">
            <span className="text-muted-foreground flex-shrink-0">{t('人设')}:</span>
            <span className="text-foreground">{profile.persona}</span>
          </div>
        )}
        {profile.tone && (
          <div className="flex gap-2">
            <span className="text-muted-foreground flex-shrink-0">{t('调性')}:</span>
            <span className="text-foreground">{profile.tone}</span>
          </div>
        )}
        {profile.target_audience && (
          <div className="flex gap-2">
            <span className="text-muted-foreground flex-shrink-0">{t('目标受众')}:</span>
            <span className="text-foreground">{profile.target_audience}</span>
          </div>
        )}
      </div>

      {/* 关键词标签云 */}
      {profile.keywords && (
        <div className="flex flex-wrap gap-1.5">
          {profile.keywords.split(',').map((kw: string, i: number) => (
            <span key={i} className="inline-flex rounded-full bg-muted/60 px-2 py-0.5 text-[10px] text-muted-foreground">
              {kw.trim()}
            </span>
          ))}
        </div>
      )}
    </div>
  )
}
