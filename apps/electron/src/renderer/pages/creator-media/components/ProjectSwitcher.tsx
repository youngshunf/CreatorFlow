import { useT } from '@/context/LocaleContext'
import {
  Select, SelectTrigger, SelectValue, SelectContent, SelectItem,
} from '@/components/ui/select'
import type { Project } from '@sprouty-ai/shared/db/types'

const PLATFORM_LABELS: Record<string, string> = {
  xiaohongshu: '小红书',
  douyin: '抖音',
  bilibili: 'B站',
  wechat: '微信',
  zhihu: '知乎',
  weibo: '微博',
  x: 'X',
}

interface ProjectSwitcherProps {
  projects: Project[]
  activeProject: Project | null
  onSwitch: (projectId: string) => void
}

export function ProjectSwitcher({ projects, activeProject, onSwitch }: ProjectSwitcherProps) {
  const t = useT()

  if (projects.length === 0) {
    return (
      <span className="text-xs text-muted-foreground">{t('暂无项目')}</span>
    )
  }

  return (
    <Select value={activeProject?.id || ''} onValueChange={onSwitch}>
      <SelectTrigger className="h-7 w-auto min-w-[140px] max-w-[220px] text-xs">
        <SelectValue placeholder={t('选择项目')} />
      </SelectTrigger>
      <SelectContent>
        {projects.map((p) => (
          <SelectItem key={p.id} value={p.id}>
            <span className="flex items-center gap-1.5">
              <span className="truncate">{p.name}</span>
              {p.platform && (
                <span className="text-muted-foreground text-[10px]">
                  {PLATFORM_LABELS[p.platform] || p.platform}
                </span>
              )}
            </span>
          </SelectItem>
        ))}
      </SelectContent>
    </Select>
  )
}
