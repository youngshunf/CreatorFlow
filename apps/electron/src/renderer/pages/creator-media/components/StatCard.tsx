import { useT } from '@/context/LocaleContext'
import type { ContentStats } from '../types'

interface StatCardProps {
  label: string
  value: number
  icon: React.ReactNode
  color?: string
}

function StatCard({ label, value, icon, color = 'text-foreground' }: StatCardProps) {
  return (
    <div className="flex items-center gap-3 rounded-lg border border-border/60 bg-background/40 px-4 py-3">
      <div className={`flex h-9 w-9 items-center justify-center rounded-md bg-muted/50 ${color}`}>
        {icon}
      </div>
      <div>
        <p className="text-xl font-semibold text-foreground">{value}</p>
        <p className="text-xs text-muted-foreground">{label}</p>
      </div>
    </div>
  )
}

interface StatCardsProps {
  stats: ContentStats
}

export function StatCards({ stats }: StatCardsProps) {
  const t = useT()
  return (
    <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
      <StatCard
        label={t('选题')}
        value={stats.idea + stats.researching}
        icon={<svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3"/></svg>}
        color="text-blue-500"
      />
      <StatCard
        label={t('进行中')}
        value={stats.scripting + stats.creating}
        icon={<svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M12 22c5.523 0 10-4.477 10-10S17.523 2 12 2 2 6.477 2 12s4.477 10 10 10z"/><path d="M12 6v6l4 2"/></svg>}
        color="text-amber-500"
      />
      <StatCard
        label={t('待发布')}
        value={stats.reviewing + stats.scheduled}
        icon={<svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M16 21v-2a4 4 0 0 0-4-4H6a4 4 0 0 0-4 4v2"/><circle cx="9" cy="7" r="4"/><path d="M22 21v-2a4 4 0 0 0-3-3.87"/><path d="M16 3.13a4 4 0 0 1 0 7.75"/></svg>}
        color="text-purple-500"
      />
      <StatCard
        label={t('已发布')}
        value={stats.published}
        icon={<svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="m9 12 2 2 4-4"/><path d="M12 22c5.523 0 10-4.477 10-10S17.523 2 12 2 2 6.477 2 12s4.477 10 10 10z"/></svg>}
        color="text-green-500"
      />
    </div>
  )
}
