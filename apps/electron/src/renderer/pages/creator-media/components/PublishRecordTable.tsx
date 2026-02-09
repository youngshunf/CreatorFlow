import { useT } from '@/context/LocaleContext'
import type { PublishRecord } from '@sprouty-ai/shared/db/types'

/** 发布状态颜色 */
const PUBLISH_STATUS_COLORS: Record<string, string> = {
  pending: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  publishing: 'bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400',
  success: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  failed: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
  deleted: 'bg-gray-100 text-gray-700 dark:bg-gray-900/30 dark:text-gray-400',
}

const PUBLISH_STATUS_LABELS: Record<string, string> = {
  pending: '待发布',
  publishing: '发布中',
  success: '成功',
  failed: '失败',
  deleted: '已删除',
}

const PLATFORM_LABELS: Record<string, string> = {
  xiaohongshu: '小红书',
  douyin: '抖音',
  bilibili: 'B站',
  wechat: '微信',
  zhihu: '知乎',
  weibo: '微博',
  x: 'X',
}

interface PublishRecordTableProps {
  records: PublishRecord[]
  maxItems?: number
}

/**
 * 发布记录排行表格 — 按互动量排序
 */
export function PublishRecordTable({ records, maxItems = 20 }: PublishRecordTableProps) {
  const t = useT()

  // 按互动量排序（views + likes + comments + favorites + shares）
  const sorted = [...records]
    .sort((a, b) => {
      const scoreA = (a.views || 0) + (a.likes || 0) * 3 + (a.comments || 0) * 5 + (a.favorites || 0) * 2 + (a.shares || 0) * 4
      const scoreB = (b.views || 0) + (b.likes || 0) * 3 + (b.comments || 0) * 5 + (b.favorites || 0) * 2 + (b.shares || 0) * 4
      return scoreB - scoreA
    })
    .slice(0, maxItems)

  if (sorted.length === 0) {
    return (
      <div className="flex items-center justify-center py-8 text-muted-foreground">
        <p className="text-sm">{t('暂无发布记录')}</p>
      </div>
    )
  }

  return (
    <div className="overflow-hidden rounded-lg border border-border/60">
      <table className="w-full text-sm">
        <thead>
          <tr className="border-b border-border/40 bg-muted/30">
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('平台')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('状态')}</th>
            <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('阅读')}</th>
            <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('点赞')}</th>
            <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('评论')}</th>
            <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('收藏')}</th>
            <th className="px-4 py-2.5 text-right font-medium text-muted-foreground">{t('分享')}</th>
            <th className="px-4 py-2.5 text-left font-medium text-muted-foreground">{t('发布时间')}</th>
          </tr>
        </thead>
        <tbody>
          {sorted.map((r) => (
            <tr key={r.id} className="border-b border-border/20 last:border-0 hover:bg-muted/20 transition-colors">
              <td className="px-4 py-2.5 text-foreground">
                {t(PLATFORM_LABELS[r.platform] || r.platform)}
              </td>
              <td className="px-4 py-2.5">
                <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-xs font-medium ${PUBLISH_STATUS_COLORS[r.status] || ''}`}>
                  {t(PUBLISH_STATUS_LABELS[r.status] || r.status)}
                </span>
              </td>
              <td className="px-4 py-2.5 text-right text-foreground tabular-nums">{(r.views || 0).toLocaleString()}</td>
              <td className="px-4 py-2.5 text-right text-foreground tabular-nums">{(r.likes || 0).toLocaleString()}</td>
              <td className="px-4 py-2.5 text-right text-foreground tabular-nums">{(r.comments || 0).toLocaleString()}</td>
              <td className="px-4 py-2.5 text-right text-foreground tabular-nums">{(r.favorites || 0).toLocaleString()}</td>
              <td className="px-4 py-2.5 text-right text-foreground tabular-nums">{(r.shares || 0).toLocaleString()}</td>
              <td className="px-4 py-2.5 text-muted-foreground text-xs">
                {r.published_at ? new Date(r.published_at).toLocaleDateString('zh-CN') : '-'}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
