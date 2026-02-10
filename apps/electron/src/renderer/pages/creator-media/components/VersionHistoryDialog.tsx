import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import type { ContentVersion, VersionStage, ChangeSource } from '@sprouty-ai/shared/db/types'

/** 阶段标签映射 */
const STAGE_LABELS: Record<VersionStage, string> = {
  script: '脚本',
  content: '内容',
  adapted: '适配',
}

/** 阶段标签颜色 */
const STAGE_COLORS: Record<VersionStage, string> = {
  script: 'bg-blue-500/15 text-blue-600 dark:text-blue-400',
  content: 'bg-green-500/15 text-green-600 dark:text-green-400',
  adapted: 'bg-purple-500/15 text-purple-600 dark:text-purple-400',
}

/** 来源标签映射 */
const SOURCE_LABELS: Record<ChangeSource, string> = {
  auto: '自动',
  user_edit: '手动',
  rollback: '回滚',
}

interface VersionHistoryDialogProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  contentId: string
  contentTitle: string | null
  onRollback?: () => void
  listContentVersions: (contentId: string) => Promise<ContentVersion[]>
  rollbackContentVersion: (contentId: string, versionNumber: number) => Promise<ContentVersion | null>
}

/**
 * 版本历史对话框 — 时间线布局，支持预览、对比、回滚
 */
export function VersionHistoryDialog({
  open,
  onOpenChange,
  contentId,
  contentTitle,
  onRollback,
  listContentVersions,
  rollbackContentVersion,
}: VersionHistoryDialogProps) {
  const t = useT()

  const [versions, setVersions] = useState<ContentVersion[]>([])
  const [loading, setLoading] = useState(false)
  const [expandedId, setExpandedId] = useState<string | null>(null)
  const [diffTargetId, setDiffTargetId] = useState<string | null>(null)
  const [rollbackConfirm, setRollbackConfirm] = useState<number | null>(null)
  const [rollbackLoading, setRollbackLoading] = useState(false)

  /** 加载版本列表 */
  const loadVersions = useCallback(async () => {
    if (!contentId) return
    setLoading(true)
    try {
      const result = await listContentVersions(contentId)
      setVersions(result)
    } catch {
      setVersions([])
    } finally {
      setLoading(false)
    }
  }, [contentId, listContentVersions])

  useEffect(() => {
    if (open && contentId) {
      loadVersions()
      setExpandedId(null)
      setDiffTargetId(null)
      setRollbackConfirm(null)
    }
  }, [open, contentId, loadVersions])

  /** 切换预览展开 */
  const toggleExpand = (id: string) => {
    setExpandedId(expandedId === id ? null : id)
    setDiffTargetId(null)
  }

  /** 切换对比模式 */
  const toggleDiff = (id: string) => {
    setDiffTargetId(diffTargetId === id ? null : id)
    setExpandedId(null)
  }

  /** 执行回滚 */
  const handleRollback = async (versionNumber: number) => {
    setRollbackLoading(true)
    try {
      await rollbackContentVersion(contentId, versionNumber)
      await loadVersions()
      setRollbackConfirm(null)
      onRollback?.()
    } finally {
      setRollbackLoading(false)
    }
  }

  /** 简单文本对比：逐行比较，标记差异 */
  const renderDiff = (currentSnapshot: string, targetSnapshot: string) => {
    const currentLines = parseSnapshot(currentSnapshot).split('\n')
    const targetLines = parseSnapshot(targetSnapshot).split('\n')
    const maxLen = Math.max(currentLines.length, targetLines.length)
    const diffLines: Array<{ line: string; type: 'same' | 'added' | 'removed' }> = []

    for (let i = 0; i < maxLen; i++) {
      const curr = currentLines[i] ?? ''
      const tgt = targetLines[i] ?? ''
      if (curr === tgt) {
        diffLines.push({ line: curr, type: 'same' })
      } else {
        if (curr) diffLines.push({ line: curr, type: 'removed' })
        if (tgt) diffLines.push({ line: tgt, type: 'added' })
      }
    }

    return (
      <div className="text-xs font-mono whitespace-pre-wrap">
        {diffLines.map((d, i) => (
          <div
            key={i}
            className={
              d.type === 'added'
                ? 'bg-green-500/10 text-green-700 dark:text-green-400'
                : d.type === 'removed'
                  ? 'bg-red-500/10 text-red-700 dark:text-red-400 line-through'
                  : 'text-muted-foreground'
            }
          >
            {d.type === 'added' ? '+ ' : d.type === 'removed' ? '- ' : '  '}
            {d.line || ' '}
          </div>
        ))}
      </div>
    )
  }

  if (!open) return null

  const latestVersion = versions[0]

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center">
      {/* 遮罩 */}
      <div className="absolute inset-0 bg-black/50" onClick={() => onOpenChange(false)} />

      {/* 对话框 */}
      <div className="relative w-full max-w-2xl rounded-lg border border-border bg-background shadow-lg mx-4">
        <div className="px-6 py-4 border-b border-border/40">
          <h2 className="text-base font-semibold text-foreground">
            {t('版本历史')} - {contentTitle || t('未命名内容')}
          </h2>
          <p className="text-xs text-muted-foreground mt-1">
            {t('共 {count} 个版本').replace('{count}', String(versions.length))}
          </p>
        </div>

        <div className="px-6 py-4 max-h-[65vh] overflow-auto">
          {loading ? (
            <div className="flex items-center justify-center py-12 text-sm text-muted-foreground">
              {t('加载中...')}
            </div>
          ) : versions.length === 0 ? (
            <div className="flex items-center justify-center py-12 text-sm text-muted-foreground">
              {t('暂无版本记录')}
            </div>
          ) : (
            <div className="relative">
              {/* 时间线竖线 */}
              <div className="absolute left-3 top-2 bottom-2 w-px bg-border/60" />

              <div className="space-y-1">
                {versions.map((version, index) => {
                  const isLatest = index === 0
                  const isExpanded = expandedId === version.id
                  const isDiffTarget = diffTargetId === version.id
                  const stage = version.stage as VersionStage
                  const source = version.change_source as ChangeSource | null

                  return (
                    <div key={version.id} className="relative pl-8">
                      {/* 时间线圆点 */}
                      <div
                        className={`absolute left-1.5 top-3 w-3 h-3 rounded-full border-2 ${
                          isLatest
                            ? 'bg-foreground border-foreground'
                            : 'bg-background border-border'
                        }`}
                      />

                      <div className="rounded-md border border-border/40 p-3 hover:bg-muted/20 transition-colors">
                        {/* 版本头部 */}
                        <div className="flex items-center justify-between gap-2">
                          <div className="flex items-center gap-2 min-w-0">
                            <span className="text-sm font-medium text-foreground shrink-0">
                              v{version.version_number}
                            </span>
                            <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-medium ${STAGE_COLORS[stage] || 'bg-muted text-muted-foreground'}`}>
                              {t(STAGE_LABELS[stage] || stage)}
                            </span>
                            {source && (
                              <span className="inline-flex items-center rounded-full bg-muted px-2 py-0.5 text-[10px] text-muted-foreground">
                                {t(SOURCE_LABELS[source] || source)}
                              </span>
                            )}
                            {isLatest && (
                              <span className="inline-flex items-center rounded-full bg-foreground/10 px-2 py-0.5 text-[10px] font-medium text-foreground">
                                {t('当前')}
                              </span>
                            )}
                          </div>
                          <span className="text-[10px] text-muted-foreground shrink-0">
                            {formatTime(version.created_at)}
                          </span>
                        </div>

                        {/* 描述 */}
                        {version.change_description && (
                          <p className="text-xs text-muted-foreground mt-1 truncate">
                            {version.change_description}
                          </p>
                        )}

                        {/* 标题 */}
                        {version.title && (
                          <p className="text-xs text-foreground/70 mt-1 truncate">
                            {version.title}
                          </p>
                        )}

                        {/* 操作按钮 */}
                        <div className="flex items-center gap-2 mt-2">
                          <button
                            type="button"
                            onClick={() => toggleExpand(version.id)}
                            className="text-[11px] text-muted-foreground hover:text-foreground transition-colors"
                          >
                            {isExpanded ? t('收起') : t('预览')}
                          </button>
                          {latestVersion && !isLatest && (
                            <>
                              <span className="text-border">|</span>
                              <button
                                type="button"
                                onClick={() => toggleDiff(version.id)}
                                className="text-[11px] text-muted-foreground hover:text-foreground transition-colors"
                              >
                                {isDiffTarget ? t('收起对比') : t('对比当前')}
                              </button>
                              <span className="text-border">|</span>
                              <button
                                type="button"
                                onClick={() => setRollbackConfirm(version.version_number)}
                                className="text-[11px] text-red-500 hover:text-red-600 transition-colors"
                              >
                                {t('回滚到此版本')}
                              </button>
                            </>
                          )}
                        </div>

                        {/* 预览内容 */}
                        {isExpanded && (
                          <div className="mt-3 p-3 rounded-md bg-muted/30 border border-border/30 text-xs text-foreground/80 whitespace-pre-wrap max-h-60 overflow-auto font-mono">
                            {parseSnapshot(version.content_snapshot)}
                          </div>
                        )}

                        {/* 对比视图 */}
                        {isDiffTarget && latestVersion && (
                          <div className="mt-3 p-3 rounded-md bg-muted/30 border border-border/30 max-h-60 overflow-auto">
                            <div className="text-[10px] text-muted-foreground mb-2">
                              {t('与当前版本 (v{n}) 对比').replace('{n}', String(latestVersion.version_number))}
                            </div>
                            {renderDiff(version.content_snapshot, latestVersion.content_snapshot)}
                          </div>
                        )}
                      </div>
                    </div>
                  )
                })}
              </div>
            </div>
          )}
        </div>

        {/* 底部 */}
        <div className="flex items-center justify-end px-6 py-4 border-t border-border/40">
          <button
            type="button"
            onClick={() => onOpenChange(false)}
            className="rounded-md border border-border/60 bg-background px-4 py-2 text-sm text-foreground hover:bg-muted/40 transition-colors"
          >
            {t('关闭')}
          </button>
        </div>
      </div>

      {/* 回滚确认对话框 */}
      {rollbackConfirm !== null && (
        <div className="fixed inset-0 z-[60] flex items-center justify-center">
          <div className="absolute inset-0 bg-black/30" onClick={() => setRollbackConfirm(null)} />
          <div className="relative w-full max-w-sm rounded-lg border border-border bg-background shadow-lg mx-4 p-6">
            <h3 className="text-sm font-semibold text-foreground mb-2">
              {t('确认回滚')}
            </h3>
            <p className="text-xs text-muted-foreground mb-4">
              {t('将创建一个新版本，内容恢复到 v{n} 的状态。此操作不会删除现有版本。').replace('{n}', String(rollbackConfirm))}
            </p>
            <div className="flex items-center justify-end gap-2">
              <button
                type="button"
                onClick={() => setRollbackConfirm(null)}
                disabled={rollbackLoading}
                className="rounded-md border border-border/60 bg-background px-3 py-1.5 text-xs text-foreground hover:bg-muted/40 transition-colors disabled:opacity-50"
              >
                {t('取消')}
              </button>
              <button
                type="button"
                onClick={() => handleRollback(rollbackConfirm)}
                disabled={rollbackLoading}
                className="rounded-md bg-red-500 px-3 py-1.5 text-xs font-medium text-white hover:bg-red-600 transition-colors disabled:opacity-50"
              >
                {rollbackLoading ? t('回滚中...') : t('确认回滚')}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

/** 解析 content_snapshot JSON，提取可读文本 */
function parseSnapshot(snapshot: string): string {
  try {
    const parsed = JSON.parse(snapshot)
    if (typeof parsed === 'string') return parsed
    if (parsed.content) return String(parsed.content)
    if (parsed.text) return String(parsed.text)
    return JSON.stringify(parsed, null, 2)
  } catch {
    return snapshot
  }
}

/** 格式化时间 */
function formatTime(isoString: string): string {
  try {
    const date = new Date(isoString)
    const now = new Date()
    const diffMs = now.getTime() - date.getTime()
    const diffMin = Math.floor(diffMs / 60000)
    const diffHour = Math.floor(diffMs / 3600000)
    const diffDay = Math.floor(diffMs / 86400000)

    if (diffMin < 1) return '刚刚'
    if (diffMin < 60) return `${diffMin} 分钟前`
    if (diffHour < 24) return `${diffHour} 小时前`
    if (diffDay < 7) return `${diffDay} 天前`

    return date.toLocaleDateString('zh-CN', {
      month: 'short',
      day: 'numeric',
      hour: '2-digit',
      minute: '2-digit',
    })
  } catch {
    return isoString
  }
}
