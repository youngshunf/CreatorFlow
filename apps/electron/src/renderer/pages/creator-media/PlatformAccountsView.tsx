import { useState, useEffect, useCallback } from 'react'
import { useT } from '@/context/LocaleContext'
import { useActiveWorkspace } from '@/context/AppShellContext'
import { useCreatorMedia } from './hooks/useCreatorMedia'
import { ProjectSwitcher } from './components/ProjectSwitcher'
import { AddPlatformAccountDialog } from './components/AddPlatformAccountDialog'
import type { PlatformAccount, Platform } from '@sprouty-ai/shared/db/types'
import { PLATFORM_MAP, PLATFORM_LIST } from '@sprouty-ai/shared/db/types'

/** 平台配置 — 从共享定义查找 */
const getPlatformConfig = (platform: string) => {
  const meta = PLATFORM_MAP[platform as Platform]
  return meta
    ? { label: meta.shortLabel, color: meta.color.split(' ')[0] }
    : { label: platform, color: 'text-foreground' }
}

/** 登录状态配置 */
const AUTH_STATUS_CONFIG: Record<string, { label: string; dotClass: string }> = {
  logged_in: { label: '已登录', dotClass: 'bg-green-500' },
  expired: { label: '登录过期', dotClass: 'bg-yellow-500' },
  not_logged_in: { label: '未登录', dotClass: 'bg-red-500' },
  error: { label: '登录异常', dotClass: 'bg-red-500' },
}

/** 认证方式标签 */
const AUTH_METHOD_LABELS: Record<string, string> = {
  cookie: 'Cookie',
  oauth: 'OAuth',
  api_key: 'API Key',
  browser_profile: '浏览器 Profile',
}

/**
 * 平台账号管理视图
 */
export default function PlatformAccountsView() {
  const t = useT()
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''
  const {
    projects, activeProject, loading, switchProject,
    checkBrowserAuth, deleteBrowserProfile,
    createPlatformAccount, updatePlatformAccount,
    launchBrowserLogin, generateBrowserFingerprint,
  } = useCreatorMedia()

  const [accounts, setAccounts] = useState<PlatformAccount[]>([])
  const [accountsLoading, setAccountsLoading] = useState(false)
  const [filterPlatform, setFilterPlatform] = useState<string>('all')
  const [filterStatus, setFilterStatus] = useState<string>('all')
  const [checkingId, setCheckingId] = useState<string | null>(null)
  const [deletingId, setDeletingId] = useState<string | null>(null)
  const [showAddDialog, setShowAddDialog] = useState(false)
  const [reloginAccount, setReloginAccount] = useState<PlatformAccount | null>(null)

  /** 加载平台账号列表 */
  const loadAccounts = useCallback(async () => {
    if (!workspaceId || !activeProject) {
      setAccounts([])
      return
    }
    setAccountsLoading(true)
    try {
      const list = await window.electronAPI.creatorMedia.platformAccounts.list(workspaceId, activeProject.id)
      setAccounts(list)
    } catch {
      setAccounts([])
    } finally {
      setAccountsLoading(false)
    }
  }, [workspaceId, activeProject])

  useEffect(() => {
    loadAccounts()
  }, [loadAccounts])

  /** 检查登录态 */
  const handleCheckAuth = useCallback(async (account: PlatformAccount) => {
    setCheckingId(account.id)
    try {
      const result = await checkBrowserAuth(account.id)
      // 刷新列表以反映最新状态
      await loadAccounts()
      if (!result.loggedIn) {
        // 可以在这里添加 toast 提示
      }
    } finally {
      setCheckingId(null)
    }
  }, [checkBrowserAuth, loadAccounts])

  /** 删除账号 */
  const handleDelete = useCallback(async (account: PlatformAccount) => {
    if (!workspaceId) return
    const confirmed = window.confirm(t(`确定要删除账号 "${account.nickname || account.platform_uid || account.platform}" 吗？此操作不可撤销。`))
    if (!confirmed) return

    setDeletingId(account.id)
    try {
      // 先删除浏览器 Profile
      await deleteBrowserProfile(account.id)
      // 再删除数据库记录
      await window.electronAPI.creatorMedia.platformAccounts.delete(workspaceId, account.id)
      await loadAccounts()
    } finally {
      setDeletingId(null)
    }
  }, [workspaceId, deleteBrowserProfile, loadAccounts, t])

  /** 过滤后的账号列表 */
  const filteredAccounts = accounts.filter((a) => {
    if (filterPlatform !== 'all' && a.platform !== filterPlatform) return false
    if (filterStatus !== 'all' && a.auth_status !== filterStatus) return false
    return true
  })

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
          <h1 className="text-base font-semibold text-foreground">{t('平台账号管理')}</h1>
          <p className="mt-0.5 text-xs text-muted-foreground">
            {activeProject ? activeProject.name : t('集中管理所有平台账号')}
          </p>
        </div>
        <div className="titlebar-no-drag flex items-center gap-2">
          <button
            type="button"
            onClick={() => {
              if (!activeProject) {
                window.alert(t('请先创建或选择一个项目'))
                return
              }
              setShowAddDialog(true)
            }}
            disabled={!activeProject}
            className="inline-flex items-center gap-1.5 rounded-md bg-foreground px-3 py-1.5 text-xs font-medium text-background hover:bg-foreground/90 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
            </svg>
            {t('添加账号')}
          </button>
          <ProjectSwitcher projects={projects} activeProject={activeProject} onSwitch={switchProject} />
        </div>
      </div>

      {/* 内容区 */}
      <div className="flex-1 overflow-auto px-6 py-6 space-y-4">
        {/* 筛选栏 */}
        {accounts.length > 0 && (
          <div className="flex items-center gap-3">
            <select
              value={filterPlatform}
              onChange={(e) => setFilterPlatform(e.target.value)}
              className="rounded-md border border-border/60 bg-background px-2.5 py-1.5 text-xs text-foreground"
            >
              <option value="all">{t('全部平台')}</option>
              {PLATFORM_LIST.map((p) => (
                <option key={p.id} value={p.id}>{t(p.shortLabel)}</option>
              ))}
            </select>
            <select
              value={filterStatus}
              onChange={(e) => setFilterStatus(e.target.value)}
              className="rounded-md border border-border/60 bg-background px-2.5 py-1.5 text-xs text-foreground"
            >
              <option value="all">{t('全部状态')}</option>
              {Object.entries(AUTH_STATUS_CONFIG).map(([id, cfg]) => (
                <option key={id} value={id}>{t(cfg.label)}</option>
              ))}
            </select>
          </div>
        )}

        {/* 账号卡片网格 */}
        {accountsLoading ? (
          <div className="flex items-center justify-center py-12">
            <p className="text-sm text-muted-foreground">{t('加载账号...')}</p>
          </div>
        ) : filteredAccounts.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
            {filteredAccounts.map((account) => (
              <AccountCard
                key={account.id}
                account={account}
                checking={checkingId === account.id}
                deleting={deletingId === account.id}
                onCheckAuth={() => handleCheckAuth(account)}
                onReLogin={() => setReloginAccount(account)}
                onDelete={() => handleDelete(account)}
              />
            ))}
          </div>
        ) : accounts.length > 0 ? (
          /* 有账号但筛选后为空 */
          <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-8 text-center">
            <p className="text-sm text-muted-foreground">{t('没有符合筛选条件的账号')}</p>
          </div>
        ) : (
          /* 空状态 */
          <EmptyState />
        )}
      </div>

      {/* 添加账号对话框 */}
      {activeProject && (
        <AddPlatformAccountDialog
          open={showAddDialog}
          onClose={() => setShowAddDialog(false)}
          onSuccess={loadAccounts}
          createPlatformAccount={createPlatformAccount}
          launchBrowserLogin={launchBrowserLogin}
          generateBrowserFingerprint={generateBrowserFingerprint}
          updatePlatformAccount={updatePlatformAccount}
          projectId={activeProject.id}
        />
      )}

      {/* 重新登录对话框 */}
      {activeProject && reloginAccount && (
        <AddPlatformAccountDialog
          open={!!reloginAccount}
          onClose={() => setReloginAccount(null)}
          onSuccess={loadAccounts}
          preselectedPlatform={reloginAccount.platform}
          existingAccountId={reloginAccount.id}
          createPlatformAccount={createPlatformAccount}
          launchBrowserLogin={launchBrowserLogin}
          generateBrowserFingerprint={generateBrowserFingerprint}
          updatePlatformAccount={updatePlatformAccount}
          projectId={activeProject.id}
        />
      )}
    </div>
  )
}

/** 账号卡片 */
function AccountCard({
  account,
  checking,
  deleting,
  onCheckAuth,
  onReLogin,
  onDelete,
}: {
  account: PlatformAccount
  checking: boolean
  deleting: boolean
  onCheckAuth: () => void
  onReLogin: () => void
  onDelete: () => void
}) {
  const t = useT()
  const platform = getPlatformConfig(account.platform)
  const authStatus = AUTH_STATUS_CONFIG[account.auth_status] || AUTH_STATUS_CONFIG.not_logged_in
  const authMethod = account.auth_method ? AUTH_METHOD_LABELS[account.auth_method] || account.auth_method : null

  return (
    <div className="rounded-lg border border-border/60 bg-background/40 p-4 space-y-3 hover:border-border transition-colors">
      {/* 平台 + 昵称 */}
      <div className="flex items-start justify-between">
        <div className="flex items-center gap-2">
          <span className={`text-sm font-semibold ${platform.color}`}>{t(platform.label)}</span>
          {account.is_primary === 1 && (
            <span className="inline-flex rounded-full bg-blue-100 dark:bg-blue-900/30 px-1.5 py-0.5 text-[10px] text-blue-700 dark:text-blue-400">
              {t('主账号')}
            </span>
          )}
        </div>
        {/* 登录状态指示灯 */}
        <div className="flex items-center gap-1.5">
          <span className={`w-2 h-2 rounded-full ${authStatus.dotClass}`} />
          <span className="text-xs text-muted-foreground">{t(authStatus.label)}</span>
        </div>
      </div>

      {/* 账号信息 */}
      <div className="space-y-1">
        {account.nickname && (
          <p className="text-sm font-medium text-foreground">{account.nickname}</p>
        )}
        {account.platform_uid && (
          <p className="text-xs text-muted-foreground">@{account.platform_uid}</p>
        )}
        {account.followers > 0 && (
          <p className="text-xs text-muted-foreground">
            {formatNumber(account.followers)} {t('粉丝')}
            {account.total_likes > 0 && ` · ${formatNumber(account.total_likes)} ${t('获赞')}`}
          </p>
        )}
      </div>

      {/* 认证方式标签 */}
      {authMethod && (
        <div>
          <span className="inline-flex rounded-full bg-muted/60 px-2 py-0.5 text-[10px] text-muted-foreground">
            {authMethod}
          </span>
        </div>
      )}

      {/* 操作按钮 */}
      <div className="flex items-center gap-2 pt-1 border-t border-border/30">
        <button
          type="button"
          onClick={onCheckAuth}
          disabled={checking}
          className="rounded-md px-2.5 py-1 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors disabled:opacity-50"
        >
          {checking ? t('检查中...') : t('检查登录态')}
        </button>
        <button
          type="button"
          onClick={onReLogin}
          className="rounded-md px-2.5 py-1 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/50 transition-colors"
        >
          {t('重新登录')}
        </button>
        <button
          type="button"
          onClick={onDelete}
          disabled={deleting}
          className="rounded-md px-2.5 py-1 text-xs text-red-500 hover:text-red-600 hover:bg-red-50 dark:hover:bg-red-900/20 transition-colors disabled:opacity-50 ml-auto"
        >
          {deleting ? t('删除中...') : t('删除')}
        </button>
      </div>
    </div>
  )
}

/** 空状态 */
function EmptyState() {
  const t = useT()
  return (
    <div className="rounded-lg border border-dashed border-border/60 bg-background/40 px-4 py-12">
      <div className="text-center space-y-3">
        <svg className="mx-auto w-10 h-10 text-muted-foreground/40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
          <path strokeLinecap="round" strokeLinejoin="round" d="M15 19.128a9.38 9.38 0 0 0 2.625.372 9.337 9.337 0 0 0 4.121-.952 4.125 4.125 0 0 0-7.533-2.493M15 19.128v-.003c0-1.113-.285-2.16-.786-3.07M15 19.128v.106A12.318 12.318 0 0 1 8.624 21c-2.331 0-4.512-.645-6.374-1.766l-.001-.109a6.375 6.375 0 0 1 11.964-3.07M12 6.375a3.375 3.375 0 1 1-6.75 0 3.375 3.375 0 0 1 6.75 0Zm8.25 2.25a2.625 2.625 0 1 1-5.25 0 2.625 2.625 0 0 1 5.25 0Z" />
        </svg>
        <div>
          <p className="text-sm font-medium text-foreground">{t('还没有添加平台账号')}</p>
          <p className="mt-1 text-xs text-muted-foreground">
            {t('添加你的社交媒体账号，管理登录态并自动发布内容')}
          </p>
        </div>
        <div className="flex items-center justify-center gap-2 pt-1">
          {PLATFORM_LIST.map(p => p.shortLabel).map((name) => (
            <span key={name} className="inline-flex items-center rounded-full bg-muted/60 px-2.5 py-0.5 text-[10px] text-muted-foreground">
              {t(name)}
            </span>
          ))}
        </div>
      </div>
    </div>
  )
}

/** 格式化数字（如 12500 → 1.2w） */
function formatNumber(n: number): string {
  if (n >= 10000) return `${(n / 10000).toFixed(1)}w`
  if (n >= 1000) return `${(n / 1000).toFixed(1)}k`
  return String(n)
}
