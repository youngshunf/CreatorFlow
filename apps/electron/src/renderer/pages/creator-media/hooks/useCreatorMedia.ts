import { useState, useEffect, useCallback, useMemo } from 'react'
import { useActiveWorkspace } from '@/context/AppShellContext'
import type { ContentStats, PublishStats } from '../types'
import type {
  Project, AccountProfile, Content, PublishRecord, ViralPattern, Competitor,
  CreateProject, UpdateProject, CreateAccountProfile, CreateContent,
  CreateViralPattern, UpdateViralPattern,
} from '@sprouty-ai/shared/db/types'

/**
 * 创作面板数据 Hook — 封装 window.electronAPI.creatorMedia.* IPC 调用
 */
export function useCreatorMedia() {
  const activeWorkspace = useActiveWorkspace()
  const workspaceId = activeWorkspace?.id || ''

  const [projects, setProjects] = useState<Project[]>([])
  const [activeProject, setActiveProject] = useState<Project | null>(null)
  const [contents, setContents] = useState<Content[]>([])
  const [profile, setProfile] = useState<AccountProfile | null>(null)
  const [publishRecords, setPublishRecords] = useState<PublishRecord[]>([])
  const [viralPatterns, setViralPatterns] = useState<ViralPattern[]>([])
  const [competitors, setCompetitors] = useState<Competitor[]>([])
  const [loading, setLoading] = useState(true)

  /** 加载项目列表和活跃项目 */
  const loadProjects = useCallback(async () => {
    if (!workspaceId) return
    try {
      const [list, active] = await Promise.all([
        window.electronAPI.creatorMedia.projects.list(workspaceId),
        window.electronAPI.creatorMedia.projects.getActive(workspaceId),
      ])
      setProjects(list)
      setActiveProject(active)
    } catch {
      setProjects([])
      setActiveProject(null)
    }
  }, [workspaceId])

  /** 加载内容列表 */
  const loadContents = useCallback(async (projectId: string) => {
    if (!workspaceId) return
    try {
      const result = await window.electronAPI.creatorMedia.contents.list(workspaceId, projectId)
      setContents(result)
    } catch {
      setContents([])
    }
  }, [workspaceId])

  /** 加载账号画像 */
  const loadProfile = useCallback(async (projectId: string) => {
    if (!workspaceId) return
    try {
      const result = await window.electronAPI.creatorMedia.profiles.get(workspaceId, projectId)
      setProfile(result)
    } catch {
      setProfile(null)
    }
  }, [workspaceId])

  /** 加载发布记录（项目下所有内容的发布记录） */
  const loadPublishRecords = useCallback(async (projectId: string) => {
    if (!workspaceId) return
    try {
      // 先获取项目下所有内容，再聚合发布记录
      const contentList = await window.electronAPI.creatorMedia.contents.list(workspaceId, projectId)
      const allRecords: PublishRecord[] = []
      for (const c of contentList) {
        try {
          const records = await window.electronAPI.creatorMedia.publishRecords.list(workspaceId, c.id)
          allRecords.push(...records)
        } catch {
          // 单条失败不影响整体
        }
      }
      setPublishRecords(allRecords)
    } catch {
      setPublishRecords([])
    }
  }, [workspaceId])

  /** 加载爆款模式 */
  const loadViralPatterns = useCallback(async (filters?: { projectId?: string }) => {
    if (!workspaceId) return
    try {
      const result = await window.electronAPI.creatorMedia.viralPatterns.list(workspaceId, filters)
      setViralPatterns(result)
    } catch {
      setViralPatterns([])
    }
  }, [workspaceId])

  /** 加载竞品 */
  const loadCompetitors = useCallback(async (projectId: string) => {
    if (!workspaceId) return
    try {
      const result = await window.electronAPI.creatorMedia.competitors.list(workspaceId, projectId)
      setCompetitors(result)
    } catch {
      setCompetitors([])
    }
  }, [workspaceId])

  /** 创建爆款模式 */
  const createViralPattern = useCallback(async (data: Omit<CreateViralPattern, 'id'>) => {
    if (!workspaceId) return null
    const id = crypto.randomUUID()
    const result = await window.electronAPI.creatorMedia.viralPatterns.create(workspaceId, { ...data, id })
    await loadViralPatterns({ projectId: activeProject?.id })
    return result
  }, [workspaceId, activeProject?.id, loadViralPatterns])

  /** 更新爆款模式 */
  const updateViralPattern = useCallback(async (id: string, data: UpdateViralPattern) => {
    if (!workspaceId) return null
    const result = await window.electronAPI.creatorMedia.viralPatterns.update(workspaceId, id, data)
    await loadViralPatterns({ projectId: activeProject?.id })
    return result
  }, [workspaceId, activeProject?.id, loadViralPatterns])

  /** 删除爆款模式 */
  const deleteViralPattern = useCallback(async (id: string) => {
    if (!workspaceId) return false
    const ok = await window.electronAPI.creatorMedia.viralPatterns.delete(workspaceId, id)
    if (ok) await loadViralPatterns({ projectId: activeProject?.id })
    return ok
  }, [workspaceId, activeProject?.id, loadViralPatterns])

  /** 切换项目 */
  const switchProject = useCallback(async (projectId: string) => {
    if (!workspaceId) return
    try {
      await window.electronAPI.creatorMedia.projects.setActive(workspaceId, projectId)
      await loadProjects()
    } catch {
      // 静默处理
    }
  }, [workspaceId, loadProjects])

  /** 创建项目，创建后自动设为活跃并刷新 */
  const createProject = useCallback(async (data: Omit<CreateProject, 'id' | 'is_active'>) => {
    if (!workspaceId) return null
    const id = crypto.randomUUID()
    const project = await window.electronAPI.creatorMedia.projects.create(workspaceId, {
      ...data,
      id,
      is_active: 1,
    })
    await window.electronAPI.creatorMedia.projects.setActive(workspaceId, project.id)
    await loadProjects()
    return project
  }, [workspaceId, loadProjects])

  /** 更新项目 */
  const updateProject = useCallback(async (projectId: string, data: UpdateProject) => {
    if (!workspaceId) return null
    const result = await window.electronAPI.creatorMedia.projects.update(workspaceId, projectId, data)
    await loadProjects()
    return result
  }, [workspaceId, loadProjects])

  /** 删除项目，删除后切换到下一个项目 */
  const deleteProject = useCallback(async (projectId: string) => {
    if (!workspaceId) return false
    const ok = await window.electronAPI.creatorMedia.projects.delete(workspaceId, projectId)
    if (ok) {
      // 刷新列表，自动切换到第一个可用项目
      const list = await window.electronAPI.creatorMedia.projects.list(workspaceId)
      setProjects(list)
      if (list.length > 0) {
        await window.electronAPI.creatorMedia.projects.setActive(workspaceId, list[0].id)
        setActiveProject(list[0])
      } else {
        setActiveProject(null)
        setContents([])
        setProfile(null)
      }
    }
    return ok
  }, [workspaceId])

  /** 更新/创建账号画像 */
  const upsertProfile = useCallback(async (data: Omit<CreateAccountProfile, 'id'>) => {
    if (!workspaceId) return null
    const result = await window.electronAPI.creatorMedia.profiles.upsert(workspaceId, {
      ...data,
      id: profile?.id || crypto.randomUUID(),
    })
    setProfile(result)
    return result
  }, [workspaceId, profile?.id])

  /** 创建内容 */
  const createContent = useCallback(async (data: Omit<CreateContent, 'id' | 'project_id'>) => {
    if (!workspaceId || !activeProject) return null
    const id = crypto.randomUUID()
    const result = await window.electronAPI.creatorMedia.contents.create(workspaceId, {
      ...data,
      id,
      project_id: activeProject.id,
    })
    await loadContents(activeProject.id)
    return result
  }, [workspaceId, activeProject, loadContents])

  /** 删除内容 */
  const deleteContent = useCallback(async (contentId: string) => {
    if (!workspaceId || !activeProject) return false
    const ok = await window.electronAPI.creatorMedia.contents.delete(workspaceId, contentId)
    if (ok) await loadContents(activeProject.id)
    return ok
  }, [workspaceId, activeProject, loadContents])

  /** 更新内容状态 */
  const updateContentStatus = useCallback(async (contentId: string, status: string) => {
    if (!workspaceId || !activeProject) return null
    const result = await window.electronAPI.creatorMedia.contents.updateStatus(workspaceId, contentId, status)
    await loadContents(activeProject.id)
    return result
  }, [workspaceId, activeProject, loadContents])

  /** 计算统计数据 */
  const stats: ContentStats = useMemo(() => {
    const base: ContentStats = {
      total: contents.length,
      idea: 0,
      researching: 0,
      scripting: 0,
      creating: 0,
      reviewing: 0,
      scheduled: 0,
      published: 0,
      archived: 0,
    }
    for (const c of contents) {
      const key = c.status as keyof Omit<ContentStats, 'total'>
      if (key in base) {
        base[key]++
      }
    }
    return base
  }, [contents])

  /** 计算发布统计数据 */
  const publishStats: PublishStats = useMemo(() => {
    const base: PublishStats = {
      totalViews: 0, totalLikes: 0, totalComments: 0,
      totalFavorites: 0, totalShares: 0, successCount: 0, failedCount: 0,
    }
    for (const r of publishRecords) {
      base.totalViews += r.views || 0
      base.totalLikes += r.likes || 0
      base.totalComments += r.comments || 0
      base.totalFavorites += r.favorites || 0
      base.totalShares += r.shares || 0
      if (r.status === 'success') base.successCount++
      if (r.status === 'failed') base.failedCount++
    }
    return base
  }, [publishRecords])

  /** 初始加载 */
  useEffect(() => {
    if (!workspaceId) {
      setLoading(false)
      return
    }
    setLoading(true)
    loadProjects().finally(() => setLoading(false))
  }, [workspaceId, loadProjects])

  /** activeProject 变化时重新加载数据 */
  useEffect(() => {
    if (activeProject) {
      loadContents(activeProject.id)
      loadProfile(activeProject.id)
      loadPublishRecords(activeProject.id)
      loadViralPatterns({ projectId: activeProject.id })
      loadCompetitors(activeProject.id)
    } else {
      setContents([])
      setProfile(null)
      setPublishRecords([])
      setViralPatterns([])
      setCompetitors([])
    }
  }, [activeProject, loadContents, loadProfile, loadPublishRecords, loadViralPatterns, loadCompetitors])

  return {
    projects,
    activeProject,
    contents,
    profile,
    stats,
    publishRecords,
    viralPatterns,
    competitors,
    publishStats,
    loading,
    switchProject,
    refresh: loadProjects,
    createProject,
    updateProject,
    deleteProject,
    upsertProfile,
    createContent,
    deleteContent,
    updateContentStatus,
    createViralPattern,
    updateViralPattern,
    deleteViralPattern,
  }
}
