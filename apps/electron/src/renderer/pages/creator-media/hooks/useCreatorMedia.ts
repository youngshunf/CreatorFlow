import { useState, useEffect, useCallback, useMemo } from "react";
import { useActiveWorkspace } from "@/context/AppShellContext";
import type { ContentStats, PublishStats } from "../types";
import type {
  Project,
  AccountProfile,
  Content,
  PublishRecord,
  ViralPattern,
  Competitor,
  PlatformAccount,
  CreateProject,
  UpdateProject,
  CreateAccountProfile,
  CreateContent,
  CreateViralPattern,
  UpdateViralPattern,
  CreatePlatformAccount,
  UpdatePlatformAccount,
  ReviewTask,
  CreateReviewTask,
  UpdateReviewTask,
  ContentVersion,
  CreateContentVersionInput,
  PublishQueueItem,
  CreatePublishQueueInput,
  Platform,
} from "@sprouty-ai/shared/db/types";
import type { BrowserFingerprint } from "@sprouty-ai/shared/services/browser-profile-manager";

/**
 * 创作面板数据 Hook — 封装 window.electronAPI.creatorMedia.* IPC 调用
 */
export function useCreatorMedia() {
  const activeWorkspace = useActiveWorkspace();
  const workspaceId = activeWorkspace?.id || "";

  const [projects, setProjects] = useState<Project[]>([]);
  const [activeProject, setActiveProject] = useState<Project | null>(null);
  const [contents, setContents] = useState<Content[]>([]);
  const [profile, setProfile] = useState<AccountProfile | null>(null);
  const [publishRecords, setPublishRecords] = useState<PublishRecord[]>([]);
  const [viralPatterns, setViralPatterns] = useState<ViralPattern[]>([]);
  const [competitors, setCompetitors] = useState<Competitor[]>([]);
  const [platformAccounts, setPlatformAccounts] = useState<PlatformAccount[]>(
    [],
  );
  const [loading, setLoading] = useState(true);

  /** 加载项目列表和活跃项目 */
  const loadProjects = useCallback(async () => {
    if (!workspaceId) return;
    try {
      const [list, active] = await Promise.all([
        window.electronAPI.creatorMedia.projects.list(workspaceId),
        window.electronAPI.creatorMedia.projects.getActive(workspaceId),
      ]);
      setProjects(list);
      setActiveProject(active);
    } catch {
      setProjects([]);
      setActiveProject(null);
    }
  }, [workspaceId]);

  /** 加载内容列表 */
  const loadContents = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        const result = await window.electronAPI.creatorMedia.contents.list(
          workspaceId,
          projectId,
        );
        setContents(result);
      } catch {
        setContents([]);
      }
    },
    [workspaceId],
  );

  /** 加载账号画像 */
  const loadProfile = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        const result = await window.electronAPI.creatorMedia.profiles.get(
          workspaceId,
          projectId,
        );
        setProfile(result);
      } catch {
        setProfile(null);
      }
    },
    [workspaceId],
  );

  /** 加载发布记录（项目下所有内容的发布记录） */
  const loadPublishRecords = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        // 先获取项目下所有内容，再聚合发布记录
        const contentList = await window.electronAPI.creatorMedia.contents.list(
          workspaceId,
          projectId,
        );
        const allRecords: PublishRecord[] = [];
        for (const c of contentList) {
          try {
            const records =
              await window.electronAPI.creatorMedia.publishRecords.list(
                workspaceId,
                c.id,
              );
            allRecords.push(...records);
          } catch {
            // 单条失败不影响整体
          }
        }
        setPublishRecords(allRecords);
      } catch {
        setPublishRecords([]);
      }
    },
    [workspaceId],
  );

  /** 加载爆款模式 */
  const loadViralPatterns = useCallback(
    async (filters?: { projectId?: string }) => {
      if (!workspaceId) return;
      try {
        const result = await window.electronAPI.creatorMedia.viralPatterns.list(
          workspaceId,
          filters,
        );
        setViralPatterns(result);
      } catch {
        setViralPatterns([]);
      }
    },
    [workspaceId],
  );

  /** 加载竞品 */
  const loadCompetitors = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        const result = await window.electronAPI.creatorMedia.competitors.list(
          workspaceId,
          projectId,
        );
        setCompetitors(result);
      } catch {
        setCompetitors([]);
      }
    },
    [workspaceId],
  );

  /** 加载平台账号 */
  const loadPlatformAccounts = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        const result =
          await window.electronAPI.creatorMedia.platformAccounts.list(
            workspaceId,
            projectId,
          );
        setPlatformAccounts(result);
      } catch {
        setPlatformAccounts([]);
      }
    },
    [workspaceId],
  );

  /** 创建平台账号 */
  const createPlatformAccount = useCallback(
    async (data: Omit<CreatePlatformAccount, "id">) => {
      if (!workspaceId || !activeProject) return null;
      const id = crypto.randomUUID();
      const result =
        await window.electronAPI.creatorMedia.platformAccounts.create(
          workspaceId,
          { ...data, id },
        );
      await loadPlatformAccounts(activeProject.id);
      return result;
    },
    [workspaceId, activeProject, loadPlatformAccounts],
  );

  /** 更新平台账号 */
  const updatePlatformAccount = useCallback(
    async (id: string, data: UpdatePlatformAccount) => {
      if (!workspaceId || !activeProject) return null;
      const result =
        await window.electronAPI.creatorMedia.platformAccounts.update(
          workspaceId,
          id,
          data,
        );
      await loadPlatformAccounts(activeProject.id);
      return result;
    },
    [workspaceId, activeProject, loadPlatformAccounts],
  );

  /** 删除平台账号 */
  const deletePlatformAccount = useCallback(
    async (id: string) => {
      if (!workspaceId || !activeProject) return false;
      const ok = await window.electronAPI.creatorMedia.platformAccounts.delete(
        workspaceId,
        id,
      );
      if (ok) await loadPlatformAccounts(activeProject.id);
      return ok;
    },
    [workspaceId, activeProject, loadPlatformAccounts],
  );

  /** 创建爆款模式 */
  const createViralPattern = useCallback(
    async (data: Omit<CreateViralPattern, "id">) => {
      if (!workspaceId) return null;
      const id = crypto.randomUUID();
      const result = await window.electronAPI.creatorMedia.viralPatterns.create(
        workspaceId,
        { ...data, id },
      );
      await loadViralPatterns({ projectId: activeProject?.id });
      return result;
    },
    [workspaceId, activeProject?.id, loadViralPatterns],
  );

  /** 更新爆款模式 */
  const updateViralPattern = useCallback(
    async (id: string, data: UpdateViralPattern) => {
      if (!workspaceId) return null;
      const result = await window.electronAPI.creatorMedia.viralPatterns.update(
        workspaceId,
        id,
        data,
      );
      await loadViralPatterns({ projectId: activeProject?.id });
      return result;
    },
    [workspaceId, activeProject?.id, loadViralPatterns],
  );

  /** 删除爆款模式 */
  const deleteViralPattern = useCallback(
    async (id: string) => {
      if (!workspaceId) return false;
      const ok = await window.electronAPI.creatorMedia.viralPatterns.delete(
        workspaceId,
        id,
      );
      if (ok) await loadViralPatterns({ projectId: activeProject?.id });
      return ok;
    },
    [workspaceId, activeProject?.id, loadViralPatterns],
  );

  /** 切换项目 */
  const switchProject = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return;
      try {
        await window.electronAPI.creatorMedia.projects.setActive(
          workspaceId,
          projectId,
        );
        await loadProjects();
      } catch {
        // 静默处理
      }
    },
    [workspaceId, loadProjects],
  );

  /** 创建项目，创建后自动设为活跃并刷新 */
  const createProject = useCallback(
    async (data: Omit<CreateProject, "id" | "is_active">) => {
      if (!workspaceId) return null;
      const id = crypto.randomUUID();
      const project = await window.electronAPI.creatorMedia.projects.create(
        workspaceId,
        {
          ...data,
          id,
          is_active: 1,
        },
      );
      await window.electronAPI.creatorMedia.projects.setActive(
        workspaceId,
        project.id,
      );
      await loadProjects();
      return project;
    },
    [workspaceId, loadProjects],
  );

  /** 更新项目 */
  const updateProject = useCallback(
    async (projectId: string, data: UpdateProject) => {
      if (!workspaceId) return null;
      const result = await window.electronAPI.creatorMedia.projects.update(
        workspaceId,
        projectId,
        data,
      );
      await loadProjects();
      return result;
    },
    [workspaceId, loadProjects],
  );

  /** 删除项目，删除后切换到下一个项目 */
  const deleteProject = useCallback(
    async (projectId: string) => {
      if (!workspaceId) return false;
      const ok = await window.electronAPI.creatorMedia.projects.delete(
        workspaceId,
        projectId,
      );
      if (ok) {
        // 刷新列表，自动切换到第一个可用项目
        const list =
          await window.electronAPI.creatorMedia.projects.list(workspaceId);
        setProjects(list);
        if (list.length > 0) {
          await window.electronAPI.creatorMedia.projects.setActive(
            workspaceId,
            list[0].id,
          );
          setActiveProject(list[0]);
        } else {
          setActiveProject(null);
          setContents([]);
          setProfile(null);
        }
      }
      return ok;
    },
    [workspaceId],
  );

  /** 更新/创建账号画像 */
  const upsertProfile = useCallback(
    async (data: Omit<CreateAccountProfile, "id">) => {
      if (!workspaceId) return null;
      const result = await window.electronAPI.creatorMedia.profiles.upsert(
        workspaceId,
        {
          ...data,
          id: profile?.id || crypto.randomUUID(),
        },
      );
      setProfile(result);
      return result;
    },
    [workspaceId, profile?.id],
  );

  /** 创建内容 */
  const createContent = useCallback(
    async (data: Omit<CreateContent, "id" | "project_id">) => {
      if (!workspaceId || !activeProject) return null;
      const id = crypto.randomUUID();
      const result = await window.electronAPI.creatorMedia.contents.create(
        workspaceId,
        {
          ...data,
          id,
          project_id: activeProject.id,
        },
      );
      await loadContents(activeProject.id);
      return result;
    },
    [workspaceId, activeProject, loadContents],
  );

  /** 删除内容 */
  const deleteContent = useCallback(
    async (contentId: string) => {
      if (!workspaceId || !activeProject) return false;
      const ok = await window.electronAPI.creatorMedia.contents.delete(
        workspaceId,
        contentId,
      );
      if (ok) await loadContents(activeProject.id);
      return ok;
    },
    [workspaceId, activeProject, loadContents],
  );

  /** 更新内容状态 */
  const updateContentStatus = useCallback(
    async (contentId: string, status: string) => {
      if (!workspaceId || !activeProject) return null;
      const result =
        await window.electronAPI.creatorMedia.contents.updateStatus(
          workspaceId,
          contentId,
          status,
        );
      await loadContents(activeProject.id);
      return result;
    },
    [workspaceId, activeProject, loadContents],
  );

  /** 更新内容轨道 */
  const updateContentTracks = useCallback(
    async (contentId: string, tracks: string) => {
      if (!workspaceId || !activeProject) return null;
      const result = await window.electronAPI.creatorMedia.contents.update(
        workspaceId,
        contentId,
        { content_tracks: tracks },
      );
      await loadContents(activeProject.id);
      return result;
    },
    [workspaceId, activeProject, loadContents],
  );

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
    };
    for (const c of contents) {
      const key = c.status as keyof Omit<ContentStats, "total">;
      if (key in base) {
        base[key]++;
      }
    }
    return base;
  }, [contents]);

  /** 计算发布统计数据 */
  const publishStats: PublishStats = useMemo(() => {
    const base: PublishStats = {
      totalViews: 0,
      totalLikes: 0,
      totalComments: 0,
      totalFavorites: 0,
      totalShares: 0,
      successCount: 0,
      failedCount: 0,
    };
    for (const r of publishRecords) {
      base.totalViews += r.views || 0;
      base.totalLikes += r.likes || 0;
      base.totalComments += r.comments || 0;
      base.totalFavorites += r.favorites || 0;
      base.totalShares += r.shares || 0;
      if (r.status === "success") base.successCount++;
      if (r.status === "failed") base.failedCount++;
    }
    return base;
  }, [publishRecords]);

  /** 初始加载 */
  useEffect(() => {
    if (!workspaceId) {
      setLoading(false);
      return;
    }
    setLoading(true);
    loadProjects().finally(() => setLoading(false));
  }, [workspaceId, loadProjects]);

  /** activeProject 变化时重新加载数据 */
  useEffect(() => {
    if (activeProject) {
      loadContents(activeProject.id);
      loadProfile(activeProject.id);
      loadPublishRecords(activeProject.id);
      loadViralPatterns({ projectId: activeProject.id });
      loadCompetitors(activeProject.id);
      loadPlatformAccounts(activeProject.id);
    } else {
      setContents([]);
      setProfile(null);
      setPublishRecords([]);
      setViralPatterns([]);
      setCompetitors([]);
      setPlatformAccounts([]);
    }
  }, [
    activeProject,
    loadContents,
    loadProfile,
    loadPublishRecords,
    loadViralPatterns,
    loadCompetitors,
    loadPlatformAccounts,
  ]);

  // ============================================================
  // 采集调度任务
  // ============================================================

  /** 获取发布记录的采集任务列表 */
  const listReviewTasks = useCallback(
    async (publishRecordId: string): Promise<ReviewTask[]> => {
      if (!workspaceId) return [];
      try {
        return await window.electronAPI.creatorMedia.reviewTasks.list(
          workspaceId,
          publishRecordId,
        );
      } catch {
        return [];
      }
    },
    [workspaceId],
  );

  /** 创建采集任务 */
  const createReviewTask = useCallback(
    async (data: CreateReviewTask): Promise<ReviewTask | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.reviewTasks.create(
        workspaceId,
        data,
      );
    },
    [workspaceId],
  );

  /** 更新采集任务 */
  const updateReviewTask = useCallback(
    async (id: string, data: UpdateReviewTask): Promise<ReviewTask | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.reviewTasks.update(
        workspaceId,
        id,
        data,
      );
    },
    [workspaceId],
  );

  /** 取消某发布记录的所有待执行采集任务 */
  const cancelReviewTasks = useCallback(
    async (publishRecordId: string): Promise<number> => {
      if (!workspaceId) return 0;
      return await window.electronAPI.creatorMedia.reviewTasks.cancel(
        workspaceId,
        publishRecordId,
      );
    },
    [workspaceId],
  );

  // ============================================================
  // 内容版本管理
  // ============================================================

  /** 获取内容的版本列表 */
  const listContentVersions = useCallback(
    async (contentId: string): Promise<ContentVersion[]> => {
      if (!workspaceId) return [];
      try {
        return await window.electronAPI.creatorMedia.contentVersions.list(
          workspaceId,
          contentId,
        );
      } catch {
        return [];
      }
    },
    [workspaceId],
  );

  /** 获取单个版本 */
  const getContentVersion = useCallback(
    async (id: string): Promise<ContentVersion | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.contentVersions.get(
        workspaceId,
        id,
      );
    },
    [workspaceId],
  );

  /** 创建版本快照 */
  const createContentVersion = useCallback(
    async (data: CreateContentVersionInput): Promise<ContentVersion | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.contentVersions.create(
        workspaceId,
        data,
      );
    },
    [workspaceId],
  );

  /** 回滚到指定版本 */
  const rollbackContentVersion = useCallback(
    async (
      contentId: string,
      versionNumber: number,
    ): Promise<ContentVersion | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.contentVersions.rollback(
        workspaceId,
        contentId,
        versionNumber,
      );
    },
    [workspaceId],
  );

  // ============================================================
  // 发布队列
  // ============================================================

  /** 获取内容的发布队列列表 */
  const listPublishQueue = useCallback(
    async (contentId: string): Promise<PublishQueueItem[]> => {
      if (!workspaceId) return [];
      try {
        return await window.electronAPI.creatorMedia.publishQueue.list(
          workspaceId,
          contentId,
        );
      } catch {
        return [];
      }
    },
    [workspaceId],
  );

  /** 入队 */
  const enqueuePublish = useCallback(
    async (data: CreatePublishQueueInput): Promise<PublishQueueItem | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.publishQueue.enqueue(
        workspaceId,
        data,
      );
    },
    [workspaceId],
  );

  /** 取消某内容的所有待处理发布任务 */
  const cancelPublishQueue = useCallback(
    async (contentId: string): Promise<number> => {
      if (!workspaceId) return 0;
      return await window.electronAPI.creatorMedia.publishQueue.cancel(
        workspaceId,
        contentId,
      );
    },
    [workspaceId],
  );

  /** 获取某平台账号的下一个待处理发布任务 */
  const getNextPublishQueueItem = useCallback(
    async (platformAccountId: string): Promise<PublishQueueItem | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.publishQueue.getNext(
        workspaceId,
        platformAccountId,
      );
    },
    [workspaceId],
  );

  // ============================================================
  // 采集调度服务状态
  // ============================================================

  /** 获取采集调度服务状态 */
  const getReviewSchedulerStatus = useCallback(async (): Promise<{
    running: boolean;
  }> => {
    return await window.electronAPI.creatorMedia.reviewScheduler.status();
  }, []);

  // ============================================================
  // 浏览器 Profile 管理
  // ============================================================

  /** 启动登录浏览器（有头模式） */
  const launchBrowserLogin = useCallback(
    async (
      platformAccountId: string,
      platform: Platform,
    ): Promise<{ success: boolean; error?: string }> => {
      if (!workspaceId) return { success: false, error: "No workspace" };
      return await window.electronAPI.creatorMedia.browser.launchLogin(
        workspaceId,
        platformAccountId,
        platform,
      );
    },
    [workspaceId],
  );

  /** 检查登录态 */
  const checkBrowserAuth = useCallback(
    async (
      platformAccountId: string,
    ): Promise<{ loggedIn: boolean; error?: string }> => {
      if (!workspaceId) return { loggedIn: false, error: "No workspace" };
      return await window.electronAPI.creatorMedia.browser.checkAuth(
        workspaceId,
        platformAccountId,
      );
    },
    [workspaceId],
  );

  /** 检查 Profile 是否存在 */
  const browserProfileExists = useCallback(
    async (platformAccountId: string): Promise<boolean> => {
      if (!workspaceId) return false;
      return await window.electronAPI.creatorMedia.browser.profileExists(
        workspaceId,
        platformAccountId,
      );
    },
    [workspaceId],
  );

  /** 删除 Profile */
  const deleteBrowserProfile = useCallback(
    async (platformAccountId: string): Promise<boolean> => {
      if (!workspaceId) return false;
      return await window.electronAPI.creatorMedia.browser.deleteProfile(
        workspaceId,
        platformAccountId,
      );
    },
    [workspaceId],
  );

  /** 生成并保存指纹 */
  const generateBrowserFingerprint = useCallback(
    async (platformAccountId: string): Promise<BrowserFingerprint | null> => {
      if (!workspaceId) return null;
      return await window.electronAPI.creatorMedia.browser.generateFingerprint(
        workspaceId,
        platformAccountId,
      );
    },
    [workspaceId],
  );

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
    refreshContents: () => activeProject ? loadContents(activeProject.id) : Promise.resolve(),
    createProject,
    updateProject,
    deleteProject,
    upsertProfile,
    createContent,
    deleteContent,
    updateContentStatus,
    updateContentTracks,
    createViralPattern,
    updateViralPattern,
    deleteViralPattern,
    // 平台账号
    platformAccounts,
    createPlatformAccount,
    updatePlatformAccount,
    deletePlatformAccount,
    refreshPlatformAccounts: loadPlatformAccounts,
    // 采集调度任务
    listReviewTasks,
    createReviewTask,
    updateReviewTask,
    cancelReviewTasks,
    // 内容版本管理
    listContentVersions,
    getContentVersion,
    createContentVersion,
    rollbackContentVersion,
    // 发布队列
    listPublishQueue,
    enqueuePublish,
    cancelPublishQueue,
    getNextPublishQueueItem,
    // 采集调度服务
    getReviewSchedulerStatus,
    // 浏览器 Profile 管理
    launchBrowserLogin,
    checkBrowserAuth,
    browserProfileExists,
    deleteBrowserProfile,
    generateBrowserFingerprint,
  };
}
