/**
 * Electron 主进程 — 自媒体创作 IPC 通道注册
 *
 * 所有 creatorMedia:* IPC handler
 */

import { ipcMain } from 'electron';
import * as fs from 'node:fs';
import * as path from 'node:path';
import * as os from 'node:os';
import type { WindowManager } from './window-manager';
import { getCreatorMediaDB } from './creator-media-db';
import { IPC_CHANNELS } from '../shared/types';
import { getWorkspaceByNameOrId, type Workspace } from '@sprouty-ai/shared/config';
import { ipcLog } from './logger';

// Repository 导入
import * as projectsRepo from '@sprouty-ai/shared/db/repositories/projects';
import * as profilesRepo from '@sprouty-ai/shared/db/repositories/profiles';
import * as platformAccountsRepo from '@sprouty-ai/shared/db/repositories/platform-accounts';
import * as competitorsRepo from '@sprouty-ai/shared/db/repositories/competitors';
import * as contentsRepo from '@sprouty-ai/shared/db/repositories/contents';
import * as publishRecordsRepo from '@sprouty-ai/shared/db/repositories/publish-records';
import * as viralPatternsRepo from '@sprouty-ai/shared/db/repositories/viral-patterns';
import * as reviewTasksRepo from '@sprouty-ai/shared/db/repositories/review-tasks';
import * as contentVersionsRepo from '@sprouty-ai/shared/db/repositories/content-versions';
import * as publishQueueRepo from '@sprouty-ai/shared/db/repositories/publish-queue';
import * as draftsRepo from '@sprouty-ai/shared/db/repositories/drafts';
import * as mediaFilesRepo from '@sprouty-ai/shared/db/repositories/media-files';
import * as hotTopicsRepo from '@sprouty-ai/shared/db/repositories/hot-topics';
import { HotTopicService } from '@sprouty-ai/shared/services/hot-topic-service';
import { TopicRecommendService } from '@sprouty-ai/shared/services/topic-recommend-service';
import { generateProjectContext } from '@sprouty-ai/shared/db';
import { createBrowserProfileManager } from '@sprouty-ai/shared/services/browser-profile-manager';
import { createFingerprintGenerator } from '@sprouty-ai/shared/services/fingerprint-generator';
import { createAntiDetectBrowser } from '@sprouty-ai/shared/services/anti-detect-browser';
import type { Platform } from '@sprouty-ai/shared/db/types';

function getWorkspaceOrThrow(workspaceId: string): Workspace {
  const workspace = getWorkspaceByNameOrId(workspaceId);
  if (!workspace) throw new Error(`Workspace not found: ${workspaceId}`);
  return workspace;
}

export function registerCreatorMediaIpc(_windowManager: WindowManager): void {
  // ============================================================
  // 项目管理
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_LIST, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.listProjects(db);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_GET, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.getProject(db, projectId) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.createProject(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_UPDATE, async (_event, workspaceId: string, projectId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.updateProject(db, projectId, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_DELETE, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.deleteProject(db, projectId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_SET_ACTIVE, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.setActiveProject(db, projectId) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_GET_ACTIVE, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return projectsRepo.getActiveProject(db) ?? null;
  });

  // ============================================================
  // 账号画像
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROFILES_GET, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return profilesRepo.getProfileByProject(db, projectId) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PROFILES_UPSERT, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return profilesRepo.upsertProfile(db, data as any);
  });

  // ============================================================
  // 平台账号
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_LIST, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return platformAccountsRepo.listByProject(db, projectId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return platformAccountsRepo.createPlatformAccount(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return platformAccountsRepo.updatePlatformAccount(db, id, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return platformAccountsRepo.deletePlatformAccount(db, id);
  });

  // ============================================================
  // 竞品
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_LIST, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return competitorsRepo.listByProject(db, projectId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return competitorsRepo.createCompetitor(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return competitorsRepo.updateCompetitor(db, id, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return competitorsRepo.deleteCompetitor(db, id);
  });

  // ============================================================
  // 内容
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_LIST, async (_event, workspaceId: string, projectId: string, filters?: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.listByProject(db, projectId, filters as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_GET, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.getContent(db, contentId) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.createContent(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_UPDATE, async (_event, workspaceId: string, contentId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.updateContent(db, contentId, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_UPDATE_STATUS, async (_event, workspaceId: string, contentId: string, status: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.updateContentStatus(db, contentId, status as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_DELETE, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentsRepo.deleteContent(db, contentId);
  });

  // ============================================================
  // 发布记录
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_LIST, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishRecordsRepo.listByContent(db, contentId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_GET, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishRecordsRepo.getPublishRecord(db, id) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishRecordsRepo.createPublishRecord(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishRecordsRepo.updatePublishRecord(db, id, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishRecordsRepo.deletePublishRecord(db, id);
  });

  // ============================================================
  // 爆款模式
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_LIST, async (_event, workspaceId: string, filters?: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return viralPatternsRepo.listViralPatterns(db, filters as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_GET, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return viralPatternsRepo.getViralPattern(db, id) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return viralPatternsRepo.createViralPattern(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return viralPatternsRepo.updateViralPattern(db, id, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return viralPatternsRepo.deleteViralPattern(db, id);
  });

  // ============================================================
  // 项目上下文
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTEXT_GET, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return generateProjectContext(db, projectId);
  });

  // ============================================================
  // 采集调度任务
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_TASKS_LIST, async (_event, workspaceId: string, publishRecordId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return reviewTasksRepo.listByPublishRecord(db, publishRecordId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_TASKS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return reviewTasksRepo.createReviewTask(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_TASKS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return reviewTasksRepo.updateReviewTask(db, id, data as any) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_TASKS_CANCEL, async (_event, workspaceId: string, publishRecordId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return reviewTasksRepo.cancelByPublishRecord(db, publishRecordId);
  });

  // ============================================================
  // 内容版本管理
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENT_VERSIONS_LIST, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentVersionsRepo.listByContent(db, contentId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENT_VERSIONS_GET, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentVersionsRepo.getVersion(db, id) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENT_VERSIONS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentVersionsRepo.createVersion(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_CONTENT_VERSIONS_ROLLBACK, async (_event, workspaceId: string, contentId: string, versionNumber: number) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return contentVersionsRepo.rollbackToVersion(db, contentId, versionNumber);
  });

  // ============================================================
  // 发布队列
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_QUEUE_LIST, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishQueueRepo.listByContent(db, contentId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_QUEUE_ENQUEUE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishQueueRepo.enqueue(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_QUEUE_CANCEL, async (_event, workspaceId: string, contentId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishQueueRepo.cancelByContent(db, contentId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_QUEUE_GET_NEXT, async (_event, workspaceId: string, platformAccountId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return publishQueueRepo.getNextInQueue(db, platformAccountId) ?? null;
  });

  // ============================================================
  // 采集调度服务（生命周期由主进程管理，此处仅暴露状态查询）
  // ============================================================

  // 注意：ReviewScheduler 的 start/stop 需要 collector 注入，
  // 实际生命周期在主进程启动时管理，此处提供状态查询接口
  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_SCHEDULER_STATUS, async () => {
    // 由外部注入的 scheduler 实例管理，此处返回占位状态
    return { running: false };
  });

  // ============================================================
  // 浏览器 Profile 管理
  // ============================================================

  const fingerprintGenerator = createFingerprintGenerator();

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_BROWSER_LAUNCH_LOGIN, async (_event, workspaceId: string, platformAccountId: string, platform: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const profileManager = createBrowserProfileManager(ws.rootPath);
    const browser = createAntiDetectBrowser({ profileManager, fingerprintGenerator });

    try {
      const result = await browser.launchForLogin(platformAccountId, platform as Platform);

      if (result.success) {
        // 登录成功，更新 platform_accounts.auth_status + 账号资料
        const db = getCreatorMediaDB(ws.rootPath);
        const updateData: Record<string, unknown> = {
          auth_status: 'logged_in',
          last_login_at: new Date().toISOString(),
        };

        // 如果采集到了账号资料，一并更新
        if (result.profile) {
          const p = result.profile;
          if (p.nickname) updateData.nickname = p.nickname;
          if (p.avatar_url) updateData.avatar_url = p.avatar_url;
          if (p.bio) updateData.bio = p.bio;
          if (p.home_url) updateData.home_url = p.home_url;
          if (p.platform_uid) updateData.platform_uid = p.platform_uid;
          if (p.followers !== undefined) updateData.followers = p.followers;
          if (p.following !== undefined) updateData.following = p.following;
          if (p.total_likes !== undefined) updateData.total_likes = p.total_likes;
          if (p.total_favorites !== undefined) updateData.total_favorites = p.total_favorites;
          if (p.total_posts !== undefined) updateData.total_posts = p.total_posts;
          updateData.metrics_updated_at = new Date().toISOString();
        }

        platformAccountsRepo.updatePlatformAccount(db, platformAccountId, updateData as any);
      }

      return result;
    } catch (error: any) {
      return { success: false, error: error.message ?? String(error) };
    }
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_BROWSER_CHECK_AUTH, async (_event, workspaceId: string, platformAccountId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const profileManager = createBrowserProfileManager(ws.rootPath);

    if (!profileManager.profileExists(platformAccountId)) {
      return { loggedIn: false, error: 'Profile not found' };
    }

    const browser = createAntiDetectBrowser({ profileManager, fingerprintGenerator });

    try {
      // 启动无头浏览器检查登录态
      const context = await browser.launch({ platformAccountId, headless: true });
      const page = context.pages()[0] ?? await context.newPage();

      // 获取平台账号信息以确定平台
      const db = getCreatorMediaDB(ws.rootPath);
      const accounts = platformAccountsRepo.listByProject(db, '');
      const account = accounts.find(a => a.id === platformAccountId);

      if (!account?.platform) {
        await browser.close();
        return { loggedIn: false, error: 'Platform account not found' };
      }

      // 访问平台首页检查是否仍然登录
      const platformUrls: Record<string, string> = {
        xiaohongshu: 'https://creator.xiaohongshu.com/',
        douyin: 'https://creator.douyin.com/',
        wechat: 'https://mp.weixin.qq.com/',
        bilibili: 'https://member.bilibili.com/',
        zhihu: 'https://www.zhihu.com/creator',
        weibo: 'https://weibo.com/',
        x: 'https://x.com/home',
      };

      const checkUrl = platformUrls[account.platform] ?? '';
      if (!checkUrl) {
        await browser.close();
        return { loggedIn: false, error: `Unknown platform: ${account.platform}` };
      }

      await page.goto(checkUrl, { waitUntil: 'domcontentloaded', timeout: 15000 });
      const finalUrl = page.url();
      await browser.close();

      // 如果被重定向到登录页，说明登录态已失效
      const loginPatterns = ['login', 'signin', 'passport', 'sso'];
      const isLoginPage = loginPatterns.some(p => finalUrl.toLowerCase().includes(p));

      const loggedIn = !isLoginPage;

      // 更新数据库中的 auth_status
      platformAccountsRepo.updatePlatformAccount(db, platformAccountId, {
        auth_status: loggedIn ? 'logged_in' : 'expired',
        last_login_check: new Date().toISOString(),
      } as any);

      return { loggedIn };
    } catch (error: any) {
      await browser.close();
      return { loggedIn: false, error: error.message ?? String(error) };
    }
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_BROWSER_PROFILE_EXISTS, async (_event, workspaceId: string, platformAccountId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const profileManager = createBrowserProfileManager(ws.rootPath);
    return profileManager.profileExists(platformAccountId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_BROWSER_DELETE_PROFILE, async (_event, workspaceId: string, platformAccountId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const profileManager = createBrowserProfileManager(ws.rootPath);
    return profileManager.deleteProfile(platformAccountId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_BROWSER_GENERATE_FINGERPRINT, async (_event, workspaceId: string, platformAccountId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const profileManager = createBrowserProfileManager(ws.rootPath);
    const fingerprint = fingerprintGenerator.generate(platformAccountId);
    profileManager.saveFingerprint(platformAccountId, fingerprint);
    return fingerprint;
  });

  // ============================================================
  // 草稿管理
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_DRAFTS_LIST, async (_event, workspaceId: string, projectId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return draftsRepo.listByProject(db, projectId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_DRAFTS_GET, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return draftsRepo.getDraft(db, id);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_DRAFTS_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return draftsRepo.createDraft(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_DRAFTS_UPDATE, async (_event, workspaceId: string, id: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return draftsRepo.updateDraft(db, id, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_DRAFTS_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return draftsRepo.deleteDraft(db, id);
  });

  // ============================================================
  // 素材管理
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_MEDIA_FILES_LIST, async (_event, workspaceId: string, projectId: string, filters?: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return mediaFilesRepo.listByProject(db, projectId, filters as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_MEDIA_FILES_GET, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return mediaFilesRepo.getMediaFile(db, id);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_MEDIA_FILES_CREATE, async (_event, workspaceId: string, data: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return mediaFilesRepo.createMediaFile(db, data as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_MEDIA_FILES_DELETE, async (_event, workspaceId: string, id: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return mediaFilesRepo.deleteMediaFile(db, id);
  });

  // ============================================================
  // 热榜
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOT_TOPICS_FETCH, async (_event, workspaceId: string, platforms?: string[]) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new HotTopicService(db);
    return await service.fetchHotTopics(platforms);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOT_TOPICS_LIST, async (_event, workspaceId: string, filters?: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new HotTopicService(db);
    return service.getHotTopics(filters as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOT_TOPICS_GET_LATEST_BATCH, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new HotTopicService(db);
    return service.getLatestBatch() ?? null;
  });

  // ============================================================
  // 选题推荐
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_LIST, async (_event, workspaceId: string, projectId: string, filters?: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new TopicRecommendService(db);
    return service.getTopics(projectId, filters as any);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_GET, async (_event, workspaceId: string, topicId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new TopicRecommendService(db);
    return service.getTopic(topicId) ?? null;
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_ADOPT, async (_event, workspaceId: string, topicId: string, projectId: string, pipelineMode?: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new TopicRecommendService(db);
    return service.adoptTopic(topicId, projectId, (pipelineMode as any) ?? 'semi-auto');
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_IGNORE, async (_event, workspaceId: string, topicId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new TopicRecommendService(db);
    return service.ignoreTopic(topicId);
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_BATCH_IGNORE, async (_event, workspaceId: string, topicIds: string[]) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    const service = new TopicRecommendService(db);
    return service.batchIgnore(topicIds);
  });

  // ============================================================
  // 选题调度配置
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPIC_SCHEDULE_GET, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    // 从项目设置读取调度配置，默认早中晚三次
    const config = db.prepare<{ value: string }>(
      "SELECT value FROM settings WHERE key = 'topic_schedule'"
    ).get();
    if (config) return JSON.parse(config.value);
    return { hours: [8, 12, 20], autoGenerate: true };
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPIC_SCHEDULE_UPDATE, async (_event, workspaceId: string, config: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    db.prepare(
      "INSERT OR REPLACE INTO settings (key, value) VALUES ('topic_schedule', ?)"
    ).run(JSON.stringify(config));
  });

  // ============================================================
  // Hooks（读写 hooks.json）
  // ============================================================

  /** 展开 ~ 为用户主目录 */
  function resolveWorkspacePath(rootPath: string): string {
    if (rootPath.startsWith('~/') || rootPath === '~') {
      return path.join(os.homedir(), rootPath.slice(1));
    }
    return rootPath;
  }

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOOKS_READ, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const hooksPath = path.join(resolveWorkspacePath(ws.rootPath), '.sprouty-ai', 'hooks.json');
    try {
      const raw = fs.readFileSync(hooksPath, 'utf-8');
      return JSON.parse(raw);
    } catch {
      return { hooks: {} };
    }
  });

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOOKS_WRITE, async (_event, workspaceId: string, config: unknown) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const resolved = resolveWorkspacePath(ws.rootPath);
    const hooksDir = path.join(resolved, '.sprouty-ai');
    const hooksPath = path.join(hooksDir, 'hooks.json');
    // 确保目录存在
    if (!fs.existsSync(hooksDir)) {
      fs.mkdirSync(hooksDir, { recursive: true });
    }
    fs.writeFileSync(hooksPath, JSON.stringify(config, null, 2), 'utf-8');
    // ConfigWatcher 会自动检测文件变更并调用 hookSystem.reloadConfig()
    return { success: true };
  });

  // 采集调度任务 — 全量只读查看（定时任务视图用）
  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_REVIEW_TASKS_LIST_ALL, async (_event, workspaceId: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const db = getCreatorMediaDB(ws.rootPath);
    return db.prepare('SELECT * FROM review_tasks ORDER BY scheduled_at DESC').all();
  });

  // ============================================================
  // Hook 执行记录（从数据库读取）
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_HOOK_EVENTS_LIST, async (_event, workspaceId: string, options?: { limit?: number; taskId?: string }) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    const resolved = resolveWorkspacePath(ws.rootPath);
    const dbPath = path.join(resolved, '.sprouty-ai', 'db', 'creator.db');

    if (!fs.existsSync(dbPath)) return [];

    const limit = options?.limit ?? 50;
    const taskId = options?.taskId;

    try {
      const { getCachedConnection } = await import('@sprouty-ai/shared/db');
      const db = getCachedConnection(dbPath);
      if (!db) return [];

      let sql = `
        SELECT * FROM scheduled_task_executions
        ${taskId ? 'WHERE task_id = ?' : ''}
        ORDER BY trigger_time DESC
        LIMIT ?
      `;

      const params = taskId ? [taskId, limit] : [limit];
      const results = db.all(sql, params);

      return results;
    } catch (err) {
      ipcLog.warn('[creator-media-ipc] 读取执行记录失败:', err);
      return [];
    }
  });

  // ============================================================
  // 选题推荐 — 读取 Markdown 详情文件
  // ============================================================

  ipcMain.handle(IPC_CHANNELS.CREATOR_MEDIA_TOPICS_READ_MD, async (_event, workspaceId: string, mdFilePath: string) => {
    const ws = getWorkspaceOrThrow(workspaceId);
    // 安全校验：确保路径在工作区目录内，防止路径穿越
    const resolved = path.resolve(ws.rootPath, mdFilePath);
    if (!resolved.startsWith(path.resolve(ws.rootPath))) {
      throw new Error('路径不在工作区目录内');
    }
    try {
      return fs.readFileSync(resolved, 'utf-8');
    } catch {
      return null;
    }
  });

  ipcLog.info('[creator-media-ipc] 已注册所有 creatorMedia IPC 通道');
}
