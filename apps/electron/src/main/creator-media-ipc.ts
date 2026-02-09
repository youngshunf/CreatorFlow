/**
 * Electron 主进程 — 自媒体创作 IPC 通道注册
 *
 * 所有 creatorMedia:* IPC handler
 */

import { ipcMain } from 'electron';
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
import { generateProjectContext } from '@sprouty-ai/shared/db';

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

  ipcLog.info('[creator-media-ipc] 已注册所有 creatorMedia IPC 通道');
}
