/**
 * Video IPC Handlers
 *
 * 视频相关 IPC 通道注册 — 直接使用 VideoRepository 操作 SQLite，
 * 不再依赖 MCP Server / ProjectManager。
 */

import { ipcMain, dialog, BrowserWindow } from "electron";
import { extname, basename } from "path";
import { statSync } from "fs";
import type {
  RenderProgress,
  AssetType,
  OutputFormat,
  QualityPreset,
} from "@sprouty-ai/video";
import { SUPPORTED_ASSET_EXTENSIONS } from "@sprouty-ai/video";
import { COMPOSITION_IDS } from "@sprouty-ai/video/compositions";
import {
  ALL_TEMPLATES,
  getTemplateById,
  getTemplatesByCategory,
} from "@sprouty-ai/video/templates";
import type { VideoTemplate, TemplateCategory } from "@sprouty-ai/video/templates";
import type { IRenderWorker } from "./render-worker";
import { RenderWorker } from "./render-worker";
import { createBunPathResolver } from "../bun-path";
import { getCreatorMediaDB } from "../creator-media-db";
import { getWorkspaceByNameOrId } from "@sprouty-ai/shared/config";
import * as videoRepo from "@sprouty-ai/shared/db/repositories/video";
import type {
  VideoProject,
  VideoProjectFull,
  VideoScene,
  VideoAsset,
  CreateVideoProject,
  CreateVideoScene,
  UpdateVideoScene,
  UpdateVideoProject,
  CreateVideoAsset,
} from "@sprouty-ai/shared/db/types";
import log from "../logger";

const ipcLog = log.scope("video:ipc");

// ============================================================================
// IPC Channel Definitions
// ============================================================================

export const VIDEO_IPC_CHANNELS = {
  // 项目管理
  CREATE_PROJECT: "video:create-project",
  LIST_PROJECTS: "video:list-projects",
  GET_PROJECT: "video:get-project",
  UPDATE_PROJECT: "video:update-project",
  DELETE_PROJECT: "video:delete-project",

  // 场景管理
  ADD_SCENE: "video:add-scene",
  UPDATE_SCENE: "video:update-scene",
  REMOVE_SCENE: "video:remove-scene",
  REORDER_SCENES: "video:reorder-scenes",
  GET_SCENES: "video:get-scenes",

  // 素材管理
  SELECT_ASSET_FILES: "video:select-asset-files",
  ADD_ASSET: "video:add-asset",
  REMOVE_ASSET: "video:remove-asset",
  LIST_ASSETS: "video:list-assets",

  // 组合和模板
  LIST_COMPOSITIONS: "video:list-compositions",
  LIST_TEMPLATES: "video:list-templates",
  GET_TEMPLATE: "video:get-template",
  CREATE_FROM_TEMPLATE: "video:create-from-template",

  // 渲染
  RENDER: "video:render",
  CANCEL_RENDER: "video:cancel-render",
  RENDER_PROGRESS: "video:render-progress",
  RENDER_COMPLETED: "video:render-completed",
  RENDER_FAILED: "video:render-failed",
} as const;

// ============================================================================
// IPC Request Types
// ============================================================================

export interface CreateProjectRequest {
  workspaceId: string;
  contentId: string;
  name: string;
  description?: string;
  width?: number;
  height?: number;
  fps?: number;
}

export interface AddSceneRequest {
  workspaceId: string;
  projectId: string;
  compositionId: string;
  name?: string;
  durationInFrames?: number;
  props?: string;
  transitionType?: string;
  transitionDuration?: number;
  transitionDirection?: string;
  insertAt?: number;
}

export interface UpdateSceneRequest {
  workspaceId: string;
  sceneId: string;
  updates: UpdateVideoScene;
}

export interface AddAssetRequest {
  workspaceId: string;
  projectId: string;
  filePath: string;
  type: AssetType;
  name?: string;
}

export interface RenderRequest {
  projectId: string;
  workspaceId: string;
  outputFormat?: OutputFormat;
  quality?: QualityPreset;
}

export interface CreateFromTemplateRequest {
  workspaceId: string;
  contentId: string;
  templateId: string;
  name: string;
  description?: string;
}

export interface CompositionInfo {
  id: string;
  name: string;
}

// ============================================================================
// Helper
// ============================================================================

function getDB(workspaceId: string) {
  const ws = getWorkspaceByNameOrId(workspaceId);
  if (!ws) throw new Error(`Workspace not found: ${workspaceId}`);
  return getCreatorMediaDB(ws.rootPath);
}

function generateId(): string {
  return crypto.randomUUID();
}

// ============================================================================
// Global Instances
// ============================================================================

let renderWorker: IRenderWorker | null = null;

function getRenderWorker(): IRenderWorker {
  if (!renderWorker) {
    const bunResolver = createBunPathResolver();
    renderWorker = new RenderWorker(bunResolver);
  }
  return renderWorker;
}

// ============================================================================
// IPC Handler Registration
// ============================================================================

export function registerVideoIpcHandlers(): void {
  ipcLog.info("注册视频 IPC handlers");

  // ============================================================================
  // 项目管理
  // ============================================================================

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.CREATE_PROJECT,
    async (_event, request: CreateProjectRequest): Promise<VideoProjectFull> => {
      ipcLog.info(`创建视频项目: ${request.name}`);
      const db = getDB(request.workspaceId);
      const data: CreateVideoProject = {
        id: generateId(),
        content_id: request.contentId,
        name: request.name,
        description: request.description,
        width: request.width ?? 1080,
        height: request.height ?? 1920,
        fps: request.fps ?? 30,
        metadata: undefined,
      };
      const project = videoRepo.createVideoProject(db, data);
      return { ...project, scenes: [], assets: [] };
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.GET_PROJECT,
    async (_event, workspaceId: string, projectId: string): Promise<VideoProjectFull | null> => {
      const db = getDB(workspaceId);
      return videoRepo.getVideoProjectFull(db, projectId) ?? null;
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.LIST_PROJECTS,
    async (_event, workspaceId: string, projectId: string): Promise<VideoProject[]> => {
      const db = getDB(workspaceId);
      return videoRepo.listVideoProjects(db, projectId);
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.UPDATE_PROJECT,
    async (_event, workspaceId: string, projectId: string, updates: UpdateVideoProject): Promise<VideoProject | null> => {
      ipcLog.info(`更新视频项目: ${projectId}`);
      const db = getDB(workspaceId);
      return videoRepo.updateVideoProject(db, projectId, updates) ?? null;
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.DELETE_PROJECT,
    async (_event, workspaceId: string, projectId: string): Promise<boolean> => {
      ipcLog.info(`删除视频项目: ${projectId}`);
      const db = getDB(workspaceId);
      return videoRepo.deleteVideoProject(db, projectId);
    },
  );

  // ============================================================================
  // 场景管理
  // ============================================================================

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.ADD_SCENE,
    async (_event, request: AddSceneRequest): Promise<{ sceneId: string }> => {
      ipcLog.info(`添加场景到项目: ${request.projectId}`);
      const db = getDB(request.workspaceId);

      // 计算 sort_order
      let sortOrder: number;
      if (request.insertAt !== undefined) {
        // 插入到指定位置，先把后面的场景 sort_order +1
        const scenes = videoRepo.listVideoScenes(db, request.projectId);
        db.transaction(() => {
          for (const s of scenes) {
            if (s.sort_order >= request.insertAt!) {
              videoRepo.updateVideoScene(db, s.id, { sort_order: s.sort_order + 1 });
            }
          }
        });
        sortOrder = request.insertAt;
      } else {
        sortOrder = videoRepo.getNextSortOrder(db, request.projectId);
      }

      const data: CreateVideoScene = {
        id: generateId(),
        project_id: request.projectId,
        composition_id: request.compositionId,
        name: request.name,
        sort_order: sortOrder,
        duration_in_frames: request.durationInFrames ?? 90,
        props: request.props ?? "{}",
        transition_type: (request.transitionType as CreateVideoScene["transition_type"]) ?? "none",
        transition_duration: request.transitionDuration ?? 0,
        transition_direction: request.transitionDirection as CreateVideoScene["transition_direction"],
      };

      const scene = videoRepo.addVideoScene(db, data);
      return { sceneId: scene.id };
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.UPDATE_SCENE,
    async (_event, request: UpdateSceneRequest): Promise<VideoScene | null> => {
      ipcLog.info(`更新场景: ${request.sceneId}`);
      const db = getDB(request.workspaceId);
      return videoRepo.updateVideoScene(db, request.sceneId, request.updates) ?? null;
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.REMOVE_SCENE,
    async (_event, workspaceId: string, projectId: string, sceneId: string): Promise<boolean> => {
      ipcLog.info(`删除场景: ${sceneId}`);
      const db = getDB(workspaceId);
      return videoRepo.removeVideoScene(db, sceneId);
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.REORDER_SCENES,
    async (_event, workspaceId: string, projectId: string, sceneIds: string[]): Promise<void> => {
      ipcLog.info(`重排场景: ${projectId}`);
      const db = getDB(workspaceId);
      videoRepo.reorderVideoScenes(db, projectId, sceneIds);
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.GET_SCENES,
    async (_event, workspaceId: string, projectId: string): Promise<VideoScene[]> => {
      const db = getDB(workspaceId);
      return videoRepo.listVideoScenes(db, projectId);
    },
  );

  // ============================================================================
  // 素材管理
  // ============================================================================

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.SELECT_ASSET_FILES,
    async (
      _event,
      filterType?: AssetType,
    ): Promise<{ filePath: string; assetType: AssetType }[]> => {
      const browserWindow = BrowserWindow.getFocusedWindow();
      const filters: Electron.FileDialogFilter[] = [];

      if (!filterType || filterType === "image") {
        filters.push({
          name: "图片",
          extensions: ["png", "jpg", "jpeg", "gif", "webp", "svg"],
        });
      }
      if (!filterType || filterType === "video") {
        filters.push({
          name: "视频",
          extensions: ["mp4", "webm", "mov"],
        });
      }
      if (!filterType || filterType === "audio") {
        filters.push({
          name: "音频",
          extensions: ["mp3", "wav", "ogg", "m4a"],
        });
      }
      if (!filterType || filterType === "font") {
        filters.push({
          name: "字体",
          extensions: ["ttf", "otf", "woff", "woff2"],
        });
      }
      if (!filterType) {
        filters.push({ name: "所有文件", extensions: ["*"] });
      }

      const result = await dialog.showOpenDialog(browserWindow!, {
        title: "选择素材文件",
        properties: ["openFile", "multiSelections"],
        filters,
      });

      if (result.canceled || result.filePaths.length === 0) return [];

      return result.filePaths.map((filePath) => {
        const ext = extname(filePath).toLowerCase();
        let detectedType: AssetType = "image";
        if (SUPPORTED_ASSET_EXTENSIONS.video.includes(ext as any)) {
          detectedType = "video";
        } else if (SUPPORTED_ASSET_EXTENSIONS.audio.includes(ext as any)) {
          detectedType = "audio";
        } else if (SUPPORTED_ASSET_EXTENSIONS.font.includes(ext as any)) {
          detectedType = "font";
        }
        return { filePath, assetType: filterType ?? detectedType };
      });
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.ADD_ASSET,
    async (_event, request: AddAssetRequest): Promise<VideoAsset> => {
      ipcLog.info(`添加素材到项目 ${request.projectId}: ${request.filePath}`);
      const db = getDB(request.workspaceId);

      // 验证素材类型
      const ext = extname(request.filePath).toLowerCase();
      const supportedExts = SUPPORTED_ASSET_EXTENSIONS[request.type];
      if (!supportedExts.includes(ext as any)) {
        throw new Error(
          `不支持的 ${request.type} 格式: ${ext}。支持的格式: ${supportedExts.join(", ")}`,
        );
      }

      let fileSize: number | undefined;
      try {
        fileSize = statSync(request.filePath).size;
      } catch {
        // 忽略
      }

      const data: CreateVideoAsset = {
        id: generateId(),
        project_id: request.projectId,
        type: request.type,
        name: request.name ?? basename(request.filePath),
        file_path: request.filePath,
        file_size: fileSize,
        metadata: undefined,
      };

      return videoRepo.addVideoAsset(db, data);
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.REMOVE_ASSET,
    async (_event, workspaceId: string, projectId: string, assetId: string): Promise<boolean> => {
      ipcLog.info(`删除素材: ${assetId}`);
      const db = getDB(workspaceId);
      return videoRepo.removeVideoAsset(db, assetId);
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.LIST_ASSETS,
    async (_event, workspaceId: string, projectId: string): Promise<VideoAsset[]> => {
      const db = getDB(workspaceId);
      return videoRepo.listVideoAssets(db, projectId);
    },
  );

  // ============================================================================
  // 组合和模板
  // ============================================================================

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.LIST_COMPOSITIONS,
    async (): Promise<CompositionInfo[]> => {
      return COMPOSITION_IDS.map((id) => ({ id, name: id }));
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.LIST_TEMPLATES,
    async (_event, category?: string): Promise<VideoTemplate[]> => {
      if (category) {
        return [...getTemplatesByCategory(category as TemplateCategory)];
      }
      return [...ALL_TEMPLATES];
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.GET_TEMPLATE,
    async (_event, templateId: string): Promise<VideoTemplate | null> => {
      return getTemplateById(templateId) ?? null;
    },
  );

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.CREATE_FROM_TEMPLATE,
    async (_event, request: CreateFromTemplateRequest): Promise<VideoProjectFull> => {
      ipcLog.info(`从模板创建项目: ${request.templateId}`);
      const template = getTemplateById(request.templateId);
      if (!template) throw new Error(`模板不存在: ${request.templateId}`);

      const db = getDB(request.workspaceId);

      // 创建项目
      const projectData: CreateVideoProject = {
        id: generateId(),
        content_id: request.contentId,
        name: request.name,
        description: request.description ?? template.description,
        width: template.defaultConfig.width,
        height: template.defaultConfig.height,
        fps: template.defaultConfig.fps,
        metadata: JSON.stringify({ templateId: template.id }),
      };
      const project = videoRepo.createVideoProject(db, projectData);

      // 创建默认场景
      const sceneData: CreateVideoScene = {
        id: generateId(),
        project_id: project.id,
        composition_id: template.compositionId,
        name: template.name,
        sort_order: 0,
        duration_in_frames: template.defaultConfig.durationInFrames,
        props: JSON.stringify(template.defaultProps),
        transition_type: "none",
        transition_duration: 0,
        transition_direction: undefined,
      };
      videoRepo.addVideoScene(db, sceneData);

      return videoRepo.getVideoProjectFull(db, project.id)!;
    },
  );

  // ============================================================================
  // 渲染（保留现有逻辑）
  // ============================================================================

  ipcMain.handle(
    VIDEO_IPC_CHANNELS.RENDER,
    async (event, request: RenderRequest): Promise<string | null> => {
      const { projectId, workspaceId, outputFormat = "mp4", quality = "standard" } = request;
      ipcLog.info(
        `开始渲染: project=${projectId}, format=${outputFormat}, quality=${quality}`,
      );

      try {
        const db = getDB(workspaceId);
        const rw = getRenderWorker();

        const project = videoRepo.getVideoProjectFull(db, projectId);
        if (!project) throw new Error(`项目不存在: ${projectId}`);
        if (!project.scenes || project.scenes.length === 0) {
          throw new Error("项目没有场景，无法渲染");
        }

        const browserWindow = BrowserWindow.fromWebContents(event.sender);

        const fileExtension = outputFormat;
        const defaultFileName = `${project.name}.${fileExtension}`;

        const dialogResult = await dialog.showSaveDialog(browserWindow!, {
          title: "保存视频",
          defaultPath: defaultFileName,
          filters: [
            {
              name: getFormatDisplayName(outputFormat),
              extensions: [fileExtension],
            },
          ],
        });

        if (dialogResult.canceled || !dialogResult.filePath) {
          ipcLog.info("用户取消保存");
          return null;
        }

        const outputPath = dialogResult.filePath;

        const onProgress = (progress: RenderProgress) => {
          event.sender.send(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, {
            projectId,
            ...progress,
          });
        };

        // 将 DB 场景转换为渲染所需格式
        const scenes = project.scenes.map((s) => ({
          compositionId: s.composition_id,
          durationInFrames: s.duration_in_frames,
          props: JSON.parse(s.props),
        }));

        const transitions = project.scenes
          .filter((s) => s.transition_type && s.transition_type !== "none")
          .map((s) => ({
            type: s.transition_type!,
            durationInFrames: s.transition_duration ?? 0,
            direction: s.transition_direction,
          }));

        const result = await rw.render({
          outputPath,
          quality,
          outputFormat,
          inputProps: { scenes, transitions } as any,
          onProgress,
        });

        event.sender.send(VIDEO_IPC_CHANNELS.RENDER_COMPLETED, {
          projectId,
          outputPath: result,
        });

        ipcLog.info(`渲染完成: ${result}`);
        return result;
      } catch (error) {
        ipcLog.error("渲染失败:", error);
        event.sender.send(VIDEO_IPC_CHANNELS.RENDER_FAILED, {
          projectId,
          error: error instanceof Error ? error.message : String(error),
        });
        throw error;
      }
    },
  );

  ipcMain.handle(VIDEO_IPC_CHANNELS.CANCEL_RENDER, async (): Promise<void> => {
    ipcLog.info("取消渲染");
    const rw = getRenderWorker();
    rw.cancel();
  });

  ipcLog.info("视频 IPC handlers 注册完成");
}

// ============================================================================
// Helper Functions
// ============================================================================

function getFormatDisplayName(format: OutputFormat): string {
  switch (format) {
    case "mp4":
      return "MP4 视频 (H.264)";
    case "webm":
      return "WebM 视频 (VP8)";
    case "gif":
      return "GIF 动画";
    default:
      return "Video";
  }
}

/**
 * 注销所有视频 IPC handlers
 */
export function unregisterVideoIpcHandlers(): void {
  ipcLog.info("注销视频 IPC handlers");
  Object.values(VIDEO_IPC_CHANNELS).forEach((channel) => {
    if (
      channel === VIDEO_IPC_CHANNELS.RENDER_PROGRESS ||
      channel === VIDEO_IPC_CHANNELS.RENDER_COMPLETED ||
      channel === VIDEO_IPC_CHANNELS.RENDER_FAILED
    ) {
      return;
    }
    ipcMain.removeHandler(channel);
  });
}

/**
 * 清理视频服务
 */
export async function cleanupVideoServices(): Promise<void> {
  ipcLog.info("清理视频服务");
  if (renderWorker) {
    renderWorker.cancel();
    renderWorker = null;
  }
}
