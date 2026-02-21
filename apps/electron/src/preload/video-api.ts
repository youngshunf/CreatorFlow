/**
 * Video API for Preload
 *
 * 暴露视频相关 IPC 调用到渲染进程。
 * 已移除 VideoServiceAPI（MCP Server 生命周期管理）和预览相关接口。
 */

import { ipcRenderer } from "electron";
import type {
  VideoProject,
  VideoProjectFull,
  VideoScene,
  VideoAsset,
  UpdateVideoProject,
  UpdateVideoScene,
  AssetType,
} from "@sprouty-ai/shared/db/types";
import type { VideoTemplate } from "@sprouty-ai/video/templates";
import type { OutputFormat, QualityPreset, RenderProgress } from "@sprouty-ai/video";

// ============================================================================
// IPC Channel Definitions (must match main process)
// ============================================================================

const VIDEO_IPC_CHANNELS = {
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
// Request Types
// ============================================================================

export interface CreateProjectRequest {
  workspaceId: string;
  projectId: string;
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
  projectId: string;
  contentId: string;
  templateId: string;
  name: string;
  description?: string;
}

export interface CompositionInfo {
  id: string;
  name: string;
}

export interface RenderProgressData extends RenderProgress {
  projectId: string;
}

export interface RenderCompletedData {
  projectId: string;
  outputPath: string;
}

export interface RenderFailedData {
  projectId: string;
  error: string;
}

// ============================================================================
// VideoAPI Interface
// ============================================================================

export interface VideoAPI {
  // 项目
  createProject(request: CreateProjectRequest): Promise<VideoProjectFull>;
  getProject(workspaceId: string, projectId: string): Promise<VideoProjectFull | null>;
  listProjects(workspaceId: string, projectId: string): Promise<VideoProject[]>;
  updateProject(workspaceId: string, projectId: string, updates: UpdateVideoProject): Promise<VideoProject | null>;
  deleteProject(workspaceId: string, projectId: string): Promise<boolean>;

  // 场景
  addScene(request: AddSceneRequest): Promise<{ sceneId: string }>;
  updateScene(request: UpdateSceneRequest): Promise<VideoScene | null>;
  removeScene(workspaceId: string, projectId: string, sceneId: string): Promise<boolean>;
  reorderScenes(workspaceId: string, projectId: string, sceneIds: string[]): Promise<void>;
  getScenes(workspaceId: string, projectId: string): Promise<VideoScene[]>;

  // 素材
  selectAssetFiles(filterType?: AssetType): Promise<{ filePath: string; assetType: AssetType }[]>;
  addAsset(request: AddAssetRequest): Promise<VideoAsset>;
  removeAsset(workspaceId: string, projectId: string, assetId: string): Promise<boolean>;
  listAssets(workspaceId: string, projectId: string): Promise<VideoAsset[]>;

  // 组合和模板
  listCompositions(): Promise<CompositionInfo[]>;
  listTemplates(category?: string): Promise<VideoTemplate[]>;
  getTemplate(templateId: string): Promise<VideoTemplate | null>;
  createFromTemplate(request: CreateFromTemplateRequest): Promise<VideoProjectFull>;

  // 渲染
  render(request: RenderRequest): Promise<string | null>;
  cancelRender(): void;
  onRenderProgress(callback: (data: RenderProgressData) => void): () => void;
  onRenderCompleted(callback: (data: RenderCompletedData) => void): () => void;
  onRenderFailed(callback: (data: RenderFailedData) => void): () => void;
}

// ============================================================================
// VideoAPI Implementation
// ============================================================================

export function createVideoAPI(): VideoAPI {
  return {
    // 项目
    createProject: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CREATE_PROJECT, request),

    getProject: (workspaceId, projectId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.GET_PROJECT, workspaceId, projectId),

    listProjects: (workspaceId, projectId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.LIST_PROJECTS, workspaceId, projectId),

    updateProject: (workspaceId, projectId, updates) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.UPDATE_PROJECT, workspaceId, projectId, updates),

    deleteProject: (workspaceId, projectId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.DELETE_PROJECT, workspaceId, projectId),

    // 场景
    addScene: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.ADD_SCENE, request),

    updateScene: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.UPDATE_SCENE, request),

    removeScene: (workspaceId, projectId, sceneId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.REMOVE_SCENE, workspaceId, projectId, sceneId),

    reorderScenes: (workspaceId, projectId, sceneIds) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.REORDER_SCENES, workspaceId, projectId, sceneIds),

    getScenes: (workspaceId, projectId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.GET_SCENES, workspaceId, projectId),

    // 素材
    selectAssetFiles: (filterType?) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.SELECT_ASSET_FILES, filterType),

    addAsset: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.ADD_ASSET, request),

    removeAsset: (workspaceId, projectId, assetId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.REMOVE_ASSET, workspaceId, projectId, assetId),

    listAssets: (workspaceId, projectId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.LIST_ASSETS, workspaceId, projectId),

    // 组合和模板
    listCompositions: () =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.LIST_COMPOSITIONS),

    listTemplates: (category?) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.LIST_TEMPLATES, category),

    getTemplate: (templateId) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.GET_TEMPLATE, templateId),

    createFromTemplate: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CREATE_FROM_TEMPLATE, request),

    // 渲染
    render: (request) =>
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.RENDER, request),

    cancelRender: () => {
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CANCEL_RENDER);
    },

    onRenderProgress: (callback) => {
      const handler = (_event: Electron.IpcRendererEvent, data: RenderProgressData) => {
        callback(data);
      };
      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);
      };
    },

    onRenderCompleted: (callback) => {
      const handler = (_event: Electron.IpcRendererEvent, data: RenderCompletedData) => {
        callback(data);
      };
      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_COMPLETED, handler);
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_COMPLETED, handler);
      };
    },

    onRenderFailed: (callback) => {
      const handler = (_event: Electron.IpcRendererEvent, data: RenderFailedData) => {
        callback(data);
      };
      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_FAILED, handler);
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_FAILED, handler);
      };
    },
  };
}

export const videoAPI = createVideoAPI();
