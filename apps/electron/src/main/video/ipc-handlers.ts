/**
 * Video IPC Handlers
 *
 * Registers IPC handlers for video service management, project management,
 * asset management, rendering, and preview operations.
 *
 * @requirements 5.1, 5.2, 8.2, 8.3, 8.4, 8.5, 8.6
 */

import { ipcMain, dialog, BrowserWindow } from 'electron';
import { extname, join } from 'path';
import { app } from 'electron';
import type {
  VideoProject,
  RenderProgress,
  AssetType,
  Asset,
  OutputFormat,
  QualityPreset,
} from '@creator-flow/video';
import { SUPPORTED_ASSET_EXTENSIONS } from '@creator-flow/video';
import type { IProjectManager, CreateProjectOptions } from './project-manager';
import type { IRenderWorker } from './render-worker';
import type { IPreviewServer, PreviewServerResult } from './preview-server';
import type { VideoServiceManager, ServiceStatus } from './service-manager';
import { createVideoServiceManager } from './service-manager';
import { createPreviewServer } from './preview-server';
import { RenderWorker } from './render-worker';
import { ProjectManager } from './project-manager';
import { createBunPathResolver } from '../bun-path';
import log from '../logger';

const ipcLog = log.scope('video:ipc');

// ============================================================================
// IPC Channel Definitions
// ============================================================================

/**
 * Video IPC channel names
 * @requirements 5.1, 5.2, 8.2, 8.3, 8.4
 */
export const VIDEO_IPC_CHANNELS = {
  // Service management
  // @requirements 5.1, 5.2
  START_VIDEO_SERVICE: 'video:start-service',
  STOP_VIDEO_SERVICE: 'video:stop-service',
  GET_SERVICE_STATUS: 'video:get-status',
  SERVICE_STATUS_CHANGED: 'video:status-changed',

  // Project management
  CREATE_PROJECT: 'video:create-project',
  LIST_PROJECTS: 'video:list-projects',
  GET_PROJECT: 'video:get-project',
  UPDATE_PROJECT: 'video:update-project',
  DELETE_PROJECT: 'video:delete-project',

  // Asset management
  ADD_ASSET: 'video:add-asset',
  REMOVE_ASSET: 'video:remove-asset',

  // Rendering
  RENDER: 'video:render',
  CANCEL_RENDER: 'video:cancel-render',
  RENDER_PROGRESS: 'video:render-progress',
  RENDER_COMPLETED: 'video:render-completed',
  RENDER_FAILED: 'video:render-failed',

  // Preview
  START_PREVIEW: 'video:start-preview',
  STOP_PREVIEW: 'video:stop-preview',
  PREVIEW_STARTED: 'video:preview-started',
  PREVIEW_STOPPED: 'video:preview-stopped',
} as const;

// ============================================================================
// IPC Request/Response Types
// ============================================================================

/**
 * Service status response for IPC
 */
export interface ServiceStatusResponse {
  running: boolean;
  pid?: number;
  startedAt?: number;
  restartCount: number;
  lastError?: string;
}

/**
 * Create project request
 */
export interface CreateProjectRequest {
  name: string;
  workspaceId: string;
  template?: string;
  config?: {
    width?: number;
    height?: number;
    fps?: number;
    durationInFrames?: number;
  };
  description?: string;
}

/**
 * Render request
 */
export interface RenderRequest {
  projectId: string;
  compositionId: string;
  outputFormat?: OutputFormat;
  quality?: QualityPreset;
}

/**
 * Add asset request
 */
export interface AddAssetRequest {
  projectId: string;
  assetPath: string;
  assetType: AssetType;
}

// ============================================================================
// Global Service Instances
// ============================================================================

/** Global video service manager instance */
let videoServiceManager: VideoServiceManager | null = null;

/** Global preview server instance */
let previewServer: IPreviewServer | null = null;

/** Global render worker instance */
let renderWorker: IRenderWorker | null = null;

/** Global project manager instance */
let projectManager: IProjectManager | null = null;

/** Track registered windows for broadcasting events */
const registeredWindows = new Set<BrowserWindow>();

/**
 * Get or create the video service manager
 */
function getVideoServiceManager(): VideoServiceManager {
  if (!videoServiceManager) {
    const bunResolver = createBunPathResolver();
    const mcpServerEntry = app.isPackaged
      ? join(process.resourcesPath, 'app', 'apps', 'mcp-video', 'src', 'index.ts')
      : join(__dirname, '..', '..', '..', '..', 'mcp-video', 'src', 'index.ts');

    videoServiceManager = createVideoServiceManager(bunResolver, {
      mcpServerEntry,
      maxRetries: 3,
      retryInterval: 2000,
      startupTimeout: 30000,
      stopTimeout: 5000,
    });

    // Set up event handlers for broadcasting status changes
    videoServiceManager.on('status-changed', (status: ServiceStatus) => {
      broadcastToAllWindows(VIDEO_IPC_CHANNELS.SERVICE_STATUS_CHANGED, status);
    });

    videoServiceManager.on('error', (error: Error) => {
      ipcLog.error('Video service error:', error);
    });

    videoServiceManager.on('log', (level: 'info' | 'error', message: string) => {
      if (level === 'error') {
        ipcLog.error(message);
      } else {
        ipcLog.info(message);
      }
    });
  }
  return videoServiceManager;
}

/**
 * Get or create the preview server
 */
function getPreviewServer(): IPreviewServer {
  if (!previewServer) {
    const bunResolver = createBunPathResolver();
    previewServer = createPreviewServer(bunResolver);
  }
  return previewServer;
}

/**
 * Get or create the render worker
 */
function getRenderWorker(): IRenderWorker {
  if (!renderWorker) {
    const bunResolver = createBunPathResolver();
    renderWorker = new RenderWorker(bunResolver);
  }
  return renderWorker;
}

/**
 * Get or create the project manager
 */
function getProjectManager(): IProjectManager {
  if (!projectManager) {
    projectManager = new ProjectManager();
  }
  return projectManager;
}

/**
 * Broadcast a message to all registered windows
 */
function broadcastToAllWindows(channel: string, data: unknown): void {
  registeredWindows.forEach((window) => {
    if (!window.isDestroyed() && !window.webContents.isDestroyed()) {
      window.webContents.send(channel, data);
    }
  });
}

// ============================================================================
// IPC Handler Registration
// ============================================================================

/**
 * Register all video-related IPC handlers
 *
 * This function registers handlers for:
 * - Service management (start/stop/status)
 * - Project management (CRUD operations)
 * - Asset management (add/remove)
 * - Rendering (render/cancel/progress)
 * - Preview (start/stop)
 *
 * @requirements 5.1, 5.2, 8.2, 8.3, 8.4, 8.5, 8.6
 */
export function registerVideoIpcHandlers(): void {
  ipcLog.info('Registering video IPC handlers');

  // ============================================================================
  // Service Management IPC Handlers
  // @requirements 5.1, 5.2
  // ============================================================================

  /**
   * Start the video service (MCP Video Server)
   * @requirements 5.1 - 使用内置 Bun 启动 MCP Video Server 作为子进程
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.START_VIDEO_SERVICE,
    async (event): Promise<ServiceStatusResponse> => {
      ipcLog.info('Starting video service');

      try {
        // Register the window for broadcasts
        const browserWindow = BrowserWindow.fromWebContents(event.sender);
        if (browserWindow) {
          registeredWindows.add(browserWindow);
          browserWindow.on('closed', () => {
            registeredWindows.delete(browserWindow);
          });
        }

        const manager = getVideoServiceManager();
        await manager.startMcpServer();
        const status = manager.getStatus();
        ipcLog.info(`Video service started: PID=${status.pid}`);
        return status;
      } catch (error) {
        ipcLog.error('Failed to start video service:', error);
        throw error;
      }
    }
  );

  /**
   * Stop the video service
   * @requirements 5.1 - 管理 MCP Video Server 生命周期
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.STOP_VIDEO_SERVICE,
    async (): Promise<ServiceStatusResponse> => {
      ipcLog.info('Stopping video service');

      try {
        const manager = getVideoServiceManager();
        await manager.stopMcpServer();
        const status = manager.getStatus();
        ipcLog.info('Video service stopped');
        return status;
      } catch (error) {
        ipcLog.error('Failed to stop video service:', error);
        throw error;
      }
    }
  );

  /**
   * Get the current video service status
   * @requirements 5.2 - 通过 stdio 传输模式与 Electron 主进程通信
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.GET_SERVICE_STATUS,
    async (): Promise<ServiceStatusResponse> => {
      ipcLog.debug('Getting video service status');

      try {
        const manager = getVideoServiceManager();
        return manager.getStatus();
      } catch (error) {
        ipcLog.error('Failed to get video service status:', error);
        throw error;
      }
    }
  );

  // ============================================================================
  // Project Management IPC Handlers
  // @requirements 8.2
  // ============================================================================

  /**
   * Create a new video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.CREATE_PROJECT,
    async (_event, request: CreateProjectRequest): Promise<VideoProject> => {
      ipcLog.info(`Creating video project: ${request.name}`);

      try {
        const pm = getProjectManager();
        const options: CreateProjectOptions = {
          name: request.name,
          workspaceId: request.workspaceId,
          template: request.template,
          config: request.config,
          description: request.description,
        };

        const project = await pm.createProject(options);
        ipcLog.info(`Video project created: ${project.id}`);
        return project;
      } catch (error) {
        ipcLog.error('Failed to create video project:', error);
        throw error;
      }
    }
  );

  /**
   * List all video projects in a workspace
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.LIST_PROJECTS,
    async (_event, workspaceId: string): Promise<VideoProject[]> => {
      ipcLog.debug(`Listing video projects for workspace: ${workspaceId}`);

      try {
        const pm = getProjectManager();
        const projects = await pm.listProjects(workspaceId);
        return projects;
      } catch (error) {
        ipcLog.error('Failed to list video projects:', error);
        throw error;
      }
    }
  );

  /**
   * Get a single video project by ID
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.GET_PROJECT,
    async (_event, projectId: string): Promise<VideoProject | null> => {
      ipcLog.debug(`Getting video project: ${projectId}`);

      try {
        const pm = getProjectManager();
        const project = await pm.getProject(projectId);
        return project;
      } catch (error) {
        ipcLog.error('Failed to get video project:', error);
        throw error;
      }
    }
  );

  /**
   * Update a video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.UPDATE_PROJECT,
    async (
      _event,
      projectId: string,
      updates: Partial<VideoProject>
    ): Promise<VideoProject> => {
      ipcLog.info(`Updating video project: ${projectId}`);

      try {
        const pm = getProjectManager();
        const project = await pm.updateProject(projectId, updates);
        ipcLog.info(`Video project updated: ${projectId}`);
        return project;
      } catch (error) {
        ipcLog.error('Failed to update video project:', error);
        throw error;
      }
    }
  );

  /**
   * Delete a video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.DELETE_PROJECT,
    async (_event, projectId: string): Promise<boolean> => {
      ipcLog.info(`Deleting video project: ${projectId}`);

      try {
        const pm = getProjectManager();
        const success = await pm.deleteProject(projectId);
        if (success) {
          ipcLog.info(`Video project deleted: ${projectId}`);
        } else {
          ipcLog.warn(`Video project not found: ${projectId}`);
        }
        return success;
      } catch (error) {
        ipcLog.error('Failed to delete video project:', error);
        throw error;
      }
    }
  );

  // ============================================================================
  // Asset Management IPC Handlers
  // @requirements 8.3, 13.3
  // ============================================================================

  /**
   * Add an asset to a video project
   * Validates asset file existence and format before adding
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.ADD_ASSET,
    async (_event, request: AddAssetRequest): Promise<Asset> => {
      const { projectId, assetPath, assetType } = request;
      ipcLog.info(`Adding ${assetType} asset to project ${projectId}: ${assetPath}`);

      try {
        const pm = getProjectManager();
        // Validate asset type and extension
        const ext = extname(assetPath).toLowerCase();
        const supportedExts = SUPPORTED_ASSET_EXTENSIONS[assetType];

        if (!supportedExts.includes(ext as any)) {
          const errorMsg = `不支持的 ${assetType} 格式: ${ext}。支持的格式: ${supportedExts.join(', ')}`;
          ipcLog.error(errorMsg);
          throw new Error(errorMsg);
        }

        const asset = await pm.addAsset(projectId, assetPath, assetType);
        ipcLog.info(`Asset added: ${asset.id}`);
        return asset;
      } catch (error) {
        ipcLog.error('Failed to add asset:', error);
        throw error;
      }
    }
  );

  /**
   * Remove an asset from a video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.REMOVE_ASSET,
    async (_event, projectId: string, assetId: string): Promise<boolean> => {
      ipcLog.info(`Removing asset ${assetId} from project ${projectId}`);

      try {
        const pm = getProjectManager();
        const success = await pm.removeAsset(projectId, assetId);
        if (success) {
          ipcLog.info(`Asset removed: ${assetId}`);
        } else {
          ipcLog.warn(`Asset not found: ${assetId}`);
        }
        return success;
      } catch (error) {
        ipcLog.error('Failed to remove asset:', error);
        throw error;
      }
    }
  );

  // ============================================================================
  // Render IPC Handlers
  // @requirements 8.4, 8.5, 8.6
  // ============================================================================

  /**
   * Render a video composition
   * Shows save dialog for output path selection
   * Pushes progress events via video:render-progress channel
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.RENDER,
    async (event, request: RenderRequest): Promise<string | null> => {
      const { projectId, compositionId, outputFormat = 'mp4', quality = 'standard' } = request;
      ipcLog.info(`Starting render: project=${projectId}, composition=${compositionId}, format=${outputFormat}, quality=${quality}`);

      try {
        const pm = getProjectManager();
        const rw = getRenderWorker();

        // Get project to validate it exists
        const project = await pm.getProject(projectId);
        if (!project) {
          throw new Error(`项目不存在: ${projectId}`);
        }

        // Validate composition exists
        const composition = project.compositions.find((c) => c.id === compositionId);
        if (!composition) {
          throw new Error(`组件不存在: ${compositionId}`);
        }

        // Get the browser window for the dialog
        const browserWindow = BrowserWindow.fromWebContents(event.sender);

        // Show save dialog for output path selection
        // @requirements 8.5 - Show system save dialog
        const fileExtension = outputFormat;
        const defaultFileName = `${project.name}_${composition.name}.${fileExtension}`;

        const dialogResult = await dialog.showSaveDialog(browserWindow!, {
          title: '保存视频',
          defaultPath: defaultFileName,
          filters: [
            {
              name: getFormatDisplayName(outputFormat),
              extensions: [fileExtension],
            },
          ],
        });

        if (dialogResult.canceled || !dialogResult.filePath) {
          ipcLog.info('Render cancelled by user (save dialog)');
          return null;
        }

        const outputPath = dialogResult.filePath;
        ipcLog.info(`Output path selected: ${outputPath}`);

        // Get project path for rendering
        const projectPath = pm.getProjectPath(projectId);

        // Create progress callback that sends events to renderer
        // @requirements 8.6 - Push render progress via IPC
        const onProgress = (progress: RenderProgress) => {
          // Send progress event to the renderer process
          event.sender.send(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, {
            projectId,
            compositionId,
            ...progress,
          });
        };

        // Start rendering
        const result = await rw.render({
          projectPath,
          compositionId,
          outputPath,
          quality,
          outputFormat,
          onProgress,
        });

        // Send completion event
        event.sender.send(VIDEO_IPC_CHANNELS.RENDER_COMPLETED, {
          projectId,
          compositionId,
          outputPath: result,
        });

        ipcLog.info(`Render completed: ${result}`);
        return result;
      } catch (error) {
        ipcLog.error('Render failed:', error);
        // Send failure event
        event.sender.send(VIDEO_IPC_CHANNELS.RENDER_FAILED, {
          projectId,
          compositionId,
          error: error instanceof Error ? error.message : String(error),
        });
        throw error;
      }
    }
  );

  /**
   * Cancel the current render operation
   */
  ipcMain.handle(VIDEO_IPC_CHANNELS.CANCEL_RENDER, async (): Promise<void> => {
    ipcLog.info('Cancelling render');

    try {
      const rw = getRenderWorker();
      rw.cancel();
      ipcLog.info('Render cancelled');
    } catch (error) {
      ipcLog.error('Failed to cancel render:', error);
      throw error;
    }
  });

  // ============================================================================
  // Preview IPC Handlers
  // @requirements 8.4
  // ============================================================================

  /**
   * Start preview server for a video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.START_PREVIEW,
    async (event, projectId: string): Promise<PreviewServerResult> => {
      ipcLog.info(`Starting preview for project: ${projectId}`);

      try {
        const pm = getProjectManager();
        const ps = getPreviewServer();

        // Get project path
        const projectPath = pm.getProjectPath(projectId);

        // Start preview server
        const result = await ps.start(projectPath);
        ipcLog.info(`Preview started at: ${result.url}`);

        // Send preview started event
        event.sender.send(VIDEO_IPC_CHANNELS.PREVIEW_STARTED, {
          projectId,
          url: result.url,
          port: result.port,
          pid: result.pid,
        });

        return result;
      } catch (error) {
        ipcLog.error('Failed to start preview:', error);
        throw error;
      }
    }
  );

  /**
   * Stop preview server for a video project
   */
  ipcMain.handle(
    VIDEO_IPC_CHANNELS.STOP_PREVIEW,
    async (event, projectId: string): Promise<void> => {
      ipcLog.info(`Stopping preview for project: ${projectId}`);

      try {
        const pm = getProjectManager();
        const ps = getPreviewServer();

        // Get project path
        const projectPath = pm.getProjectPath(projectId);

        // Stop preview server
        await ps.stop(projectPath);
        ipcLog.info('Preview stopped');

        // Send preview stopped event
        event.sender.send(VIDEO_IPC_CHANNELS.PREVIEW_STOPPED, {
          projectId,
        });
      } catch (error) {
        ipcLog.error('Failed to stop preview:', error);
        throw error;
      }
    }
  );

  ipcLog.info('Video IPC handlers registered successfully');
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Get display name for output format
 */
function getFormatDisplayName(format: OutputFormat): string {
  switch (format) {
    case 'mp4':
      return 'MP4 视频 (H.264)';
    case 'webm':
      return 'WebM 视频 (VP8)';
    case 'gif':
      return 'GIF 动画';
    default:
      return 'Video';
  }
}

/**
 * Unregister all video IPC handlers
 * Call this when cleaning up the video service
 */
export function unregisterVideoIpcHandlers(): void {
  ipcLog.info('Unregistering video IPC handlers');

  Object.values(VIDEO_IPC_CHANNELS).forEach((channel) => {
    // Skip event channels (they don't have handlers)
    if (
      channel === VIDEO_IPC_CHANNELS.RENDER_PROGRESS ||
      channel === VIDEO_IPC_CHANNELS.RENDER_COMPLETED ||
      channel === VIDEO_IPC_CHANNELS.RENDER_FAILED ||
      channel === VIDEO_IPC_CHANNELS.SERVICE_STATUS_CHANGED ||
      channel === VIDEO_IPC_CHANNELS.PREVIEW_STARTED ||
      channel === VIDEO_IPC_CHANNELS.PREVIEW_STOPPED
    ) {
      return;
    }
    ipcMain.removeHandler(channel);
  });

  ipcLog.info('Video IPC handlers unregistered');
}

/**
 * Cleanup all video services
 * Call this when the application is shutting down
 */
export async function cleanupVideoServices(): Promise<void> {
  ipcLog.info('Cleaning up video services');

  try {
    // Stop video service manager
    if (videoServiceManager) {
      await videoServiceManager.stopAll();
      videoServiceManager = null;
    }

    // Stop all preview servers
    if (previewServer) {
      await previewServer.stopAll();
      previewServer = null;
    }

    // Cancel any ongoing renders
    if (renderWorker) {
      renderWorker.cancel();
      renderWorker = null;
    }

    // Clear project manager
    projectManager = null;

    // Clear registered windows
    registeredWindows.clear();

    ipcLog.info('Video services cleaned up');
  } catch (error) {
    ipcLog.error('Error cleaning up video services:', error);
  }
}

/**
 * Get the video service manager instance (for external use)
 */
export function getVideoServiceManagerInstance(): VideoServiceManager | null {
  return videoServiceManager;
}
