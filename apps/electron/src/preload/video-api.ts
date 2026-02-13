/**
 * Video API for Preload
 *
 * Exposes video-related IPC calls to the renderer process via contextBridge.
 * Provides type-safe VideoAPI interface for video project management,
 * asset management, rendering, and preview operations.
 *
 * Also provides VideoServiceAPI for video service lifecycle management
 * (start/stop/status) using the internal Bun runtime.
 *
 * @requirements 5.1, 15.1, 15.2, 15.3, 15.4
 */

import { ipcRenderer } from 'electron';
import type {
  VideoProject,
  Asset,
  AssetType,
  RenderProgress,
  OutputFormat,
  QualityPreset,
} from '@sprouty-ai/video';

// ============================================================================
// IPC Channel Definitions (must match main process)
// ============================================================================

/**
 * Video IPC channel names
 * These must match the channels defined in main/video/ipc-handlers.ts
 */
const VIDEO_IPC_CHANNELS = {
  // Service management
  // @requirements 5.1
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
// Request Types
// ============================================================================

/**
 * Service status response from video service manager
 * @requirements 5.1
 */
export interface ServiceStatusResponse {
  /** Whether the service is running */
  running: boolean;
  /** Process ID if running */
  pid?: number;
  /** Timestamp when service started */
  startedAt?: number;
  /** Number of times service has been restarted */
  restartCount: number;
  /** Last error message if any */
  lastError?: string;
}

/**
 * Render completed event data
 */
export interface RenderCompletedEvent {
  /** Project ID that was rendered */
  projectId: string;
  /** Composition ID that was rendered */
  compositionId: string;
  /** Output file path */
  outputPath: string;
}

/**
 * Render failed event data
 */
export interface RenderFailedEvent {
  /** Project ID that failed */
  projectId: string;
  /** Composition ID that failed */
  compositionId: string;
  /** Error message */
  error: string;
}

/**
 * Options for creating a new video project
 */
export interface CreateProjectOptions {
  /** Project name */
  name: string;
  /** Workspace ID where the project will be created */
  workspaceId: string;
  /** Optional template ID to use */
  template?: string;
  /** Optional video configuration overrides */
  config?: {
    width?: number;
    height?: number;
    fps?: number;
    durationInFrames?: number;
  };
  /** Optional project description */
  description?: string;
}

/**
 * Options for rendering a video
 */
export interface RenderOptions {
  /** Project ID to render */
  projectId: string;
  /** Composition ID to render */
  compositionId: string;
  /** Output format (default: 'mp4') */
  outputFormat?: OutputFormat;
  /** Quality preset (default: 'standard') */
  quality?: QualityPreset;
}

/**
 * Extended render progress with project context
 */
export interface RenderProgressEvent extends RenderProgress {
  /** Project ID being rendered */
  projectId: string;
  /** Composition ID being rendered */
  compositionId: string;
}

// ============================================================================
// VideoServiceAPI Interface
// ============================================================================

/**
 * Video Service API interface for service lifecycle management
 * Provides methods to start/stop the video service and monitor its status.
 * @requirements 5.1
 */
export interface VideoServiceAPI {
  /**
   * Start the video service (MCP Video Server)
   * @returns Service status after starting
   */
  startService(): Promise<ServiceStatusResponse>;

  /**
   * Stop the video service
   * @returns Service status after stopping
   */
  stopService(): Promise<ServiceStatusResponse>;

  /**
   * Get current service status
   * @returns Current service status
   */
  getStatus(): Promise<ServiceStatusResponse>;

  /**
   * Start preview server for a video project
   * @param projectId - Project ID to preview
   * @returns Preview URL and port
   */
  startPreview(projectId: string): Promise<{ url: string; port: number }>;

  /**
   * Stop preview server for a video project
   * @param projectId - Project ID to stop preview for
   */
  stopPreview(projectId: string): Promise<void>;

  /**
   * Render a video composition
   * @param request - Render options
   * @returns Output file path or null if cancelled
   */
  render(request: RenderOptions): Promise<string | null>;

  /**
   * Cancel the current render operation
   */
  cancelRender(): void;

  /**
   * Subscribe to service status change events
   * @param callback - Callback function for status updates
   * @returns Cleanup function to unsubscribe
   */
  onStatusChanged(callback: (status: ServiceStatusResponse) => void): () => void;

  /**
   * Subscribe to render progress events
   * @param callback - Callback function for progress updates
   * @returns Cleanup function to unsubscribe
   */
  onRenderProgress(callback: (progress: RenderProgressEvent) => void): () => void;

  /**
   * Subscribe to render completed events
   * @param callback - Callback function for completion
   * @returns Cleanup function to unsubscribe
   */
  onRenderCompleted(callback: (event: RenderCompletedEvent) => void): () => void;

  /**
   * Subscribe to render failed events
   * @param callback - Callback function for failures
   * @returns Cleanup function to unsubscribe
   */
  onRenderFailed(callback: (event: RenderFailedEvent) => void): () => void;
}

// ============================================================================
// VideoAPI Interface
// ============================================================================

/**
 * Video API interface exposed to renderer process
 * @requirements 15.2, 15.3
 */
export interface VideoAPI {
  // ============================================================================
  // Project Management
  // ============================================================================

  /**
   * Create a new video project
   * @param options - Project creation options
   * @returns The created video project
   */
  createProject(options: CreateProjectOptions): Promise<VideoProject>;

  /**
   * List all video projects in a workspace
   * @param workspaceId - Workspace ID to list projects from
   * @returns Array of video projects
   */
  listProjects(workspaceId: string): Promise<VideoProject[]>;

  /**
   * Get a single video project by ID
   * @param projectId - Project ID to retrieve
   * @returns The video project or null if not found
   */
  getProject(projectId: string): Promise<VideoProject | null>;

  /**
   * Update a video project
   * @param projectId - Project ID to update
   * @param updates - Partial project updates
   * @returns The updated video project
   */
  updateProject(
    projectId: string,
    updates: Partial<VideoProject>
  ): Promise<VideoProject>;

  /**
   * Delete a video project
   * @param projectId - Project ID to delete
   * @returns True if deleted successfully
   */
  deleteProject(projectId: string): Promise<boolean>;

  // ============================================================================
  // Asset Management
  // ============================================================================

  /**
   * Add an asset to a video project
   * @param projectId - Project ID to add asset to
   * @param assetPath - Path to the asset file
   * @param assetType - Type of the asset
   * @returns The created asset
   */
  addAsset(
    projectId: string,
    assetPath: string,
    assetType: AssetType
  ): Promise<Asset>;

  /**
   * Remove an asset from a video project
   * @param projectId - Project ID to remove asset from
   * @param assetId - Asset ID to remove
   * @returns True if removed successfully
   */
  removeAsset(projectId: string, assetId: string): Promise<boolean>;

  // ============================================================================
  // Rendering
  // ============================================================================

  /**
   * Render a video composition
   * Shows a save dialog for output path selection.
   * @param options - Render options
   * @returns Output file path or null if cancelled
   */
  render(options: RenderOptions): Promise<string | null>;

  /**
   * Cancel the current render operation
   */
  cancelRender(): void;

  /**
   * Subscribe to render progress events
   * @param callback - Callback function for progress updates
   * @returns Cleanup function to unsubscribe
   * @requirements 15.4
   */
  onRenderProgress(callback: (progress: RenderProgressEvent) => void): () => void;

  // ============================================================================
  // Preview
  // ============================================================================

  /**
   * Start preview server for a video project
   * @param projectId - Project ID to preview
   * @returns Preview URL and port
   */
  startPreview(projectId: string): Promise<{ url: string; port: number }>;

  /**
   * Stop preview server for a video project
   * @param projectId - Project ID to stop preview for
   */
  stopPreview(projectId: string): Promise<void>;
}

// ============================================================================
// VideoAPI Implementation
// ============================================================================

/**
 * Create the VideoAPI implementation
 * @returns VideoAPI object with all methods implemented
 * @requirements 15.1, 15.2, 15.3, 15.4
 */
export function createVideoAPI(): VideoAPI {
  return {
    // ============================================================================
    // Project Management
    // ============================================================================

    createProject: (options: CreateProjectOptions): Promise<VideoProject> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CREATE_PROJECT, options);
    },

    listProjects: (workspaceId: string): Promise<VideoProject[]> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.LIST_PROJECTS, workspaceId);
    },

    getProject: (projectId: string): Promise<VideoProject | null> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.GET_PROJECT, projectId);
    },

    updateProject: (
      projectId: string,
      updates: Partial<VideoProject>
    ): Promise<VideoProject> => {
      return ipcRenderer.invoke(
        VIDEO_IPC_CHANNELS.UPDATE_PROJECT,
        projectId,
        updates
      );
    },

    deleteProject: (projectId: string): Promise<boolean> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.DELETE_PROJECT, projectId);
    },

    // ============================================================================
    // Asset Management
    // ============================================================================

    addAsset: (
      projectId: string,
      assetPath: string,
      assetType: AssetType
    ): Promise<Asset> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.ADD_ASSET, {
        projectId,
        assetPath,
        assetType,
      });
    },

    removeAsset: (projectId: string, assetId: string): Promise<boolean> => {
      return ipcRenderer.invoke(
        VIDEO_IPC_CHANNELS.REMOVE_ASSET,
        projectId,
        assetId
      );
    },

    // ============================================================================
    // Rendering
    // ============================================================================

    render: (options: RenderOptions): Promise<string | null> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.RENDER, options);
    },

    cancelRender: (): void => {
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CANCEL_RENDER);
    },

    /**
     * Subscribe to render progress events
     * @requirements 15.4
     */
    onRenderProgress: (
      callback: (progress: RenderProgressEvent) => void
    ): (() => void) => {
      const handler = (
        _event: Electron.IpcRendererEvent,
        progress: RenderProgressEvent
      ) => {
        callback(progress);
      };

      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);

      // Return cleanup function
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);
      };
    },

    // ============================================================================
    // Preview
    // ============================================================================

    startPreview: (
      projectId: string
    ): Promise<{ url: string; port: number }> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.START_PREVIEW, projectId);
    },

    stopPreview: (projectId: string): Promise<void> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.STOP_PREVIEW, projectId);
    },
  };
}

/**
 * Video API instance for use in preload script
 */
export const videoAPI = createVideoAPI();

// ============================================================================
// VideoServiceAPI Implementation
// ============================================================================

/**
 * Create the VideoServiceAPI implementation
 * @returns VideoServiceAPI object with all methods implemented
 * @requirements 5.1
 */
export function createVideoServiceAPI(): VideoServiceAPI {
  return {
    // ============================================================================
    // Service Lifecycle
    // ============================================================================

    startService: (): Promise<ServiceStatusResponse> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.START_VIDEO_SERVICE);
    },

    stopService: (): Promise<ServiceStatusResponse> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.STOP_VIDEO_SERVICE);
    },

    getStatus: (): Promise<ServiceStatusResponse> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.GET_SERVICE_STATUS);
    },

    // ============================================================================
    // Preview
    // ============================================================================

    startPreview: (
      projectId: string
    ): Promise<{ url: string; port: number }> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.START_PREVIEW, projectId);
    },

    stopPreview: (projectId: string): Promise<void> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.STOP_PREVIEW, projectId);
    },

    // ============================================================================
    // Rendering
    // ============================================================================

    render: (request: RenderOptions): Promise<string | null> => {
      return ipcRenderer.invoke(VIDEO_IPC_CHANNELS.RENDER, request);
    },

    cancelRender: (): void => {
      ipcRenderer.invoke(VIDEO_IPC_CHANNELS.CANCEL_RENDER);
    },

    // ============================================================================
    // Event Listeners
    // ============================================================================

    onStatusChanged: (
      callback: (status: ServiceStatusResponse) => void
    ): (() => void) => {
      const handler = (
        _event: Electron.IpcRendererEvent,
        status: ServiceStatusResponse
      ) => {
        callback(status);
      };

      ipcRenderer.on(VIDEO_IPC_CHANNELS.SERVICE_STATUS_CHANGED, handler);

      // Return cleanup function
      return () => {
        ipcRenderer.removeListener(
          VIDEO_IPC_CHANNELS.SERVICE_STATUS_CHANGED,
          handler
        );
      };
    },

    onRenderProgress: (
      callback: (progress: RenderProgressEvent) => void
    ): (() => void) => {
      const handler = (
        _event: Electron.IpcRendererEvent,
        progress: RenderProgressEvent
      ) => {
        callback(progress);
      };

      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);

      // Return cleanup function
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_PROGRESS, handler);
      };
    },

    onRenderCompleted: (
      callback: (event: RenderCompletedEvent) => void
    ): (() => void) => {
      const handler = (
        _event: Electron.IpcRendererEvent,
        completedEvent: RenderCompletedEvent
      ) => {
        callback(completedEvent);
      };

      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_COMPLETED, handler);

      // Return cleanup function
      return () => {
        ipcRenderer.removeListener(
          VIDEO_IPC_CHANNELS.RENDER_COMPLETED,
          handler
        );
      };
    },

    onRenderFailed: (
      callback: (event: RenderFailedEvent) => void
    ): (() => void) => {
      const handler = (
        _event: Electron.IpcRendererEvent,
        failedEvent: RenderFailedEvent
      ) => {
        callback(failedEvent);
      };

      ipcRenderer.on(VIDEO_IPC_CHANNELS.RENDER_FAILED, handler);

      // Return cleanup function
      return () => {
        ipcRenderer.removeListener(VIDEO_IPC_CHANNELS.RENDER_FAILED, handler);
      };
    },
  };
}

/**
 * Video Service API instance for use in preload script
 * @requirements 5.1
 */
export const videoServiceAPI = createVideoServiceAPI();
