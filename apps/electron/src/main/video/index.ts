/**
 * Video Service Module
 *
 * 导出视频 IPC handlers、渲染工作器及相关类型。
 * 已移除 MCP Server 子进程管理和预览服务器依赖。
 */

// ============================================================================
// IPC Handlers
// ============================================================================

export {
  registerVideoIpcHandlers,
  unregisterVideoIpcHandlers,
  cleanupVideoServices,
  VIDEO_IPC_CHANNELS,
  type CreateProjectRequest,
  type AddSceneRequest,
  type UpdateSceneRequest,
  type AddAssetRequest,
  type RenderRequest,
  type CreateFromTemplateRequest,
  type CompositionInfo,
} from './ipc-handlers';

// ============================================================================
// Types and Interfaces - Render Worker
// ============================================================================

export type {
  RenderWorkerOptions,
  RenderError,
  IRenderWorker,
  RenderScriptQualityPreset,
} from './render-worker';

// ============================================================================
// Enums and Classes - Render Worker
// ============================================================================

export {
  RenderErrorType,
  RenderWorker,
} from './render-worker';
