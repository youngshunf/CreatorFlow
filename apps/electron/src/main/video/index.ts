/**
 * Video Service Module
 *
 * 导出视频服务管理器、预览服务器及相关类型。
 * 提供 MCP Video Server 子进程管理、预览服务器管理、服务状态跟踪和错误处理功能。
 *
 * @requirements 2.1 - 使用内置 Bun 运行时启动 MCP Video Server
 * @requirements 3.1, 3.2 - 使用内置 Bun 运行时启动预览服务器
 * @requirements 5.1, 5.2 - IPC 通道和处理器
 *
 * @example
 * ```typescript
 * import {
 *   createVideoServiceManager,
 *   createPreviewServer,
 *   registerVideoIpcHandlers,
 *   cleanupVideoServices,
 *   VIDEO_IPC_CHANNELS,
 *   VideoServiceError,
 *   VideoServiceManagerError,
 *   PreviewServerError,
 *   PreviewServerException,
 *   type VideoServiceConfig,
 *   type ServiceStatus,
 *   type VideoServiceManager,
 *   type PreviewServerConfig,
 *   type PreviewServerResult,
 *   type IPreviewServer,
 * } from './video';
 *
 * // 注册 IPC 处理器
 * registerVideoIpcHandlers();
 *
 * // 创建视频服务管理器
 * const manager = createVideoServiceManager(bunResolver, {
 *   mcpServerEntry: '/path/to/mcp-video/src/index.ts',
 * });
 * await manager.startMcpServer();
 *
 * // 创建预览服务器
 * const previewServer = createPreviewServer(bunResolver);
 * const result = await previewServer.start('/path/to/project');
 * console.log(`Preview at ${result.url} (PID: ${result.pid})`);
 *
 * // 清理服务
 * await cleanupVideoServices();
 * ```
 */

// ============================================================================
// IPC Handlers
// ============================================================================

export {
  /** 注册视频 IPC 处理器 */
  registerVideoIpcHandlers,
  /** 注销视频 IPC 处理器 */
  unregisterVideoIpcHandlers,
  /** 清理所有视频服务 */
  cleanupVideoServices,
  /** 获取视频服务管理器实例 */
  getVideoServiceManagerInstance,
  /** 视频 IPC 通道定义 */
  VIDEO_IPC_CHANNELS,
  /** 服务状态响应类型 */
  type ServiceStatusResponse,
} from './ipc-handlers';

// ============================================================================
// Types and Interfaces - Service Manager
// ============================================================================

export type {
  /** 视频服务配置 */
  VideoServiceConfig,
  /** 服务状态 */
  ServiceStatus,
  /** 视频服务管理器接口 */
  VideoServiceManager,
} from './service-manager';

// ============================================================================
// Types and Interfaces - Preview Server
// ============================================================================

export type {
  /** 预览服务器配置 */
  PreviewServerConfig,
  /** 预览服务器结果 */
  PreviewServerResult,
  /** 预览服务器接口 */
  IPreviewServer,
} from './preview-server';

// ============================================================================
// Enums and Classes - Service Manager
// ============================================================================

export {
  /** 视频服务错误类型枚举 */
  VideoServiceError,
  /** 视频服务错误类 */
  VideoServiceManagerError,
} from './service-manager';

// ============================================================================
// Enums and Classes - Preview Server
// ============================================================================

export {
  /** 预览服务器错误类型枚举 */
  PreviewServerError,
  /** 预览服务器错误类 */
  PreviewServerException,
  /** 预览服务器类 */
  PreviewServer,
} from './preview-server';

// ============================================================================
// Factory Functions
// ============================================================================

export {
  /** 创建视频服务管理器 */
  createVideoServiceManager,
} from './service-manager';

export {
  /** 创建预览服务器 */
  createPreviewServer,
} from './preview-server';

// ============================================================================
// Types and Interfaces - Render Worker
// ============================================================================

export type {
  /** 渲染工作器选项 */
  RenderWorkerOptions,
  /** 渲染错误 */
  RenderError,
  /** 渲染工作器接口 */
  IRenderWorker,
  /** 渲染脚本质量预设 */
  RenderScriptQualityPreset,
} from './render-worker';

// ============================================================================
// Enums and Classes - Render Worker
// ============================================================================

export {
  /** 渲染错误类型枚举 */
  RenderErrorType,
  /** 渲染工作器类 */
  RenderWorker,
} from './render-worker';
