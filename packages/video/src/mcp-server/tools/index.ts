/**
 * MCP 工具导出
 *
 * 导出所有 MCP 工具供服务器注册
 */

import type { FastMCP } from 'fastmcp';

// 导出项目管理工具
export {
  registerProjectTools,
  handleCreateProject,
  handleListProjects,
  handleGetProject,
  handleUpdateProject,
  handleDeleteProject,
  CreateProjectInputSchema,
  ListProjectsInputSchema,
  GetProjectInputSchema,
  UpdateProjectInputSchema,
  DeleteProjectInputSchema,
} from './project';

// 导出素材管理工具
export {
  registerAssetTools,
  handleAddAsset,
  handleRemoveAsset,
  handleListAssets,
  AddAssetInputSchema,
  RemoveAssetInputSchema,
  ListAssetsInputSchema,
} from './asset';

// 导出组合管理工具
export {
  registerCompositionTools,
  handleAddComposition,
  handleUpdateComposition,
  handleRemoveComposition,
  AddCompositionInputSchema,
  UpdateCompositionInputSchema,
  RemoveCompositionInputSchema,
} from './composition';

// 导出渲染工具
export {
  registerRenderTools,
  handleRender,
  RenderInputSchema,
} from './render';

// 导出预览工具
export {
  registerPreviewTools,
  handlePreviewStart,
  handlePreviewStop,
  PreviewStartInputSchema,
  PreviewStopInputSchema,
} from './preview';

// 导出模板工具
export {
  registerTemplateTools,
  handleListTemplates,
  handleGetTemplate,
  ListTemplatesInputSchema,
  GetTemplateInputSchema,
} from './template';

// 导入注册函数
import { registerProjectTools } from './project';
import { registerAssetTools } from './asset';
import { registerCompositionTools } from './composition';
import { registerRenderTools } from './render';
import { registerPreviewTools } from './preview';
import { registerTemplateTools } from './template';

/**
 * 注册所有 MCP 工具到 FastMCP 服务器
 *
 * @param mcp FastMCP 服务器实例
 */
export function registerAllTools(mcp: FastMCP): void {
  console.log('[MCP Video Server] Registering tools...');

  // 注册项目管理工具
  registerProjectTools(mcp);
  console.log('[MCP Video Server] - Project tools registered');

  // 注册素材管理工具
  registerAssetTools(mcp);
  console.log('[MCP Video Server] - Asset tools registered');

  // 注册组合管理工具
  registerCompositionTools(mcp);
  console.log('[MCP Video Server] - Composition tools registered');

  // 注册渲染工具
  registerRenderTools(mcp);
  console.log('[MCP Video Server] - Render tools registered');

  // 注册预览工具
  registerPreviewTools(mcp);
  console.log('[MCP Video Server] - Preview tools registered');

  // 注册模板工具
  registerTemplateTools(mcp);
  console.log('[MCP Video Server] - Template tools registered');

  console.log('[MCP Video Server] All tools registered successfully');
}

/**
 * 工具列表（用于文档和调试）
 */
export const TOOL_LIST = [
  // 项目管理
  { name: 'video_create_project', category: 'project', description: '创建新的视频项目' },
  { name: 'video_list_projects', category: 'project', description: '列出工作区中的所有视频项目' },
  { name: 'video_get_project', category: 'project', description: '获取视频项目详情' },
  { name: 'video_update_project', category: 'project', description: '更新视频项目' },
  { name: 'video_delete_project', category: 'project', description: '删除视频项目' },

  // 素材管理
  { name: 'video_add_asset', category: 'asset', description: '添加素材到视频项目' },
  { name: 'video_remove_asset', category: 'asset', description: '从视频项目移除素材' },
  { name: 'video_list_assets', category: 'asset', description: '列出项目中的所有素材' },

  // 组合管理
  { name: 'video_add_composition', category: 'composition', description: '添加组合到视频项目' },
  { name: 'video_update_composition', category: 'composition', description: '更新视频组合' },
  { name: 'video_remove_composition', category: 'composition', description: '从视频项目移除组合' },

  // 渲染
  { name: 'video_render', category: 'render', description: '渲染视频到文件' },

  // 预览
  { name: 'video_preview_start', category: 'preview', description: '启动视频预览服务器' },
  { name: 'video_preview_stop', category: 'preview', description: '停止视频预览服务器' },

  // 模板
  { name: 'video_list_templates', category: 'template', description: '列出所有可用的视频模板' },
  { name: 'video_get_template', category: 'template', description: '获取视频模板详情' },
] as const;

/**
 * 工具名称类型
 */
export type ToolName = (typeof TOOL_LIST)[number]['name'];

/**
 * 工具分类类型
 */
export type ToolCategory = (typeof TOOL_LIST)[number]['category'];
