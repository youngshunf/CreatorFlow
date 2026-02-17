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

// 导出场景管理工具
export {
  registerSceneTools,
  handleAddScene,
  handleUpdateScene,
  handleRemoveScene,
  handleReorderScenes,
  AddSceneInputSchema,
  UpdateSceneInputSchema,
  RemoveSceneInputSchema,
  ReorderScenesInputSchema,
} from './scene';

// 导出过渡效果工具
export {
  registerTransitionTools,
  handleSetTransitions,
  SetTransitionsInputSchema,
} from './transition';

// 导出组合列表工具
export {
  registerCompositionListTools,
  handleListAvailableCompositions,
  ListAvailableCompositionsInputSchema,
} from './composition-list';

// 导出渲染工具
export {
  registerRenderTools,
  handleRender,
  RenderInputSchema,
} from './render';

// 导出模板工具
export {
  registerTemplateTools,
  handleListTemplates,
  handleGetTemplate,
  ListTemplatesInputSchema,
  GetTemplateInputSchema,
} from './template';

// 导出素材发现工具
export {
  registerAssetDiscoveryTools,
  handleListAvailableAssets,
  ListAvailableAssetsInputSchema,
} from './asset-discovery';

// 导出代码验证工具
export {
  registerCodeValidationTools,
  handleValidateComposition,
  ValidateCompositionInputSchema,
} from './code-validation';

// 导出渲染状态工具
export {
  registerRenderStatusTools,
  handleGetRenderStatus,
  GetRenderStatusInputSchema,
  updateRenderStatus,
  clearRenderStatus,
} from './render-status';

// 导入注册函数
import { registerProjectTools } from './project';
import { registerAssetTools } from './asset';
import { registerSceneTools } from './scene';
import { registerTransitionTools } from './transition';
import { registerCompositionListTools } from './composition-list';
import { registerRenderTools } from './render';
import { registerTemplateTools } from './template';
import { registerAssetDiscoveryTools } from './asset-discovery';
import { registerCodeValidationTools } from './code-validation';
import { registerRenderStatusTools } from './render-status';

/**
 * 注册所有 MCP 工具到 FastMCP 服务器
 *
 * @param mcp FastMCP 服务器实例
 */
export function registerAllTools(mcp: FastMCP): void {
  console.error('[MCP Video Server] Registering tools...');

  // 注册项目管理工具
  registerProjectTools(mcp);
  console.error('[MCP Video Server] - Project tools registered');

  // 注册素材管理工具
  registerAssetTools(mcp);
  console.error('[MCP Video Server] - Asset tools registered');

  // 注册场景管理工具
  registerSceneTools(mcp);
  console.error('[MCP Video Server] - Scene tools registered');

  // 注册过渡效果工具
  registerTransitionTools(mcp);
  console.error('[MCP Video Server] - Transition tools registered');

  // 注册组合列表工具
  registerCompositionListTools(mcp);
  console.error('[MCP Video Server] - Composition list tools registered');

  // 注册渲染工具
  registerRenderTools(mcp);
  console.error('[MCP Video Server] - Render tools registered');

  // 注册模板工具
  registerTemplateTools(mcp);
  console.error('[MCP Video Server] - Template tools registered');

  // 注册素材发现工具
  registerAssetDiscoveryTools(mcp);
  console.error('[MCP Video Server] - Asset discovery tools registered');

  // 注册代码验证工具
  registerCodeValidationTools(mcp);
  console.error('[MCP Video Server] - Code validation tools registered');

  // 注册渲染状态工具
  registerRenderStatusTools(mcp);
  console.error('[MCP Video Server] - Render status tools registered');

  console.error('[MCP Video Server] All tools registered successfully');
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

  // 场景管理
  { name: 'video_add_scene', category: 'scene', description: '向项目添加场景片段' },
  { name: 'video_update_scene', category: 'scene', description: '更新场景参数' },
  { name: 'video_remove_scene', category: 'scene', description: '移除场景片段' },
  { name: 'video_reorder_scenes', category: 'scene', description: '重排场景顺序' },

  // 过渡效果
  { name: 'video_set_transitions', category: 'transition', description: '设置场景间过渡效果' },

  // 组合列表
  { name: 'video_list_available_compositions', category: 'composition', description: '列出所有可用的内置组合及其参数' },

  // 渲染
  { name: 'video_render', category: 'render', description: '渲染视频到文件' },

  // 模板
  { name: 'video_list_templates', category: 'template', description: '列出所有可用的视频模板' },
  { name: 'video_get_template', category: 'template', description: '获取视频模板详情' },

  // 素材发现
  { name: 'video_list_available_assets', category: 'asset-discovery', description: '列出工作区中可用的素材文件' },

  // 代码验证
  { name: 'video_validate_composition', category: 'code-validation', description: '验证 Remotion 组合代码的语法和正确性' },

  // 渲染状态
  { name: 'video_get_render_status', category: 'render-status', description: '获取视频渲染的当前状态和进度' },
] as const;

/**
 * 工具名称类型
 */
export type ToolName = (typeof TOOL_LIST)[number]['name'];

/**
 * 工具分类类型
 */
export type ToolCategory = (typeof TOOL_LIST)[number]['category'];
