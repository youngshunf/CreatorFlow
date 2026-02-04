/**
 * 组合管理工具
 *
 * MCP 工具：video_add_composition, video_update_composition, video_remove_composition
 *
 * @requirements 7.1, 7.2, 7.3, 7.4, 7.5
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import {
  type Composition,
  createSuccessResponse,
} from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_add_composition 输入 Schema
 */
export const AddCompositionInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  name: z.string().describe('组合名称'),
  code: z.string().describe('组合代码（React 组件代码或文件路径）'),
  props: z.record(z.string(), z.any()).optional().describe('组合属性'),
});

/**
 * video_update_composition 输入 Schema
 */
export const UpdateCompositionInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  compositionId: z.string().describe('组合ID'),
  name: z.string().optional().describe('新的组合名称'),
  code: z.string().optional().describe('新的组合代码'),
  props: z.record(z.string(), z.any()).optional().describe('新的组合属性'),
});

/**
 * video_remove_composition 输入 Schema
 */
export const RemoveCompositionInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  compositionId: z.string().describe('组合ID'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 添加组合到视频项目
 *
 * @requirements 7.1 - 添加新组合到项目
 * @requirements 7.2 - 验证组合代码引用有效组件
 */
async function handleAddComposition(
  input: z.infer<typeof AddCompositionInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const composition = await store.addComposition(input.projectId, {
      name: input.name,
      code: input.code,
      props: input.props,
    });
    return JSON.stringify(createSuccessResponse(composition));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 更新视频组合
 *
 * @requirements 7.3 - 更新组合的属性和代码
 * @requirements 7.5 - 组合不存在时返回错误
 */
async function handleUpdateComposition(
  input: z.infer<typeof UpdateCompositionInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const composition = await store.updateComposition(
      input.projectId,
      input.compositionId,
      {
        name: input.name,
        code: input.code,
        props: input.props,
      }
    );
    return JSON.stringify(createSuccessResponse(composition));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 从视频项目移除组合
 *
 * @requirements 7.4 - 从项目移除组合
 * @requirements 7.5 - 组合不存在时返回错误
 */
async function handleRemoveComposition(
  input: z.infer<typeof RemoveCompositionInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const success = await store.removeComposition(
      input.projectId,
      input.compositionId
    );
    return JSON.stringify(createSuccessResponse({ removed: success }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册组合管理工具到 FastMCP 服务器
 */
export function registerCompositionTools(mcp: FastMCP): void {
  // video_add_composition
  mcp.addTool({
    name: 'video_add_composition',
    description: `添加组合到视频项目。组合是视频的基本构建块，定义了视频内容和动画。

组合代码可以是：
- React 组件代码字符串
- 指向组件文件的路径

组合属性（props）用于配置组合的行为和外观。`,
    parameters: AddCompositionInputSchema,
    execute: handleAddComposition,
  });

  // video_update_composition
  mcp.addTool({
    name: 'video_update_composition',
    description: '更新视频组合的名称、代码或属性。只更新提供的字段，其他字段保持不变。',
    parameters: UpdateCompositionInputSchema,
    execute: handleUpdateComposition,
  });

  // video_remove_composition
  mcp.addTool({
    name: 'video_remove_composition',
    description: '从视频项目移除组合。此操作不可撤销。',
    parameters: RemoveCompositionInputSchema,
    execute: handleRemoveComposition,
  });
}

// ============================================================================
// 导出
// ============================================================================

export {
  handleAddComposition,
  handleUpdateComposition,
  handleRemoveComposition,
};
