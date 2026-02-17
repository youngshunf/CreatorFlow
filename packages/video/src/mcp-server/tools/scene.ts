/**
 * 场景管理工具
 *
 * MCP 工具：video_add_scene, video_update_scene, video_remove_scene, video_reorder_scenes
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import { createSuccessResponse } from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

export const AddSceneInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  compositionId: z.string().describe('内置组合 ID（如 "TitleAnimation"、"Slideshow" 等）'),
  durationInFrames: z.number().int().positive().describe('场景时长（帧数），例如 30fps 下 90 帧 = 3 秒'),
  props: z.record(z.string(), z.any()).optional().describe('传给组合组件的参数'),
  name: z.string().optional().describe('场景名称（可选，默认使用组合 ID）'),
  insertAt: z.number().int().nonnegative().optional().describe('插入位置索引（可选，默认追加到末尾）'),
});

export const UpdateSceneInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  sceneId: z.string().describe('场景ID'),
  props: z.record(z.string(), z.any()).optional().describe('更新的组合参数'),
  durationInFrames: z.number().int().positive().optional().describe('更新的时长（帧数）'),
  compositionId: z.string().optional().describe('更换组合类型'),
  name: z.string().optional().describe('更新场景名称'),
});

export const RemoveSceneInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  sceneId: z.string().describe('场景ID'),
});

export const ReorderScenesInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  sceneIds: z.array(z.string()).describe('新的场景 ID 顺序（必须包含所有场景 ID）'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

async function handleAddScene(
  input: z.infer<typeof AddSceneInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const scene = await store.addScene(input.projectId, {
      name: input.name ?? input.compositionId,
      compositionId: input.compositionId,
      durationInFrames: input.durationInFrames,
      props: input.props,
      insertAt: input.insertAt,
    });

    // 返回场景和项目最新状态
    const project = await store.getProject(input.projectId);
    return JSON.stringify(createSuccessResponse({
      scene,
      projectState: project ? {
        scenes: project.scenes,
        transitions: project.transitions,
      } : undefined,
    }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

async function handleUpdateScene(
  input: z.infer<typeof UpdateSceneInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const scene = await store.updateScene(input.projectId, input.sceneId, {
      props: input.props,
      durationInFrames: input.durationInFrames,
      compositionId: input.compositionId,
      name: input.name,
    });
    return JSON.stringify(createSuccessResponse(scene));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

async function handleRemoveScene(
  input: z.infer<typeof RemoveSceneInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const success = await store.removeScene(input.projectId, input.sceneId);

    const project = await store.getProject(input.projectId);
    return JSON.stringify(createSuccessResponse({
      removed: success,
      projectState: project ? {
        scenes: project.scenes,
        transitions: project.transitions,
      } : undefined,
    }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

async function handleReorderScenes(
  input: z.infer<typeof ReorderScenesInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const project = await store.reorderScenes(input.projectId, input.sceneIds);
    return JSON.stringify(createSuccessResponse({
      scenes: project.scenes,
      transitions: project.transitions,
    }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

export function registerSceneTools(mcp: FastMCP): void {
  mcp.addTool({
    name: 'video_add_scene',
    description: `向视频项目添加一个场景片段。

场景是视频的基本构建块，每个场景引用一个内置组合（如 TitleAnimation、Slideshow 等），
并通过 props 配置其内容和外观。

使用 video_list_available_compositions 查看所有可用的组合及其参数说明。`,
    parameters: AddSceneInputSchema,
    execute: handleAddScene,
  });

  mcp.addTool({
    name: 'video_update_scene',
    description: '更新场景片段的参数、时长或组合类型。只更新提供的字段，其他字段保持不变。',
    parameters: UpdateSceneInputSchema,
    execute: handleUpdateScene,
  });

  mcp.addTool({
    name: 'video_remove_scene',
    description: '从视频项目移除场景片段。关联的过渡效果会自动调整。',
    parameters: RemoveSceneInputSchema,
    execute: handleRemoveScene,
  });

  mcp.addTool({
    name: 'video_reorder_scenes',
    description: '重新排列场景顺序。需要传入所有场景 ID 的新顺序。',
    parameters: ReorderScenesInputSchema,
    execute: handleReorderScenes,
  });
}

export {
  handleAddScene,
  handleUpdateScene,
  handleRemoveScene,
  handleReorderScenes,
};
