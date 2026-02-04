/**
 * 预览工具
 *
 * MCP 工具：video_preview_start, video_preview_stop
 *
 * @requirements 6.1, 6.3, 6.4, 6.5
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import { previewServerManager } from '../services/preview-server';
import {
  type PreviewInstance,
  createSuccessResponse,
} from '../types';
import { toErrorResponse, createProjectNotFoundError } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_preview_start 输入 Schema
 */
export const PreviewStartInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  compositionId: z.string().optional().describe('组合ID（可选，指定默认显示的组合）'),
});

/**
 * video_preview_stop 输入 Schema
 */
export const PreviewStopInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 启动视频预览服务器
 *
 * @requirements 6.1 - 启动预览服务器并返回 URL
 * @requirements 6.4 - 启动失败时返回错误
 * @requirements 6.5 - 防止同一项目重复启动服务器
 */
async function handlePreviewStart(
  input: z.infer<typeof PreviewStartInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);

    // 获取项目
    const project = await store.getProject(input.projectId);
    if (!project) {
      throw createProjectNotFoundError(input.projectId);
    }

    // 获取项目路径
    const projectPath = store.getProjectPath(project.name);

    // 启动预览服务器
    const instance = await previewServerManager.start({
      projectPath,
      projectId: input.projectId,
      compositionId: input.compositionId,
    });

    return JSON.stringify(createSuccessResponse(instance));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 停止视频预览服务器
 *
 * @requirements 6.3 - 优雅关闭服务器并释放端口
 */
async function handlePreviewStop(
  input: z.infer<typeof PreviewStopInputSchema>
): Promise<string> {
  try {
    // 停止预览服务器
    const success = await previewServerManager.stop(input.projectId);

    return JSON.stringify(
      createSuccessResponse({
        stopped: success,
        projectId: input.projectId,
      })
    );
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册预览工具到 FastMCP 服务器
 */
export function registerPreviewTools(mcp: FastMCP): void {
  // video_preview_start
  mcp.addTool({
    name: 'video_preview_start',
    description: `启动视频预览服务器。

预览服务器基于 Remotion Studio，提供：
- 实时预览视频组合
- 热重载支持（修改代码后自动刷新）
- 时间轴控制
- 组合切换

如果项目已有活跃的预览服务器，将返回现有服务器的 URL。
返回预览 URL，可在浏览器中打开查看。`,
    parameters: PreviewStartInputSchema,
    execute: handlePreviewStart,
  });

  // video_preview_stop
  mcp.addTool({
    name: 'video_preview_stop',
    description: '停止视频预览服务器，释放端口资源。',
    parameters: PreviewStopInputSchema,
    execute: handlePreviewStop,
  });
}

// ============================================================================
// 导出
// ============================================================================

export { handlePreviewStart, handlePreviewStop };
