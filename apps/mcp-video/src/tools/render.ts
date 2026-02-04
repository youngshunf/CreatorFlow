/**
 * 渲染工具
 *
 * MCP 工具：video_render
 *
 * @requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import { RenderEngine } from '../services/render-engine';
import {
  type RenderResult,
  type OutputFormat,
  type QualityPreset,
  createSuccessResponse,
} from '../types';
import { toErrorResponse, createProjectNotFoundError } from '../types/errors';
import { getOutputPath } from '../utils/paths';
import { join } from 'path';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_render 输入 Schema
 */
export const RenderInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  compositionId: z.string().describe('组合ID'),
  outputFormat: z
    .enum(['mp4', 'webm', 'gif'])
    .optional()
    .default('mp4')
    .describe('输出格式'),
  quality: z
    .enum(['draft', 'standard', 'high'])
    .optional()
    .default('standard')
    .describe('质量预设'),
  outputPath: z.string().optional().describe('输出文件路径（可选，默认自动生成）'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 渲染视频到文件
 *
 * @requirements 5.1 - 启动渲染指定的组合
 * @requirements 5.2 - 支持 MP4、WebM、GIF 格式
 * @requirements 5.3 - 支持质量预设
 * @requirements 5.4 - 提供进度更新
 * @requirements 5.5 - 返回输出文件路径和渲染统计
 * @requirements 5.6 - 渲染失败时返回错误信息
 */
async function handleRender(
  input: z.infer<typeof RenderInputSchema>
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

    // 确定输出路径
    let outputPath = input.outputPath;
    if (!outputPath) {
      // 生成带时间戳的文件名
      const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);
      const fileName = `${input.compositionId}_${timestamp}.${input.outputFormat}`;
      outputPath = join(getOutputPath(input.workspacePath, project.name), fileName);
    }

    // 创建渲染引擎
    const renderEngine = RenderEngine.create();

    // 执行渲染
    // 注意：由于 MCP 工具是同步返回的，我们无法实时推送进度
    // 进度信息会在日志中输出
    const result = await renderEngine.render(
      {
        projectPath,
        compositionId: input.compositionId,
        outputFormat: input.outputFormat as OutputFormat,
        quality: input.quality as QualityPreset,
        outputPath,
      },
      (progress) => {
        // 进度回调 - 在日志中输出
        console.log(
          `[video_render] Progress: ${progress.status} - ${progress.progress}%`
        );
      }
    );

    return JSON.stringify(createSuccessResponse(result));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册渲染工具到 FastMCP 服务器
 */
export function registerRenderTools(mcp: FastMCP): void {
  // video_render
  mcp.addTool({
    name: 'video_render',
    description: `渲染视频到文件。

支持的输出格式：
- mp4: H.264 编码，最广泛支持
- webm: VP8 编码，适合网页
- gif: 动图格式，适合短循环动画

质量预设：
- draft: 快速渲染，较低质量，适合预览
- standard: 平衡质量和速度，适合大多数场景
- high: 最高质量，渲染较慢，适合最终输出

渲染完成后返回输出文件路径和统计信息（耗时、文件大小）。`,
    parameters: RenderInputSchema,
    execute: handleRender,
  });
}

// ============================================================================
// 导出
// ============================================================================

export { handleRender };
