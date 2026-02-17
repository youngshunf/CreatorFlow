/**
 * 渲染工具
 *
 * MCP 工具：video_render
 * 使用统一入口 SceneComposer 渲染，从项目读取 scenes + transitions 作为 inputProps
 *
 * @requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6
 */

import { z } from "zod";
import type { FastMCP } from "fastmcp";
import { ProjectStore } from "../services/project-store";
import { RenderEngine } from "../services/render-engine";
import {
  type RenderResult,
  type OutputFormat,
  type QualityPreset,
  createSuccessResponse,
} from "../types";
import { toErrorResponse, createProjectNotFoundError } from "../types/errors";
import { getOutputPath } from "../utils/paths";
import { join } from "path";

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_render 输入 Schema
 */
export const RenderInputSchema = z.object({
  workspacePath: z.string().describe("工作区根路径"),
  projectId: z.string().describe("项目ID"),
  outputFormat: z
    .enum(["mp4", "webm", "gif"])
    .optional()
    .default("mp4")
    .describe("输出格式"),
  quality: z
    .enum(["draft", "standard", "high"])
    .optional()
    .default("standard")
    .describe("质量预设"),
  outputPath: z
    .string()
    .optional()
    .describe("输出文件路径（可选，默认自动生成）"),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 渲染视频到文件
 *
 * 从项目读取 scenes + transitions，使用 SceneComposer 统一渲染。
 * 不再需要 compositionId 参数（固定为 SceneComposer）。
 */
async function handleRender(
  input: z.infer<typeof RenderInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);

    // 获取项目
    const project = await store.getProject(input.projectId);
    if (!project) {
      throw createProjectNotFoundError(input.projectId);
    }

    // 确定输出路径
    let outputPath = input.outputPath;
    if (!outputPath) {
      const timestamp = new Date()
        .toISOString()
        .replace(/[:.]/g, "-")
        .slice(0, 19);
      const fileName = `SceneComposer_${timestamp}.${input.outputFormat}`;
      outputPath = join(
        getOutputPath(input.workspacePath, project.name),
        fileName,
      );
    }

    // 创建渲染引擎
    const renderEngine = RenderEngine.create();

    // 从项目读取 scenes + transitions 作为 inputProps
    const inputProps = {
      scenes: project.scenes,
      transitions: project.transitions,
    };

    // 执行渲染
    const result = await renderEngine.render(
      {
        outputFormat: input.outputFormat as OutputFormat,
        quality: input.quality as QualityPreset,
        inputProps,
        outputPath,
      },
      (progress) => {
        console.error(
          `[video_render] Progress: ${progress.status} - ${progress.progress}%`,
        );
      },
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
    name: "video_render",
    description: `渲染视频到文件。自动从项目读取场景和过渡效果，使用 SceneComposer 统一渲染。

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
