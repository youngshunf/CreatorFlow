/**
 * 项目管理工具
 *
 * MCP 工具：video_create_project, video_list_projects, video_get_project,
 * video_update_project, video_delete_project
 *
 * @requirements 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 2.1, 2.2, 2.3, 2.4, 3.1, 3.2, 3.3, 3.4
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import {
  type VideoProject,
  type ProjectSummary,
  createSuccessResponse,
} from '../types';
import { toErrorResponse, MCPError } from '../types/errors';
import { getTemplateById } from '../../templates';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_create_project 输入 Schema
 */
export const CreateProjectInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  name: z.string().describe('项目名称'),
  template: z.string().optional().describe('模板ID'),
  width: z.number().optional().default(1920).describe('视频宽度（像素）'),
  height: z.number().optional().default(1080).describe('视频高度（像素）'),
  fps: z.number().optional().default(30).describe('帧率'),
  durationInSeconds: z.number().describe('视频时长（秒）'),
  description: z.string().optional().describe('项目描述'),
});

/**
 * video_list_projects 输入 Schema
 */
export const ListProjectsInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
});

/**
 * video_get_project 输入 Schema
 */
export const GetProjectInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
});

/**
 * video_update_project 输入 Schema
 */
export const UpdateProjectInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  name: z.string().optional().describe('新的项目名称'),
  description: z.string().optional().describe('新的项目描述'),
  config: z
    .object({
      width: z.number().optional().describe('视频宽度'),
      height: z.number().optional().describe('视频高度'),
      fps: z.number().optional().describe('帧率'),
      durationInFrames: z.number().optional().describe('总帧数'),
    })
    .optional()
    .describe('视频配置'),
});

/**
 * video_delete_project 输入 Schema
 */
export const DeleteProjectInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 创建视频项目
 *
 * @requirements 1.1 - 创建新项目并返回项目详情
 * @requirements 1.2 - 支持模板参数
 * @requirements 1.3 - 使用默认值
 * @requirements 1.4 - 验证项目名称
 * @requirements 1.5 - 持久化到文件系统
 * @requirements 1.6 - 需要 workspacePath 参数
 * @requirements 1.7 - 创建用户友好的目录结构
 */
async function handleCreateProject(
  input: z.infer<typeof CreateProjectInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);

    // 如果指定了模板，获取模板配置
    let templateConfig: {
      width?: number;
      height?: number;
      fps?: number;
      durationInFrames?: number;
      compositions?: Array<{ name: string; code: string; props?: Record<string, any> }>;
    } = {};

    if (input.template) {
      const template = getTemplateById(input.template);
      if (template) {
        templateConfig = {
          width: template.defaultConfig.width,
          height: template.defaultConfig.height,
          fps: template.defaultConfig.fps,
          durationInFrames: template.defaultConfig.durationInFrames,
          compositions: [
            {
              name: template.name,
              code: template.compositionCode,
              props: template.defaultProps,
            },
          ],
        };
      }
    }

    // 创建项目
    const project = await store.createProject({
      name: input.name,
      template: input.template,
      width: input.width ?? templateConfig.width ?? 1920,
      height: input.height ?? templateConfig.height ?? 1080,
      fps: input.fps ?? templateConfig.fps ?? 30,
      durationInSeconds: input.durationInSeconds,
      description: input.description,
    });

    // 如果有模板组合，添加到项目
    if (templateConfig.compositions) {
      for (const comp of templateConfig.compositions) {
        await store.addComposition(project.id, {
          name: comp.name,
          code: comp.code,
          props: comp.props,
        });
      }
      // 重新获取项目以包含组合
      const updatedProject = await store.getProject(project.id);
      if (updatedProject) {
        return JSON.stringify(createSuccessResponse(updatedProject));
      }
    }

    return JSON.stringify(createSuccessResponse(project));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 列出工作区中的所有视频项目
 *
 * @requirements 2.1 - 返回项目列表
 * @requirements 2.4 - 按更新时间降序排列
 */
async function handleListProjects(
  input: z.infer<typeof ListProjectsInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const projects = await store.listProjects();
    return JSON.stringify(createSuccessResponse(projects));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 获取视频项目详情
 *
 * @requirements 2.2 - 返回完整项目详情
 * @requirements 2.3 - 项目不存在时返回错误
 */
async function handleGetProject(
  input: z.infer<typeof GetProjectInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const project = await store.getProject(input.projectId);

    if (!project) {
      throw new MCPError(
        'PROJECT_NOT_FOUND',
        `项目不存在: ${input.projectId}`,
        { field: 'projectId', received: input.projectId }
      );
    }

    return JSON.stringify(createSuccessResponse(project));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 更新视频项目
 *
 * @requirements 3.1 - 更新指定字段
 * @requirements 3.2 - 自动更新 updatedAt
 * @requirements 3.4 - 项目不存在时返回错误
 */
async function handleUpdateProject(
  input: z.infer<typeof UpdateProjectInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const project = await store.updateProject(input.projectId, {
      name: input.name,
      description: input.description,
      config: input.config,
    });
    return JSON.stringify(createSuccessResponse(project));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 删除视频项目
 *
 * @requirements 3.3 - 删除项目和所有关联文件
 * @requirements 3.4 - 项目不存在时返回错误
 */
async function handleDeleteProject(
  input: z.infer<typeof DeleteProjectInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const success = await store.deleteProject(input.projectId);
    return JSON.stringify(createSuccessResponse({ deleted: success }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册项目管理工具到 FastMCP 服务器
 */
export function registerProjectTools(mcp: FastMCP): void {
  // video_create_project
  mcp.addTool({
    name: 'video_create_project',
    description: '创建新的视频项目。支持指定模板、分辨率、帧率和时长。项目将保存在工作区的"视频创作"目录下。',
    parameters: CreateProjectInputSchema,
    execute: handleCreateProject,
  });

  // video_list_projects
  mcp.addTool({
    name: 'video_list_projects',
    description: '列出工作区中的所有视频项目。返回项目摘要列表，按更新时间降序排列。',
    parameters: ListProjectsInputSchema,
    execute: handleListProjects,
  });

  // video_get_project
  mcp.addTool({
    name: 'video_get_project',
    description: '获取视频项目的完整详情，包括配置、组合和素材列表。',
    parameters: GetProjectInputSchema,
    execute: handleGetProject,
  });

  // video_update_project
  mcp.addTool({
    name: 'video_update_project',
    description: '更新视频项目的名称、描述或配置。只更新提供的字段，其他字段保持不变。',
    parameters: UpdateProjectInputSchema,
    execute: handleUpdateProject,
  });

  // video_delete_project
  mcp.addTool({
    name: 'video_delete_project',
    description: '删除视频项目及其所有关联文件（素材、组合、输出等）。此操作不可撤销。',
    parameters: DeleteProjectInputSchema,
    execute: handleDeleteProject,
  });
}

// ============================================================================
// 导出
// ============================================================================

export {
  handleCreateProject,
  handleListProjects,
  handleGetProject,
  handleUpdateProject,
  handleDeleteProject,
};
