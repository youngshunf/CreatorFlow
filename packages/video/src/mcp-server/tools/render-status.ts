/**
 * 渲染状态工具
 *
 * MCP 工具：video_get_render_status
 *
 * 获取视频渲染的当前状态和进度
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { existsSync, readFileSync } from 'fs';
import { join } from 'path';
import { createSuccessResponse } from '../types';
import { toErrorResponse } from '../types/errors';
import { getProjectPath, getOutputPath } from '../utils/paths';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_get_render_status 输入 Schema
 */
export const GetRenderStatusInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectName: z.string().describe('项目名称'),
  renderId: z.string().optional().describe('渲染任务 ID（可选）'),
});

/**
 * 渲染状态类型
 */
type RenderStatus = 'idle' | 'rendering' | 'completed' | 'failed';

/**
 * 渲染状态结果接口
 */
interface RenderStatusResult {
  status: RenderStatus;
  progress: number; // 0-1
  currentFrame?: number;
  totalFrames?: number;
  outputPath?: string;
  error?: string;
  estimatedTimeRemaining?: number; // 秒
  startTime?: string;
  endTime?: string;
}

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 获取视频渲染状态
 */
async function handleGetRenderStatus(
  input: z.infer<typeof GetRenderStatusInputSchema>
): Promise<string> {
  try {
    const { workspacePath, projectName, renderId } = input;

    // 获取项目路径
    const projectPath = getProjectPath(workspacePath, projectName);
    if (!existsSync(projectPath)) {
      throw new Error(`项目不存在: ${projectName}`);
    }

    // 读取渲染状态文件
    const statusFilePath = join(projectPath, '.render-status.json');

    let result: RenderStatusResult;

    if (existsSync(statusFilePath)) {
      // 读取状态文件
      const statusData = JSON.parse(readFileSync(statusFilePath, 'utf-8'));

      // 如果指定了 renderId，筛选特定任务
      if (renderId && statusData.renderId !== renderId) {
        result = {
          status: 'idle',
          progress: 0,
        };
      } else {
        result = statusData;

        // 计算预估剩余时间
        if (result.status === 'rendering' && result.progress > 0) {
          const elapsed = Date.now() - new Date(result.startTime!).getTime();
          const totalEstimated = elapsed / result.progress;
          result.estimatedTimeRemaining = Math.round((totalEstimated - elapsed) / 1000);
        }
      }
    } else {
      // 没有状态文件，检查输出目录是否有视频文件
      const outputPath = getOutputPath(workspacePath, projectName);
      const hasOutput = existsSync(outputPath);

      result = {
        status: hasOutput ? 'completed' : 'idle',
        progress: hasOutput ? 1 : 0,
      };
    }

    return JSON.stringify(createSuccessResponse(result));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 更新渲染状态（内部使用）
 *
 * 这个函数由渲染引擎调用来更新状态
 */
export function updateRenderStatus(
  workspacePath: string,
  projectName: string,
  status: Partial<RenderStatusResult>
): void {
  try {
    const projectPath = getProjectPath(workspacePath, projectName);
    const statusFilePath = join(projectPath, '.render-status.json');

    // 读取现有状态
    let currentStatus: RenderStatusResult = {
      status: 'idle',
      progress: 0,
    };

    if (existsSync(statusFilePath)) {
      currentStatus = JSON.parse(readFileSync(statusFilePath, 'utf-8'));
    }

    // 合并新状态
    const newStatus = {
      ...currentStatus,
      ...status,
    };

    // 写入状态文件
    const fs = require('fs');
    fs.writeFileSync(statusFilePath, JSON.stringify(newStatus, null, 2), 'utf-8');
  } catch (error) {
    console.error('更新渲染状态失败:', error);
  }
}

/**
 * 清除渲染状态（内部使用）
 */
export function clearRenderStatus(workspacePath: string, projectName: string): void {
  try {
    const projectPath = getProjectPath(workspacePath, projectName);
    const statusFilePath = join(projectPath, '.render-status.json');

    if (existsSync(statusFilePath)) {
      const fs = require('fs');
      fs.unlinkSync(statusFilePath);
    }
  } catch (error) {
    console.error('清除渲染状态失败:', error);
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册渲染状态工具到 FastMCP 服务器
 */
export function registerRenderStatusTools(mcp: FastMCP): void {
  mcp.addTool({
    name: 'video_get_render_status',
    description: '获取视频渲染的当前状态和进度',
    parameters: GetRenderStatusInputSchema,
    execute: handleGetRenderStatus,
  });
}

// 导出处理函数供测试使用
export { handleGetRenderStatus };
