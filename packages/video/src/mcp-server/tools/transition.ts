/**
 * 过渡效果工具
 *
 * MCP 工具：video_set_transitions
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore } from '../services/project-store';
import { TransitionTypeEnum, TransitionDirectionEnum, createSuccessResponse } from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

const TransitionItemSchema = z.object({
  type: TransitionTypeEnum.describe('过渡类型: fade | slide | wipe | flip | clock-wipe | none'),
  durationInFrames: z.number().int().positive().describe('过渡时长（帧数）'),
  direction: TransitionDirectionEnum.optional().describe('过渡方向（仅 slide/wipe/flip 需要）: from-left | from-right | from-top | from-bottom'),
});

export const SetTransitionsInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  transitions: z.array(TransitionItemSchema).describe('过渡效果列表（长度必须等于场景数 - 1）'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

async function handleSetTransitions(
  input: z.infer<typeof SetTransitionsInputSchema>,
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const project = await store.setTransitions(input.projectId, input.transitions);
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

export function registerTransitionTools(mcp: FastMCP): void {
  mcp.addTool({
    name: 'video_set_transitions',
    description: `设置场景间的过渡效果。

过渡效果列表的长度必须等于场景数 - 1（即每两个相邻场景之间一个过渡）。

可用的过渡类型：
- fade: 淡入淡出
- slide: 滑动（需指定 direction）
- wipe: 擦除（需指定 direction）
- flip: 翻转（需指定 direction）
- clock-wipe: 时钟擦除
- none: 无过渡（直接切换）

方向选项：from-left, from-right, from-top, from-bottom`,
    parameters: SetTransitionsInputSchema,
    execute: handleSetTransitions,
  });
}

export { handleSetTransitions };
