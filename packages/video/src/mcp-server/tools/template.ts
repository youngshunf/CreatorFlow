/**
 * 模板工具
 *
 * MCP 工具：video_list_templates
 *
 * @requirements 10.1, 10.2, 10.3
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import {
  ALL_TEMPLATES,
  TEMPLATES_BY_CATEGORY,
  getTemplateById,
  getTemplatesByCategory,
  searchTemplatesByTags,
  type VideoTemplate,
  type TemplateCategory,
} from '../../templates';
import { createSuccessResponse } from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_list_templates 输入 Schema
 */
export const ListTemplatesInputSchema = z.object({
  category: z
    .enum(['social-media', 'marketing', 'tutorial'])
    .optional()
    .describe('模板分类（可选，不指定则返回所有模板）'),
  tags: z
    .array(z.string())
    .optional()
    .describe('标签过滤（可选，返回包含任一标签的模板）'),
});

/**
 * video_get_template 输入 Schema
 */
export const GetTemplateInputSchema = z.object({
  templateId: z.string().describe('模板ID'),
});

// ============================================================================
// 模板摘要类型
// ============================================================================

/**
 * 模板摘要（用于列表展示）
 */
interface TemplateSummary {
  id: string;
  name: string;
  description: string;
  category: string;
  aspectRatio?: string;
  tags?: string[];
  defaultConfig: {
    width: number;
    height: number;
    fps: number;
    durationInFrames: number;
  };
}

/**
 * 将完整模板转换为摘要
 */
function toTemplateSummary(template: VideoTemplate): TemplateSummary {
  return {
    id: template.id,
    name: template.name,
    description: template.description,
    category: template.category,
    aspectRatio: template.aspectRatio,
    tags: template.tags,
    defaultConfig: template.defaultConfig,
  };
}

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 列出所有可用的视频模板
 *
 * @requirements 10.1 - 返回所有可用模板及元数据
 * @requirements 10.3 - 支持 social-media、marketing、tutorial 分类
 */
async function handleListTemplates(
  input: z.infer<typeof ListTemplatesInputSchema>
): Promise<string> {
  try {
    let templates: readonly VideoTemplate[];

    // 根据分类过滤
    if (input.category) {
      templates = getTemplatesByCategory(input.category as TemplateCategory);
    } else if (input.tags && input.tags.length > 0) {
      // 根据标签过滤
      templates = searchTemplatesByTags(input.tags);
    } else {
      // 返回所有模板
      templates = ALL_TEMPLATES;
    }

    // 转换为摘要格式
    const summaries = templates.map(toTemplateSummary);

    // 按分类分组
    const grouped = {
      'social-media': summaries.filter((t) => t.category === 'social-media'),
      marketing: summaries.filter((t) => t.category === 'marketing'),
      tutorial: summaries.filter((t) => t.category === 'tutorial'),
    };

    return JSON.stringify(
      createSuccessResponse({
        total: summaries.length,
        templates: summaries,
        byCategory: grouped,
      })
    );
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 获取模板详情
 *
 * @requirements 10.2 - 返回模板的完整配置
 */
async function handleGetTemplate(
  input: z.infer<typeof GetTemplateInputSchema>
): Promise<string> {
  try {
    const template = getTemplateById(input.templateId);

    if (!template) {
      return JSON.stringify({
        success: false,
        error: {
          code: 'TEMPLATE_NOT_FOUND',
          message: `模板不存在: ${input.templateId}`,
          details: { field: 'templateId', received: input.templateId },
        },
      });
    }

    return JSON.stringify(createSuccessResponse(template));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册模板工具到 FastMCP 服务器
 */
export function registerTemplateTools(mcp: FastMCP): void {
  // video_list_templates
  mcp.addTool({
    name: 'video_list_templates',
    description: `列出所有可用的视频模板。

模板分类：
- social-media: 社交媒体模板（TikTok、Instagram Reels、YouTube Shorts 等）
- marketing: 营销推广模板（产品展示、促销广告等）
- tutorial: 教程模板（步骤教程、讲解视频、技巧分享等）

可以按分类或标签过滤模板。返回模板的基本信息和默认配置。`,
    parameters: ListTemplatesInputSchema,
    execute: handleListTemplates,
  });

  // video_get_template (额外工具，方便获取单个模板详情)
  mcp.addTool({
    name: 'video_get_template',
    description: '获取视频模板的完整详情，包括默认配置、默认属性和组合代码。',
    parameters: GetTemplateInputSchema,
    execute: handleGetTemplate,
  });
}

// ============================================================================
// 导出
// ============================================================================

export { handleListTemplates, handleGetTemplate };
