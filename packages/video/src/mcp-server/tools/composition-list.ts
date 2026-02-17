/**
 * 组合列表工具
 *
 * MCP 工具：video_list_available_compositions
 *
 * 列出所有内置组合及其 props schema，供 Agent 选择使用
 */

import { z } from "zod";
import type { FastMCP } from "fastmcp";
import { COMPOSITION_IDS, type CompositionId } from "../../compositions";
import { createSuccessResponse } from "../types";
import { toErrorResponse } from "../types/errors";

// ============================================================================
// 组合注册表
// ============================================================================

interface CompositionInfo {
  id: CompositionId;
  name: string;
  description: string;
  category: "基础" | "社交媒体" | "教程" | "营销";
  propsSchema: Record<
    string,
    {
      type: string;
      description: string;
      required?: boolean;
      default?: any;
    }
  >;
  defaultProps: Record<string, any>;
}

/**
 * 所有内置组合的注册信息
 */
const COMPOSITION_REGISTRY: CompositionInfo[] = [
  {
    id: "TitleAnimation",
    name: "标题动画",
    description:
      "显示标题和副标题的动画效果，支持 fade/slide/scale/spring 四种风格",
    category: "基础",
    propsSchema: {
      title: { type: "string", description: "主标题文字", default: "Welcome" },
      subtitle: {
        type: "string",
        description: "副标题文字",
        default: "Your subtitle here",
      },
      animationStyle: {
        type: "enum(fade|slide|scale|spring)",
        description: "动画风格",
        default: "spring",
      },
      titleFontSize: {
        type: "number",
        description: "标题字号（像素）",
        default: 72,
      },
      subtitleFontSize: {
        type: "number",
        description: "副标题字号（像素）",
        default: 36,
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
      logo: { type: "string", description: "Logo 图片路径" },
    },
    defaultProps: {
      title: "Welcome",
      subtitle: "Your subtitle here",
      animationStyle: "spring",
      titleFontSize: 72,
      subtitleFontSize: 36,
    },
  },
  {
    id: "Slideshow",
    name: "图片幻灯片",
    description: "图片轮播展示，支持 fade/slide/zoom/crossfade 过渡效果",
    category: "基础",
    propsSchema: {
      items: {
        type: "array",
        description: "幻灯片列表 [{ title, description?, image }]",
        required: true,
      },
      animationStyle: {
        type: "enum(fade|slide|zoom)",
        description: "切换动画",
        default: "fade",
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {
      animationStyle: "fade",
    },
  },
  {
    id: "DataVisualization",
    name: "数据图表",
    description: "数据可视化图表，支持 bar/line/pie/donut 四种图表类型",
    category: "基础",
    propsSchema: {
      chartType: {
        type: "enum(bar|line|pie|donut)",
        description: "图表类型",
        default: "bar",
      },
      title: { type: "string", description: "图表标题" },
      items: {
        type: "array",
        description: "数据点 [{ title, description?, value }]",
        required: true,
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {
      chartType: "bar",
    },
  },
  {
    id: "ProductShowcase",
    name: "产品展示",
    description: "产品展示页面，支持 centered/split/features-grid 三种布局",
    category: "基础",
    propsSchema: {
      title: { type: "string", description: "产品名称" },
      subtitle: { type: "string", description: "产品描述" },
      items: {
        type: "array",
        description: "功能特性 [{ title, description, icon? }]",
      },
      logo: { type: "string", description: "产品图片路径" },
      animationStyle: {
        type: "enum(fade|slide|scale)",
        description: "布局风格",
        default: "fade",
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {
      animationStyle: "fade",
    },
  },
  {
    id: "SocialMediaVertical",
    name: "竖版社交媒体",
    description:
      "9:16 竖版社交媒体视频（TikTok、Instagram Reels、YouTube Shorts）",
    category: "社交媒体",
    propsSchema: {
      title: { type: "string", description: "主标题" },
      subtitle: { type: "string", description: "副标题" },
      items: { type: "array", description: "内容项列表" },
      cta: { type: "object", description: "行动号召 { text, url? }" },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
      logo: { type: "string", description: "Logo 图片路径" },
    },
    defaultProps: {},
  },
  {
    id: "SocialMediaSquare",
    name: "方形社交媒体",
    description: "1:1 方形社交媒体视频（Instagram Feed、Facebook）",
    category: "社交媒体",
    propsSchema: {
      title: { type: "string", description: "主标题" },
      subtitle: { type: "string", description: "副标题" },
      items: { type: "array", description: "内容项列表" },
      cta: { type: "object", description: "行动号召 { text, url? }" },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
      logo: { type: "string", description: "Logo 图片路径" },
    },
    defaultProps: {},
  },
  {
    id: "StepByStepTutorial",
    name: "分步教程",
    description: "分步骤教程视频，每个步骤带标题、描述和可选图片",
    category: "教程",
    propsSchema: {
      title: { type: "string", description: "教程标题" },
      items: {
        type: "array",
        description: "步骤列表 [{ title, description, image? }]",
        required: true,
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {},
  },
  {
    id: "Explainer",
    name: "概念讲解",
    description: "概念讲解视频，适合解释复杂概念或流程",
    category: "教程",
    propsSchema: {
      title: { type: "string", description: "讲解标题" },
      subtitle: { type: "string", description: "副标题" },
      items: {
        type: "array",
        description: "讲解要点 [{ title, description, icon? }]",
        required: true,
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {},
  },
  {
    id: "Tips",
    name: "技巧清单",
    description: "技巧/建议清单视频，逐条展示要点",
    category: "教程",
    propsSchema: {
      title: { type: "string", description: "清单标题" },
      items: {
        type: "array",
        description: "技巧列表 [{ title, description, icon? }]",
        required: true,
      },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
    },
    defaultProps: {},
  },
  {
    id: "ProductMarketing",
    name: "产品营销",
    description: "产品营销推广视频，展示产品特性和卖点",
    category: "营销",
    propsSchema: {
      title: { type: "string", description: "产品名称" },
      subtitle: { type: "string", description: "产品标语" },
      items: {
        type: "array",
        description: "产品特性 [{ title, description, icon? }]",
      },
      cta: { type: "object", description: "行动号召 { text, url? }" },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
      logo: { type: "string", description: "产品图片路径" },
    },
    defaultProps: {},
  },
  {
    id: "PromoAd",
    name: "促销广告",
    description: "促销广告视频，突出优惠信息和行动号召",
    category: "营销",
    propsSchema: {
      title: { type: "string", description: "促销标题" },
      subtitle: { type: "string", description: "促销描述" },
      cta: { type: "object", description: "行动号召 { text, url? }" },
      colors: {
        type: "object",
        description: "颜色配置 { primary, secondary, background, text }",
      },
      logo: { type: "string", description: "Logo 图片路径" },
    },
    defaultProps: {},
  },
];

// ============================================================================
// Zod Schemas
// ============================================================================

export const ListAvailableCompositionsInputSchema = z
  .object({})
  .describe("无需输入参数");

// ============================================================================
// 工具处理函数
// ============================================================================

async function handleListAvailableCompositions(): Promise<string> {
  try {
    // 按分类分组
    const byCategory: Record<string, CompositionInfo[]> = {};
    for (const comp of COMPOSITION_REGISTRY) {
      if (!byCategory[comp.category]) {
        byCategory[comp.category] = [];
      }
      byCategory[comp.category]!.push(comp);
    }

    return JSON.stringify(
      createSuccessResponse({
        total: COMPOSITION_REGISTRY.length,
        compositions: COMPOSITION_REGISTRY,
        byCategory,
      }),
    );
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

export function registerCompositionListTools(mcp: FastMCP): void {
  mcp.addTool({
    name: "video_list_available_compositions",
    description: `列出所有可用的内置组合及其参数说明。

每个组合是一个可复用的视频"积木块"，Agent 通过 video_add_scene 将组合添加到项目中，
并通过 props 配置其内容和外观。

返回每个组合的 ID、名称、描述、分类、参数 schema 和默认值。`,
    parameters: ListAvailableCompositionsInputSchema,
    execute: handleListAvailableCompositions,
  });
}

export { handleListAvailableCompositions };
