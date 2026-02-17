/**
 * Video Creation Tools for AI Agents
 *
 * This module defines video creation tools that can be used by AI agents
 * to create, edit, preview, and render videos programmatically.
 *
 * @requirements 10.1, 10.2, 10.3, 10.4, 10.5, 10.6
 */

import { z } from "zod";
import type { VideoProject, VideoConfig, Asset, AssetType } from "../types";
import type {
  VideoTemplate,
  TemplateProps,
  ColorScheme,
} from "../templates/types";

// ============================================================================
// Tool Input Schemas
// ============================================================================

/**
 * Input schema for video_create_project tool
 * @requirements 10.1
 */
export const VideoCreateProjectInputSchema = z.object({
  /** Project name */
  name: z.string().describe("Name of the video project"),
  /** Template to use (optional) */
  template: z
    .enum([
      "blank",
      "title-animation",
      "slideshow",
      "product-showcase",
      "social-media",
      "social-media-vertical",
      "social-media-square",
      "marketing",
      "product-marketing",
      "promo-ad",
      "tutorial",
      "step-by-step",
      "explainer",
      "tips",
    ])
    .optional()
    .describe("Template to use for the project"),
  /** Video width in pixels */
  width: z
    .number()
    .int()
    .positive()
    .optional()
    .default(1920)
    .describe("Video width in pixels"),
  /** Video height in pixels */
  height: z
    .number()
    .int()
    .positive()
    .optional()
    .default(1080)
    .describe("Video height in pixels"),
  /** Frames per second */
  fps: z
    .number()
    .int()
    .positive()
    .optional()
    .default(30)
    .describe("Frames per second"),
  /** Duration in seconds */
  durationInSeconds: z
    .number()
    .positive()
    .describe("Duration of the video in seconds"),
  /** Project description */
  description: z.string().optional().describe("Optional project description"),
});

export type VideoCreateProjectInput = z.infer<
  typeof VideoCreateProjectInputSchema
>;

/**
 * Input schema for video_generate_composition tool
 * @requirements 10.2
 */
export const VideoGenerateCompositionInputSchema = z.object({
  /** Project ID */
  projectId: z.string().describe("ID of the video project"),
  /** Description of what to generate */
  description: z
    .string()
    .describe("Natural language description of the video content"),
  /** Visual style */
  style: z
    .enum([
      "modern",
      "minimal",
      "playful",
      "corporate",
      "cinematic",
      "vibrant",
      "nature",
    ])
    .optional()
    .default("modern")
    .describe("Visual style for the composition"),
  /** Color scheme */
  colorScheme: z
    .object({
      primary: z.string().describe("Primary color (hex)"),
      secondary: z.string().describe("Secondary color (hex)"),
      background: z.string().describe("Background color (hex)"),
    })
    .optional()
    .describe("Custom color scheme"),
});

export type VideoGenerateCompositionInput = z.infer<
  typeof VideoGenerateCompositionInputSchema
>;

/**
 * Input schema for video_update_composition tool
 * @requirements 10.3
 */
export const VideoUpdateCompositionInputSchema = z.object({
  /** Project ID */
  projectId: z.string().describe("ID of the video project"),
  /** Composition ID */
  compositionId: z.string().describe("ID of the composition to update"),
  /** Changes to apply */
  changes: z
    .string()
    .describe("Natural language description of changes to make"),
  /** New props to merge */
  props: z
    .record(z.string(), z.any())
    .optional()
    .describe("Props to update on the composition"),
});

export type VideoUpdateCompositionInput = z.infer<
  typeof VideoUpdateCompositionInputSchema
>;

/**
 * Input schema for video_preview tool
 * @requirements 10.4
 */
export const VideoPreviewInputSchema = z.object({
  /** Project ID */
  projectId: z.string().describe("ID of the video project to preview"),
  /** Composition ID (optional, defaults to first) */
  compositionId: z
    .string()
    .optional()
    .describe("ID of specific composition to preview"),
  /** Start frame (optional) */
  startFrame: z
    .number()
    .int()
    .min(0)
    .optional()
    .describe("Frame to start preview from"),
});

export type VideoPreviewInput = z.infer<typeof VideoPreviewInputSchema>;

/**
 * Input schema for video_render tool
 * @requirements 10.5
 */
export const VideoRenderInputSchema = z.object({
  /** Project ID */
  projectId: z.string().describe("ID of the video project to render"),
  /** Composition ID */
  compositionId: z.string().describe("ID of the composition to render"),
  /** Output format */
  outputFormat: z
    .enum(["mp4", "webm", "gif"])
    .optional()
    .default("mp4")
    .describe("Output video format"),
  /** Quality preset */
  quality: z
    .enum(["draft", "standard", "high"])
    .optional()
    .default("standard")
    .describe("Quality preset for rendering"),
  /** Output path (optional, will prompt user if not provided) */
  outputPath: z.string().optional().describe("Path to save the rendered video"),
});

export type VideoRenderInput = z.infer<typeof VideoRenderInputSchema>;

/**
 * Input schema for video_add_asset tool
 * @requirements 10.6
 */
export const VideoAddAssetInputSchema = z.object({
  /** Project ID */
  projectId: z.string().describe("ID of the video project"),
  /** Path to the asset file */
  assetPath: z.string().describe("Path to the asset file to add"),
  /** Type of asset */
  assetType: z
    .enum(["image", "video", "audio", "font"])
    .describe("Type of the asset"),
  /** Optional custom name */
  name: z.string().optional().describe("Custom name for the asset"),
});

export type VideoAddAssetInput = z.infer<typeof VideoAddAssetInputSchema>;

// ============================================================================
// Tool Definitions
// ============================================================================

/**
 * Tool definition interface
 */
export interface VideoTool<TInput, TOutput> {
  /** Tool name */
  name: string;
  /** Tool description */
  description: string;
  /** Input schema */
  inputSchema: z.ZodType<TInput>;
  /** Execute the tool */
  execute: (input: TInput) => Promise<TOutput>;
}

/**
 * Video tool names
 */
export const VIDEO_TOOL_NAMES = {
  CREATE_PROJECT: "video_create_project",
  GENERATE_COMPOSITION: "video_generate_composition",
  UPDATE_COMPOSITION: "video_update_composition",
  PREVIEW: "video_preview",
  RENDER: "video_render",
  ADD_ASSET: "video_add_asset",
} as const;

/**
 * Tool definitions for AI agents
 */
export const VIDEO_TOOLS = {
  /**
   * Create a new video project
   * @requirements 10.1
   */
  [VIDEO_TOOL_NAMES.CREATE_PROJECT]: {
    name: VIDEO_TOOL_NAMES.CREATE_PROJECT,
    description:
      "Create a new video project with optional template. Use this to start a new video creation workflow.",
    inputSchema: VideoCreateProjectInputSchema,
    parameters: {
      type: "object",
      properties: {
        name: { type: "string", description: "Name of the video project" },
        template: {
          type: "string",
          enum: [
            "blank",
            "title-animation",
            "slideshow",
            "product-showcase",
            "social-media",
            "social-media-vertical",
            "social-media-square",
            "marketing",
            "product-marketing",
            "promo-ad",
            "tutorial",
            "step-by-step",
            "explainer",
            "tips",
          ],
          description: "Template to use for the project",
        },
        width: {
          type: "number",
          description: "Video width in pixels (default: 1920)",
        },
        height: {
          type: "number",
          description: "Video height in pixels (default: 1080)",
        },
        fps: { type: "number", description: "Frames per second (default: 30)" },
        durationInSeconds: {
          type: "number",
          description: "Duration of the video in seconds",
        },
        description: {
          type: "string",
          description: "Optional project description",
        },
      },
      required: ["name", "durationInSeconds"],
    },
  },

  /**
   * Generate a composition from description
   * @requirements 10.2
   */
  [VIDEO_TOOL_NAMES.GENERATE_COMPOSITION]: {
    name: VIDEO_TOOL_NAMES.GENERATE_COMPOSITION,
    description:
      "Generate a video composition from a natural language description. The AI will create appropriate animations and visuals.",
    inputSchema: VideoGenerateCompositionInputSchema,
    parameters: {
      type: "object",
      properties: {
        projectId: { type: "string", description: "ID of the video project" },
        description: {
          type: "string",
          description: "Natural language description of the video content",
        },
        style: {
          type: "string",
          enum: [
            "modern",
            "minimal",
            "playful",
            "corporate",
            "cinematic",
            "vibrant",
            "nature",
          ],
          description: "Visual style for the composition",
        },
        colorScheme: {
          type: "object",
          properties: {
            primary: { type: "string", description: "Primary color (hex)" },
            secondary: { type: "string", description: "Secondary color (hex)" },
            background: {
              type: "string",
              description: "Background color (hex)",
            },
          },
          description: "Custom color scheme",
        },
      },
      required: ["projectId", "description"],
    },
  },

  /**
   * Update an existing composition
   * @requirements 10.3
   */
  [VIDEO_TOOL_NAMES.UPDATE_COMPOSITION]: {
    name: VIDEO_TOOL_NAMES.UPDATE_COMPOSITION,
    description:
      "Update an existing video composition. Can modify text, colors, animations, or other properties.",
    inputSchema: VideoUpdateCompositionInputSchema,
    parameters: {
      type: "object",
      properties: {
        projectId: { type: "string", description: "ID of the video project" },
        compositionId: {
          type: "string",
          description: "ID of the composition to update",
        },
        changes: {
          type: "string",
          description: "Natural language description of changes to make",
        },
        props: {
          type: "object",
          description: "Props to update on the composition",
        },
      },
      required: ["projectId", "compositionId", "changes"],
    },
  },

  /**
   * Preview a video project
   * @requirements 10.4
   */
  [VIDEO_TOOL_NAMES.PREVIEW]: {
    name: VIDEO_TOOL_NAMES.PREVIEW,
    description:
      "Start a preview server for a video project. Returns a URL that can be opened in a browser.",
    inputSchema: VideoPreviewInputSchema,
    parameters: {
      type: "object",
      properties: {
        projectId: {
          type: "string",
          description: "ID of the video project to preview",
        },
        compositionId: {
          type: "string",
          description: "ID of specific composition to preview",
        },
        startFrame: {
          type: "number",
          description: "Frame to start preview from",
        },
      },
      required: ["projectId"],
    },
  },

  /**
   * Render a video
   * @requirements 10.5
   */
  [VIDEO_TOOL_NAMES.RENDER]: {
    name: VIDEO_TOOL_NAMES.RENDER,
    description:
      "Render a video composition to a file. Supports MP4, WebM, and GIF formats with different quality presets.",
    inputSchema: VideoRenderInputSchema,
    parameters: {
      type: "object",
      properties: {
        projectId: {
          type: "string",
          description: "ID of the video project to render",
        },
        compositionId: {
          type: "string",
          description: "ID of the composition to render",
        },
        outputFormat: {
          type: "string",
          enum: ["mp4", "webm", "gif"],
          description: "Output video format (default: mp4)",
        },
        quality: {
          type: "string",
          enum: ["draft", "standard", "high"],
          description: "Quality preset for rendering (default: standard)",
        },
        outputPath: {
          type: "string",
          description: "Path to save the rendered video",
        },
      },
      required: ["projectId", "compositionId"],
    },
  },

  /**
   * Add an asset to a project
   * @requirements 10.6
   */
  [VIDEO_TOOL_NAMES.ADD_ASSET]: {
    name: VIDEO_TOOL_NAMES.ADD_ASSET,
    description:
      "Add an asset (image, video, audio, or font) to a video project for use in compositions.",
    inputSchema: VideoAddAssetInputSchema,
    parameters: {
      type: "object",
      properties: {
        projectId: { type: "string", description: "ID of the video project" },
        assetPath: {
          type: "string",
          description: "Path to the asset file to add",
        },
        assetType: {
          type: "string",
          enum: ["image", "video", "audio", "font"],
          description: "Type of the asset",
        },
        name: { type: "string", description: "Custom name for the asset" },
      },
      required: ["projectId", "assetPath", "assetType"],
    },
  },
} as const;

/**
 * Get all video tool definitions
 */
export function getVideoTools() {
  return Object.values(VIDEO_TOOLS);
}

/**
 * Get a specific video tool by name
 */
export function getVideoTool(name: string) {
  return VIDEO_TOOLS[name as keyof typeof VIDEO_TOOLS];
}
