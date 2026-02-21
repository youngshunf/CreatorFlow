/**
 * Session-Scoped Video Tools
 *
 * Agent 可用的视频创作工具，直接操作 SQLite 数据库。
 * 运行在 Electron main 进程中（通过 createSdkMcpServer in-process）。
 */

import { tool } from '@anthropic-ai/claude-agent-sdk';
import { z } from 'zod';
import { getCachedConnection, type CreatorMediaDB } from '../db/connection.ts';
import { getWorkspaceDbPath } from '../workspaces/storage.ts';
import * as videoRepo from '../db/repositories/video.ts';
import * as contentsRepo from '../db/repositories/contents.ts';
import * as projectsRepo from '../db/repositories/projects.ts';
import {
  generateContentDirPath,
  ensureContentDirs,
} from '../utils/content-paths.ts';
import type {
  CreateVideoProject,
  CreateVideoScene,
  CreateVideoAsset,
  UpdateVideoScene,
  VideoProjectFull,
  CreateContent,
} from '../db/types.ts';
import { basename } from 'path';
import { debug } from '../utils/debug.ts';

// ============================================================
// Helper
// ============================================================

function getDB(workspaceRootPath: string): CreatorMediaDB {
  const dbPath = getWorkspaceDbPath(workspaceRootPath);
  const db = getCachedConnection(dbPath);
  if (!db) {
    throw new Error(
      `数据库连接未初始化。请确保工作区已正确加载。(workspace: ${workspaceRootPath})`,
    );
  }
  return db;
}

function generateId(): string {
  return crypto.randomUUID();
}

function textResult(text: string) {
  return { content: [{ type: 'text' as const, text }] };
}

function jsonResult(data: unknown) {
  return textResult(JSON.stringify(data, null, 2));
}

function errorResult(message: string) {
  return { content: [{ type: 'text' as const, text: `Error: ${message}` }], isError: true };
}

/** 确保视频项目有关联的 content 记录，没有则自动创建并生成目录结构 */
function ensureContentForVideoProject(
  db: CreatorMediaDB,
  projectId: string,
  contentId: string,
  title: string,
  workspaceRootPath: string,
): string {
  if (contentId && contentId.trim() !== '') return contentId;

  // 获取项目名称
  let projectName = title;
  const project = projectsRepo.getProject(db, projectId);
  if (project) projectName = project.name;

  // 获取下一个序号
  const nextNumber = contentsRepo.getNextContentNumber(db, projectId);
  const contentDirPath = generateContentDirPath(projectName, nextNumber, title);

  const newContentId = generateId();
  const contentData: CreateContent = {
    id: newContentId,
    project_id: projectId,
    title,
    status: 'creating',
    target_platforms: null,
    pipeline_mode: 'manual',
    content_dir_path: contentDirPath,
    viral_pattern_id: null,
    metadata: null,
    content_tracks: 'video',
  };
  contentsRepo.createContent(db, contentData);

  // 创建文件系统目录
  ensureContentDirs(workspaceRootPath, contentDirPath, 'video');

  return newContentId;
}

// ============================================================
// Composition Info (static, no DB needed)
// ============================================================

const COMPOSITION_INFO: Record<string, { name: string; description: string }> = {
  TitleAnimation: { name: '标题动画', description: '大标题 + 副标题入场动画，适合视频开头' },
  Slideshow: { name: '幻灯片', description: '多张图片轮播，支持过渡效果' },
  DataVisualization: { name: '数据可视化', description: '图表和数据展示动画' },
  ProductShowcase: { name: '产品展示', description: '产品特性展示，支持多个特性卡片' },
  SocialMediaVertical: { name: '社交媒体竖版', description: '9:16 竖版，适合抖音/快手/Reels' },
  SocialMediaSquare: { name: '社交媒体方形', description: '1:1 方形，适合 Instagram/微博' },
  StepByStepTutorial: { name: '分步教程', description: '逐步教学，每步有标题和说明' },
  Explainer: { name: '解说视频', description: '概念解释，支持多个要点展示' },
  Tips: { name: '技巧分享', description: '多条技巧/建议展示' },
  ProductMarketing: { name: '产品营销', description: '产品推广视频，突出卖点' },
  PromoAd: { name: '促销广告', description: '短促销广告，强调 CTA' },
};

// ============================================================
// Tool Factories
// ============================================================

export function createVideoListCompositionsTool() {
  return tool(
    'video_list_compositions',
    `列出所有可用的视频组合（Composition）及其说明。

组合是视频的基本构建块，每个场景使用一个组合来渲染。
调用此工具了解有哪些组合可用，然后在 video_add_scene 中使用 compositionId。`,
    {},
    async () => {
      const compositions = Object.entries(COMPOSITION_INFO).map(([id, info]) => ({
        id,
        ...info,
      }));
      return jsonResult(compositions);
    },
  );
}

export function createVideoCreateProjectTool(workspaceRootPath: string) {
  return tool(
    'video_create_project',
    `创建一个新的视频项目，关联到指定的内容（content）。

每个视频项目包含多个场景（scenes），场景按顺序排列并可设置过渡效果。
创建后使用 video_add_scene 添加场景。

如果不传 contentId，会自动创建一个 content 记录并关联。
默认分辨率 1080x1920（竖版 9:16），30fps。`,
    {
      projectId: z.string().describe('项目 ID（projects 表的 id，用于创建 content 记录）'),
      contentId: z.string().optional().describe('关联的内容 ID（contents 表的 id，不传则自动创建）'),
      name: z.string().describe('项目名称'),
      description: z.string().optional().describe('项目描述'),
      width: z.number().optional().describe('视频宽度（默认 1080）'),
      height: z.number().optional().describe('视频高度（默认 1920）'),
      fps: z.number().optional().describe('帧率（默认 30）'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        const contentId = ensureContentForVideoProject(
          db,
          args.projectId,
          args.contentId ?? '',
          args.name,
          workspaceRootPath,
        );

        const data: CreateVideoProject = {
          id: generateId(),
          content_id: contentId,
          name: args.name,
          description: args.description ?? null,
          width: args.width ?? 1080,
          height: args.height ?? 1920,
          fps: args.fps ?? 30,
          metadata: null,
        };
        const project = videoRepo.createVideoProject(db, data);

        // 回填 content metadata
        contentsRepo.updateContent(db, contentId, {
          metadata: JSON.stringify({ video_project_id: project.id }),
        });

        const full: VideoProjectFull = { ...project, scenes: [], assets: [] };
        return jsonResult(full);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoGetProjectTool(workspaceRootPath: string) {
  return tool(
    'video_get_project',
    `获取视频项目的完整信息，包括所有场景和素材。

返回项目基本信息、场景列表（按顺序）和素材列表。`,
    {
      projectId: z.string().describe('视频项目 ID'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        const project = videoRepo.getVideoProjectFull(db, args.projectId);
        if (!project) return errorResult(`项目不存在: ${args.projectId}`);
        return jsonResult(project);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoAddSceneTool(workspaceRootPath: string) {
  return tool(
    'video_add_scene',
    `向视频项目添加一个场景。

场景使用 compositionId 指定渲染组件（先用 video_list_compositions 查看可用组合）。
props 是 JSON 字符串，包含该组合的参数（如标题、颜色、图片等）。
可选设置过渡效果（fade/slide/wipe/flip/clock-wipe）。

场景默认追加到末尾，也可通过 insertAt 指定插入位置。`,
    {
      projectId: z.string().describe('视频项目 ID'),
      compositionId: z.string().describe('组合 ID（如 TitleAnimation, Slideshow 等）'),
      name: z.string().optional().describe('场景名称'),
      durationInFrames: z.number().optional().describe('场景时长（帧数，默认 90，即 30fps 下 3 秒）'),
      props: z.string().optional().describe('组合参数 JSON 字符串（默认 "{}"）'),
      transitionType: z.enum(['none', 'fade', 'slide', 'wipe', 'flip', 'clock-wipe']).optional()
        .describe('过渡效果类型（默认 none）'),
      transitionDuration: z.number().optional().describe('过渡时长（帧数，默认 0）'),
      transitionDirection: z.enum(['from-left', 'from-right', 'from-top', 'from-bottom']).optional()
        .describe('过渡方向（仅 slide/wipe/flip 有效）'),
      insertAt: z.number().optional().describe('插入位置索引（从 0 开始，默认追加到末尾）'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);

        // 验证 compositionId
        if (!COMPOSITION_INFO[args.compositionId]) {
          const validIds = Object.keys(COMPOSITION_INFO).join(', ');
          return errorResult(`无效的 compositionId: ${args.compositionId}。可用: ${validIds}`);
        }

        // 计算 sort_order
        let sortOrder: number;
        if (args.insertAt !== undefined) {
          const scenes = videoRepo.listVideoScenes(db, args.projectId);
          db.transaction(() => {
            for (const s of scenes) {
              if (s.sort_order >= args.insertAt!) {
                videoRepo.updateVideoScene(db, s.id, { sort_order: s.sort_order + 1 });
              }
            }
          });
          sortOrder = args.insertAt;
        } else {
          sortOrder = videoRepo.getNextSortOrder(db, args.projectId);
        }

        const data: CreateVideoScene = {
          id: generateId(),
          project_id: args.projectId,
          composition_id: args.compositionId,
          name: args.name ?? null,
          sort_order: sortOrder,
          duration_in_frames: args.durationInFrames ?? 90,
          props: args.props ?? '{}',
          transition_type: args.transitionType ?? 'none',
          transition_duration: args.transitionDuration ?? 0,
          transition_direction: args.transitionDirection ?? null,
        };

        const scene = videoRepo.addVideoScene(db, data);
        return jsonResult({ sceneId: scene.id, scene });
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoUpdateSceneTool(workspaceRootPath: string) {
  return tool(
    'video_update_scene',
    `更新视频场景的属性。

可更新：名称、时长、组合参数（props）、过渡效果等。
只传需要更新的字段即可。`,
    {
      sceneId: z.string().describe('场景 ID'),
      name: z.string().optional().describe('新名称'),
      durationInFrames: z.number().optional().describe('新时长（帧数）'),
      props: z.string().optional().describe('新的组合参数 JSON 字符串'),
      compositionId: z.string().optional().describe('更换组合 ID'),
      transitionType: z.enum(['none', 'fade', 'slide', 'wipe', 'flip', 'clock-wipe']).optional()
        .describe('过渡效果类型'),
      transitionDuration: z.number().optional().describe('过渡时长（帧数）'),
      transitionDirection: z.enum(['from-left', 'from-right', 'from-top', 'from-bottom']).optional()
        .describe('过渡方向'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        const updates: UpdateVideoScene = {};
        if (args.name !== undefined) updates.name = args.name;
        if (args.durationInFrames !== undefined) updates.duration_in_frames = args.durationInFrames;
        if (args.props !== undefined) updates.props = args.props;
        if (args.compositionId !== undefined) updates.composition_id = args.compositionId;
        if (args.transitionType !== undefined) updates.transition_type = args.transitionType;
        if (args.transitionDuration !== undefined) updates.transition_duration = args.transitionDuration;
        if (args.transitionDirection !== undefined) updates.transition_direction = args.transitionDirection;

        const scene = videoRepo.updateVideoScene(db, args.sceneId, updates);
        if (!scene) return errorResult(`场景不存在: ${args.sceneId}`);
        return jsonResult(scene);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoRemoveSceneTool(workspaceRootPath: string) {
  return tool(
    'video_remove_scene',
    `从视频项目中删除一个场景。`,
    {
      sceneId: z.string().describe('要删除的场景 ID'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        const removed = videoRepo.removeVideoScene(db, args.sceneId);
        if (!removed) return errorResult(`场景不存在: ${args.sceneId}`);
        return textResult(`场景 ${args.sceneId} 已删除`);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoReorderScenesTool(workspaceRootPath: string) {
  return tool(
    'video_reorder_scenes',
    `重新排列视频项目中场景的顺序。

传入场景 ID 数组，按新顺序排列。数组中的第一个场景排在最前面。`,
    {
      projectId: z.string().describe('视频项目 ID'),
      sceneIds: z.array(z.string()).describe('按新顺序排列的场景 ID 数组'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        videoRepo.reorderVideoScenes(db, args.projectId, args.sceneIds);
        const scenes = videoRepo.listVideoScenes(db, args.projectId);
        return jsonResult({ message: '场景已重新排序', scenes });
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoAddAssetTool(workspaceRootPath: string) {
  return tool(
    'video_add_asset',
    `向视频项目添加素材文件（图片、视频、音频、字体）。

素材添加后可在场景的 props 中引用其 file_path。`,
    {
      projectId: z.string().describe('视频项目 ID'),
      filePath: z.string().describe('素材文件的绝对路径'),
      type: z.enum(['image', 'video', 'audio', 'font']).describe('素材类型'),
      name: z.string().optional().describe('素材名称（默认使用文件名）'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);

        const data: CreateVideoAsset = {
          id: generateId(),
          project_id: args.projectId,
          type: args.type,
          name: args.name ?? basename(args.filePath),
          file_path: args.filePath,
          file_size: null,
          metadata: null,
        };

        const asset = videoRepo.addVideoAsset(db, data);
        return jsonResult(asset);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoListTemplatesTool() {
  return tool(
    'video_list_templates',
    `列出所有可用的视频模板。

模板是预配置的视频项目，包含默认场景和参数。
可按分类筛选：social-media（社交媒体）、marketing（营销）、tutorial（教程）。

使用 video_create_from_template 从模板快速创建项目。`,
    {
      category: z.enum(['social-media', 'marketing', 'tutorial']).optional()
        .describe('按分类筛选（可选）'),
    },
    async (args) => {
      // 动态导入避免在非 Electron 环境下报错
      try {
        const { ALL_TEMPLATES, getTemplatesByCategory } = await import('@sprouty-ai/video/templates');
        const templates = args.category
          ? [...getTemplatesByCategory(args.category)]
          : [...ALL_TEMPLATES];
        return jsonResult(templates.map(t => ({
          id: t.id,
          name: t.name,
          description: t.description,
          category: t.category,
          compositionId: t.compositionId,
          aspectRatio: t.aspectRatio,
          defaultConfig: t.defaultConfig,
          tags: t.tags,
        })));
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

export function createVideoCreateFromTemplateTool(workspaceRootPath: string) {
  return tool(
    'video_create_from_template',
    `从模板快速创建视频项目。

先用 video_list_templates 查看可用模板，然后传入 templateId 创建。
模板会自动设置分辨率、帧率，并创建一个默认场景。
如果不传 contentId，会自动创建一个 content 记录并关联。`,
    {
      projectId: z.string().describe('项目 ID（projects 表的 id，用于创建 content 记录）'),
      contentId: z.string().optional().describe('关联的内容 ID（不传则自动创建）'),
      templateId: z.string().describe('模板 ID（从 video_list_templates 获取）'),
      name: z.string().describe('项目名称'),
      description: z.string().optional().describe('项目描述（默认使用模板描述）'),
    },
    async (args) => {
      try {
        const db = getDB(workspaceRootPath);
        const { getTemplateById } = await import('@sprouty-ai/video/templates');

        const template = getTemplateById(args.templateId);
        if (!template) return errorResult(`模板不存在: ${args.templateId}`);

        const contentId = ensureContentForVideoProject(
          db,
          args.projectId,
          args.contentId ?? '',
          args.name,
        );

        // 创建项目
        const projectData: CreateVideoProject = {
          id: generateId(),
          content_id: contentId,
          name: args.name,
          description: args.description ?? template.description ?? null,
          width: template.defaultConfig.width,
          height: template.defaultConfig.height,
          fps: template.defaultConfig.fps,
          metadata: JSON.stringify({ templateId: template.id }),
        };
        const project = videoRepo.createVideoProject(db, projectData);

        // 回填 content metadata
        contentsRepo.updateContent(db, contentId, {
          metadata: JSON.stringify({ video_project_id: project.id }),
        });

        // 创建默认场景
        const sceneData: CreateVideoScene = {
          id: generateId(),
          project_id: project.id,
          composition_id: template.compositionId,
          name: template.name,
          sort_order: 0,
          duration_in_frames: template.defaultConfig.durationInFrames,
          props: JSON.stringify(template.defaultProps),
          transition_type: 'none',
          transition_duration: 0,
          transition_direction: null,
        };
        videoRepo.addVideoScene(db, sceneData);

        const full = videoRepo.getVideoProjectFull(db, project.id)!;
        return jsonResult(full);
      } catch (e) {
        return errorResult(e instanceof Error ? e.message : String(e));
      }
    },
  );
}

// ============================================================
// Export all video tools
// ============================================================

/**
 * 创建所有视频相关的 session-scoped tools
 */
export function createVideoTools(workspaceRootPath: string) {
  return [
    createVideoListCompositionsTool(),
    createVideoCreateProjectTool(workspaceRootPath),
    createVideoGetProjectTool(workspaceRootPath),
    createVideoAddSceneTool(workspaceRootPath),
    createVideoUpdateSceneTool(workspaceRootPath),
    createVideoRemoveSceneTool(workspaceRootPath),
    createVideoReorderScenesTool(workspaceRootPath),
    createVideoAddAssetTool(workspaceRootPath),
    createVideoListTemplatesTool(),
    createVideoCreateFromTemplateTool(workspaceRootPath),
  ];
}
