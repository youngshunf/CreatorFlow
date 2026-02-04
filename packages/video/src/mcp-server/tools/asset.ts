/**
 * 素材管理工具
 *
 * MCP 工具：video_add_asset, video_remove_asset, video_list_assets
 *
 * @requirements 4.1, 4.2, 4.3, 4.4, 4.5
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { ProjectStore, type AssetsByType } from '../services/project-store';
import {
  type Asset,
  type AssetType,
  createSuccessResponse,
} from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_add_asset 输入 Schema
 */
export const AddAssetInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  assetPath: z.string().describe('素材文件路径（绝对路径）'),
  assetType: z.enum(['image', 'video', 'audio', 'font']).describe('素材类型'),
  name: z.string().optional().describe('素材名称（可选，默认使用文件名）'),
});

/**
 * video_remove_asset 输入 Schema
 */
export const RemoveAssetInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
  assetId: z.string().describe('素材ID'),
  deleteFile: z.boolean().optional().default(false).describe('是否删除实际文件'),
});

/**
 * video_list_assets 输入 Schema
 */
export const ListAssetsInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  projectId: z.string().describe('项目ID'),
});

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 添加素材到视频项目
 *
 * @requirements 4.1 - 复制文件到项目素材目录并注册
 * @requirements 4.2 - 验证文件类型
 * @requirements 4.3 - 文件不存在或格式不支持时返回错误
 */
async function handleAddAsset(
  input: z.infer<typeof AddAssetInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const asset = await store.addAsset(input.projectId, {
      sourcePath: input.assetPath,
      assetType: input.assetType as AssetType,
      name: input.name,
    });
    return JSON.stringify(createSuccessResponse(asset));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 从视频项目移除素材
 *
 * @requirements 4.4 - 移除素材注册，可选删除文件
 */
async function handleRemoveAsset(
  input: z.infer<typeof RemoveAssetInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const success = await store.removeAsset(
      input.projectId,
      input.assetId,
      input.deleteFile
    );
    return JSON.stringify(createSuccessResponse({ removed: success }));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 列出项目中的所有素材
 *
 * @requirements 4.5 - 返回按类型分组的素材列表
 */
async function handleListAssets(
  input: z.infer<typeof ListAssetsInputSchema>
): Promise<string> {
  try {
    const store = ProjectStore.create(input.workspacePath);
    const assets = await store.listAssets(input.projectId);
    return JSON.stringify(createSuccessResponse(assets));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册素材管理工具到 FastMCP 服务器
 */
export function registerAssetTools(mcp: FastMCP): void {
  // video_add_asset
  mcp.addTool({
    name: 'video_add_asset',
    description: `添加素材到视频项目。支持的格式：
- 图片: png, jpg, jpeg, gif, webp, svg
- 视频: mp4, webm, mov
- 音频: mp3, wav, ogg, m4a
- 字体: ttf, otf, woff, woff2

素材文件将被复制到项目的"素材"目录下。`,
    parameters: AddAssetInputSchema,
    execute: handleAddAsset,
  });

  // video_remove_asset
  mcp.addTool({
    name: 'video_remove_asset',
    description: '从视频项目移除素材。可选择是否同时删除实际文件。',
    parameters: RemoveAssetInputSchema,
    execute: handleRemoveAsset,
  });

  // video_list_assets
  mcp.addTool({
    name: 'video_list_assets',
    description: '列出项目中的所有素材，按类型（image、video、audio、font）分组返回。',
    parameters: ListAssetsInputSchema,
    execute: handleListAssets,
  });
}

// ============================================================================
// 导出
// ============================================================================

export {
  handleAddAsset,
  handleRemoveAsset,
  handleListAssets,
};
