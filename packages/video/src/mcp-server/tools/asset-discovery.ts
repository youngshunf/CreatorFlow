/**
 * 素材发现工具
 *
 * MCP 工具：video_list_available_assets
 *
 * 列出工作区中可用的素材文件，帮助 Agent 发现可以使用的素材
 */

import { z } from "zod";
import type { FastMCP } from "fastmcp";
import { readdirSync, statSync } from "fs";
import { join, extname, basename } from "path";
import { existsSync } from "fs";
import { createSuccessResponse } from "../types";
import { toErrorResponse } from "../types/errors";
import { SUPPORTED_ASSET_EXTENSIONS } from "../../types";

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_list_available_assets 输入 Schema
 */
export const ListAvailableAssetsInputSchema = z.object({
  workspacePath: z.string().describe("工作区根路径"),
  assetType: z
    .enum(["image", "video", "audio", "font", "all"])
    .optional()
    .describe("素材类型筛选（可选）"),
  searchPattern: z.string().optional().describe('搜索模式（如 "*.png"）'),
  maxResults: z.number().optional().default(100).describe("最大返回数量"),
});

/**
 * 素材信息接口
 */
interface AssetInfo {
  path: string;
  name: string;
  type: "image" | "video" | "audio" | "font" | "unknown";
  size: number;
  dimensions?: {
    width: number;
    height: number;
  };
}

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 列出工作区中可用的素材文件
 */
async function handleListAvailableAssets(
  rawInput: z.input<typeof ListAvailableAssetsInputSchema>,
): Promise<string> {
  try {
    const {
      workspacePath,
      assetType,
      searchPattern,
      maxResults = 100,
    } = rawInput;

    // 检查工作区是否存在
    if (!existsSync(workspacePath)) {
      throw new Error(`工作区不存在: ${workspacePath}`);
    }

    // 搜索素材文件
    const assets: AssetInfo[] = [];
    const seenPaths = new Set<string>(); // 用于去重
    const searchDirs = [
      join(workspacePath, "assets"),
      join(workspacePath, "素材"),
      workspacePath, // 也搜索根目录
    ];

    for (const dir of searchDirs) {
      if (existsSync(dir)) {
        scanDirectory(
          dir,
          assets,
          assetType,
          searchPattern,
          maxResults,
          seenPaths,
        );
      }
    }

    // 限制返回数量
    const limitedAssets = assets.slice(0, maxResults);

    return JSON.stringify(
      createSuccessResponse({
        assets: limitedAssets,
        total: assets.length,
        returned: limitedAssets.length,
      }),
    );
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 递归扫描目录查找素材
 */
function scanDirectory(
  dir: string,
  assets: AssetInfo[],
  assetType?: "image" | "video" | "audio" | "font" | "all",
  searchPattern?: string,
  maxResults?: number,
  seenPaths?: Set<string>,
): void {
  // 如果已达到最大数量，停止扫描
  if (maxResults && assets.length >= maxResults) {
    return;
  }

  try {
    const entries = readdirSync(dir, { withFileTypes: true });

    for (const entry of entries) {
      // 如果已达到最大数量，停止
      if (maxResults && assets.length >= maxResults) {
        break;
      }

      const fullPath = join(dir, entry.name);

      // 跳过隐藏文件和特定目录
      if (entry.name.startsWith(".") || entry.name === "node_modules") {
        continue;
      }

      if (entry.isDirectory()) {
        // 递归扫描子目录（限制深度为 3 层）
        const depth = fullPath.split("/").length - dir.split("/").length;
        if (depth < 3) {
          scanDirectory(
            fullPath,
            assets,
            assetType,
            searchPattern,
            maxResults,
            seenPaths,
          );
        }
      } else if (entry.isFile()) {
        // 检查文件是否匹配
        const assetInfo = getAssetInfo(fullPath);
        if (assetInfo && matchesFilter(assetInfo, assetType, searchPattern)) {
          // 使用 Set 去重
          if (!seenPaths || !seenPaths.has(fullPath)) {
            assets.push(assetInfo);
            if (seenPaths) {
              seenPaths.add(fullPath);
            }
          }
        }
      }
    }
  } catch (error) {
    // 忽略无法访问的目录
    console.warn(`无法扫描目录 ${dir}:`, error);
  }
}

/**
 * 获取素材信息
 */
function getAssetInfo(filePath: string): AssetInfo | null {
  try {
    const ext = extname(filePath).toLowerCase();
    const type = getAssetType(ext);

    // 如果不是支持的素材类型，返回 null
    if (type === "unknown") {
      return null;
    }

    const stats = statSync(filePath);

    const assetInfo: AssetInfo = {
      path: filePath,
      name: basename(filePath),
      type,
      size: stats.size,
    };

    // TODO: 如果是图片或视频，可以获取尺寸信息
    // 这需要额外的库（如 sharp, ffprobe）
    // 暂时不实现，避免增加依赖

    return assetInfo;
  } catch (error) {
    return null;
  }
}

/**
 * 根据扩展名判断素材类型
 */
function getAssetType(ext: string): AssetInfo["type"] {
  if (SUPPORTED_ASSET_EXTENSIONS.image.includes(ext as any)) {
    return "image";
  }
  if (SUPPORTED_ASSET_EXTENSIONS.video.includes(ext as any)) {
    return "video";
  }
  if (SUPPORTED_ASSET_EXTENSIONS.audio.includes(ext as any)) {
    return "audio";
  }
  if (SUPPORTED_ASSET_EXTENSIONS.font.includes(ext as any)) {
    return "font";
  }
  return "unknown";
}

/**
 * 检查素材是否匹配筛选条件
 */
function matchesFilter(
  asset: AssetInfo,
  assetType?: "image" | "video" | "audio" | "font" | "all",
  searchPattern?: string,
): boolean {
  // 类型筛选
  if (assetType && assetType !== "all" && asset.type !== assetType) {
    return false;
  }

  // 搜索模式筛选
  if (searchPattern) {
    const pattern = searchPattern.replace(/\*/g, ".*").replace(/\?/g, ".");
    const regex = new RegExp(pattern, "i");
    if (!regex.test(asset.name)) {
      return false;
    }
  }

  return true;
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册素材发现工具到 FastMCP 服务器
 */
export function registerAssetDiscoveryTools(mcp: FastMCP): void {
  mcp.addTool({
    name: "video_list_available_assets",
    description: "列出工作区中可用的素材文件",
    parameters: ListAvailableAssetsInputSchema,
    execute: handleListAvailableAssets,
  });
}

// 导出处理函数供测试使用
export { handleListAvailableAssets };
