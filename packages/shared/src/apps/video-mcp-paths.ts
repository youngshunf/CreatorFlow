/**
 * 视频 MCP 服务器路径解析
 *
 * 处理开发环境和打包后应用的路径差异
 */

import { join } from 'path';
import { existsSync } from 'fs';

/**
 * 获取 Bun 可执行文件路径
 *
 * @param isPackaged - 是否为打包后的应用
 * @returns Bun 可执行文件路径
 */
export function getBunPath(isPackaged: boolean): string {
  if (isPackaged) {
    // 打包后，使用系统安装的 Bun
    // 假设 Bun 已安装在系统 PATH 中
    return 'bun';
  } else {
    // 开发环境，使用系统 Bun
    return 'bun';
  }
}

/**
 * 获取视频 MCP 服务器入口文件路径
 *
 * @param isPackaged - 是否为打包后的应用
 * @param resourcesPath - 打包后的资源路径（app.getPath('resources')）
 * @returns MCP 服务器入口文件的绝对路径
 */
export function getVideoMcpServerPath(isPackaged: boolean, resourcesPath?: string): string {
  if (isPackaged) {
    // 打包后的路径
    // 资源被打包到 app.asar 或 resources/app/ 目录
    if (!resourcesPath) {
      throw new Error('resourcesPath is required for packaged app');
    }

    // 尝试多个可能的路径
    const possiblePaths = [
      // 方案 1: resources/app/packages/video/src/mcp-server/index.ts
      join(resourcesPath, 'app', 'packages', 'video', 'src', 'mcp-server', 'index.ts'),
      // 方案 2: resources/app.asar/packages/video/src/mcp-server/index.ts
      join(resourcesPath, 'app.asar', 'packages', 'video', 'src', 'mcp-server', 'index.ts'),
      // 方案 3: resources/packages/video/src/mcp-server/index.ts
      join(resourcesPath, 'packages', 'video', 'src', 'mcp-server', 'index.ts'),
    ];

    // 返回第一个存在的路径
    for (const path of possiblePaths) {
      if (existsSync(path)) {
        return path;
      }
    }

    // 如果都不存在，返回默认路径（可能在 asar 中）
    return possiblePaths[0];
  } else {
    // 开发环境路径
    // 从当前文件位置推导项目根目录
    // 当前文件: packages/shared/src/apps/video-mcp-paths.ts
    // 目标文件: packages/video/src/mcp-server/index.ts

    // 使用 __dirname 或相对路径
    // 这里假设从 shared 包到 video 包的相对路径
    const devPath = join(__dirname, '..', '..', '..', 'video', 'src', 'mcp-server', 'index.ts');

    if (existsSync(devPath)) {
      return devPath;
    }

    // 备用方案：使用绝对路径（需要在构建时替换）
    throw new Error(`Video MCP server not found at ${devPath}`);
  }
}

/**
 * 解析配置中的路径占位符
 *
 * 将配置文件中的占位符替换为实际路径：
 * - {{BUN_PATH}} -> Bun 可执行文件路径
 * - {{VIDEO_MCP_SERVER_PATH}} -> MCP 服务器入口文件路径
 *
 * @param configValue - 配置值（可能包含占位符）
 * @param isPackaged - 是否为打包后的应用
 * @param resourcesPath - 打包后的资源路径
 * @returns 解析后的配置值
 */
export function resolveConfigPath(
  configValue: string | string[],
  isPackaged: boolean,
  resourcesPath?: string
): string | string[] {
  const resolve = (value: string): string => {
    return value
      .replace('{{BUN_PATH}}', getBunPath(isPackaged))
      .replace('{{VIDEO_MCP_SERVER_PATH}}', getVideoMcpServerPath(isPackaged, resourcesPath));
  };

  if (Array.isArray(configValue)) {
    return configValue.map(resolve);
  } else {
    return resolve(configValue);
  }
}

/**
 * 解析数据源配置中的路径占位符
 *
 * @param sourceConfig - 数据源配置对象
 * @param isPackaged - 是否为打包后的应用
 * @param resourcesPath - 打包后的资源路径
 * @returns 解析后的配置对象
 */
export function resolveSourceConfigPaths(
  sourceConfig: any,
  isPackaged: boolean,
  resourcesPath?: string
): any {
  const resolved = { ...sourceConfig };

  // 解析 MCP 配置中的路径
  if (resolved.mcp) {
    if (resolved.mcp.command) {
      resolved.mcp.command = resolveConfigPath(resolved.mcp.command, isPackaged, resourcesPath) as string;
    }
    if (resolved.mcp.args) {
      resolved.mcp.args = resolveConfigPath(resolved.mcp.args, isPackaged, resourcesPath) as string[];
    }
  }

  return resolved;
}
