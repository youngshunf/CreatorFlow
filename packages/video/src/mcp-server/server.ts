/**
 * FastMCP 服务器配置
 *
 * 配置 MCP 服务器实例，注册所有工具，支持 stdio 和 HTTP 传输模式
 *
 * @requirements 8.1 - 支持 stdio 传输模式
 * @requirements 8.2 - 支持 HTTP 传输模式
 * @requirements 8.3 - 工具具有正确的 JSON Schema 定义
 * @requirements 8.5 - 启动时记录可用工具和传输配置
 */

import { FastMCP } from 'fastmcp';
import { registerAllTools, TOOL_LIST } from './tools';

// ============================================================================
// 类型定义
// ============================================================================

/**
 * 服务器配置
 */
export interface ServerConfig {
  /** 传输模式：stdio 用于桌面客户端，http 用于远程部署 */
  transport: 'stdio' | 'http';
  /** HTTP 模式的主机地址 */
  host?: string;
  /** HTTP 模式的端口 */
  port?: number;
  /** HTTP 模式的端点路径 */
  endpoint?: `/${string}`;
}

/**
 * 默认服务器配置
 */
export const DEFAULT_CONFIG: ServerConfig = {
  transport: 'stdio',
  host: '0.0.0.0',
  port: 3000,
  endpoint: '/mcp',
};

// ============================================================================
// 服务器实例
// ============================================================================

/**
 * 服务器版本号
 */
export const SERVER_VERSION = '0.1.0' as const;

/**
 * 创建 FastMCP 服务器实例
 *
 * @returns FastMCP 服务器实例
 */
export function createServer(): FastMCP {
  const mcp = new FastMCP({
    name: 'creator-flow-video',
    version: SERVER_VERSION,
    instructions: 'CreatorFlow Video MCP Server - 提供视频创作能力，包括项目管理、素材管理、视频渲染和实时预览。',
  });

  // 注册所有工具
  registerAllTools(mcp);

  return mcp;
}

/**
 * 打印服务器启动信息
 *
 * @param config 服务器配置
 * @requirements 8.5 - 启动时记录可用工具和传输配置
 */
function printStartupInfo(config: ServerConfig): void {
  console.log('');
  console.log('╔════════════════════════════════════════════════════════════════╗');
  console.log('║           CreatorFlow Video MCP Server                         ║');
  console.log('╚════════════════════════════════════════════════════════════════╝');
  console.log('');

  // 打印传输配置
  console.log('📡 Transport Configuration:');
  console.log(`   Mode: ${config.transport.toUpperCase()}`);
  if (config.transport === 'http') {
    console.log(`   Host: ${config.host}`);
    console.log(`   Port: ${config.port}`);
    console.log(`   Endpoint: ${config.endpoint}`);
    console.log(`   URL:  http://${config.host}:${config.port}${config.endpoint}`);
  }
  console.log('');

  // 打印可用工具
  console.log('🔧 Available Tools:');

  // 按分类分组
  const toolsByCategory: Record<string, typeof TOOL_LIST[number][]> = {};
  for (const tool of TOOL_LIST) {
    const category = tool.category;
    if (!toolsByCategory[category]) {
      toolsByCategory[category] = [];
    }
    toolsByCategory[category]!.push(tool);
  }

  // 分类显示名称映射
  const categoryNames: Record<string, string> = {
    project: '📁 项目管理',
    asset: '🎨 素材管理',
    composition: '🎬 组合管理',
    render: '🎥 视频渲染',
    preview: '👁️ 实时预览',
    template: '📋 模板管理',
  };

  for (const [category, tools] of Object.entries(toolsByCategory)) {
    console.log(`   ${categoryNames[category] || category}:`);
    for (const tool of tools) {
      console.log(`     - ${tool.name}: ${tool.description}`);
    }
  }

  console.log('');
  console.log(`📊 Total: ${TOOL_LIST.length} tools registered`);
  console.log('');
  console.log('═══════════════════════════════════════════════════════════════════');
  console.log('');
}

/**
 * 启动 MCP 服务器
 *
 * @param config 服务器配置
 * @requirements 8.1 - 支持 stdio 传输模式
 * @requirements 8.2 - 支持 HTTP 传输模式
 */
export async function runServer(config: ServerConfig = DEFAULT_CONFIG): Promise<void> {
  // 打印启动信息
  printStartupInfo(config);

  // 创建服务器实例
  const mcp = createServer();

  // 根据传输模式启动服务器
  if (config.transport === 'stdio') {
    console.log('[MCP Video Server] Starting in stdio mode...');
    await mcp.start({
      transportType: 'stdio',
    });
  } else {
    const port = config.port ?? DEFAULT_CONFIG.port!;
    const endpoint = config.endpoint ?? DEFAULT_CONFIG.endpoint!;
    console.log(`[MCP Video Server] Starting HTTP server on ${config.host}:${port}${endpoint}...`);
    await mcp.start({
      transportType: 'httpStream',
      httpStream: {
        port,
        endpoint,
      },
    });
  }

  console.log('[MCP Video Server] Server started successfully');
}

// ============================================================================
// 导出
// ============================================================================

export { TOOL_LIST };
