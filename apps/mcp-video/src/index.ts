#!/usr/bin/env bun
/**
 * MCP Video Server 入口文件
 *
 * 支持两种传输模式：
 * - stdio: 用于桌面 MCP 客户端（如 Claude Desktop、Kiro）
 * - http: 用于远程/云端部署
 *
 * 使用方式：
 * - stdio 模式: bun run start
 * - http 模式: bun run start --transport http --port 3000
 *
 * 命令行参数：
 * --transport <stdio|http>  传输模式（默认: stdio）
 * --host <host>             HTTP 模式的主机地址（默认: 0.0.0.0）
 * --port <port>             HTTP 模式的端口（默认: 3000）
 * --endpoint <path>         HTTP 模式的端点路径（默认: /mcp）
 * --help                    显示帮助信息
 *
 * @requirements 8.5 - 启动时打印启动信息
 */

import { runServer, DEFAULT_CONFIG, type ServerConfig } from './server';

// ============================================================================
// 命令行参数解析
// ============================================================================

/**
 * 显示帮助信息
 */
function showHelp(): void {
  console.log(`
CreatorFlow Video MCP Server

Usage:
  bun run src/index.ts [options]

Options:
  --transport <mode>   Transport mode: stdio or http (default: stdio)
  --host <host>        Host address for HTTP mode (default: 0.0.0.0)
  --port <port>        Port for HTTP mode (default: 3000)
  --endpoint <path>    Endpoint path for HTTP mode (default: /mcp)
  --help               Show this help message

Examples:
  # Start in stdio mode (for Claude Desktop, Kiro)
  bun run src/index.ts

  # Start in HTTP mode on port 3000
  bun run src/index.ts --transport http --port 3000

  # Start in HTTP mode on custom host and port
  bun run src/index.ts --transport http --host 127.0.0.1 --port 8080

Configuration for Claude Desktop (claude_desktop_config.json):
  {
    "mcpServers": {
      "creator-flow-video": {
        "command": "bun",
        "args": ["run", "/path/to/CreatorFlow/apps/mcp-video/src/index.ts"]
      }
    }
  }

Configuration for Kiro:
  Add to your MCP server configuration with stdio transport.
`);
}

/**
 * 解析命令行参数
 *
 * @returns 服务器配置，如果显示帮助则返回 null
 */
function parseArgs(): ServerConfig | null {
  const args = process.argv.slice(2);
  const config: ServerConfig = { ...DEFAULT_CONFIG };

  for (let i = 0; i < args.length; i++) {
    const arg = args[i];

    if (arg === undefined) {
      continue;
    }

    switch (arg) {
      case '--help':
      case '-h':
        showHelp();
        return null;

      case '--transport':
      case '-t': {
        const transport = args[++i];
        if (transport !== 'stdio' && transport !== 'http') {
          console.error(`Error: Invalid transport mode "${transport}". Use "stdio" or "http".`);
          process.exit(1);
        }
        config.transport = transport;
        break;
      }

      case '--host': {
        const host = args[++i];
        if (!host) {
          console.error('Error: --host requires a value');
          process.exit(1);
        }
        config.host = host;
        break;
      }

      case '--port':
      case '-p': {
        const portStr = args[++i];
        if (!portStr) {
          console.error('Error: --port requires a value');
          process.exit(1);
        }
        const port = parseInt(portStr, 10);
        if (isNaN(port) || port < 1 || port > 65535) {
          console.error(`Error: Invalid port "${portStr}". Must be a number between 1 and 65535.`);
          process.exit(1);
        }
        config.port = port;
        break;
      }

      case '--endpoint': {
        const endpoint = args[++i];
        if (!endpoint) {
          console.error('Error: --endpoint requires a value');
          process.exit(1);
        }
        if (!endpoint.startsWith('/')) {
          console.error('Error: --endpoint must start with /');
          process.exit(1);
        }
        config.endpoint = endpoint as `/${string}`;
        break;
      }

      default:
        if (arg.startsWith('-')) {
          console.error(`Error: Unknown option "${arg}". Use --help for usage information.`);
          process.exit(1);
        }
    }
  }

  return config;
}

// ============================================================================
// 主入口
// ============================================================================

/**
 * 主函数
 */
async function main(): Promise<void> {
  // 解析命令行参数
  const config = parseArgs();

  // 如果显示帮助，直接退出
  if (config === null) {
    process.exit(0);
  }

  try {
    // 启动服务器
    await runServer(config);
  } catch (error) {
    console.error('[MCP Video Server] Failed to start server:', error);
    process.exit(1);
  }
}

// 运行主函数
main().catch((error) => {
  console.error('[MCP Video Server] Unexpected error:', error);
  process.exit(1);
});
