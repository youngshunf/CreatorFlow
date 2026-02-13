# 新增 MCP 服务器指南

本文档以 `video-mcp` 为参考，说明如何在 CreatorFlow 中新增一个 MCP 服务器。

## 目录结构

假设新增一个名为 `my-tool` 的 MCP 服务器：

```
packages/my-tool/
├── package.json
├── my-tool-mcp-source-config.json    # 源配置预设
└── src/
    └── mcp-server/
        ├── index.ts                   # 入口文件（CLI 参数解析 + 启动）
        ├── server.ts                  # FastMCP 实例创建 + 传输配置
        ├── tools/
        │   ├── index.ts              # 工具注册总入口
        │   └── my-feature.ts         # 具体工具实现
        ├── types/
        │   ├── index.ts              # 响应类型
        │   └── errors.ts             # 错误类型
        └── utils/
            └── stdio-compat-patch.ts  # stdio 协议兼容补丁（可从 video 复制）
```

## 步骤

### 1. 创建包和 package.json

```json
{
  "name": "@creator-flow/my-tool",
  "version": "0.1.0",
  "type": "module",
  "main": "src/index.ts",
  "bin": {
    "creator-flow-my-tool-mcp": "./src/mcp-server/index.ts"
  },
  "exports": {
    ".": "./src/index.ts",
    "./mcp-server": "./src/mcp-server/index.ts"
  },
  "dependencies": {
    "fastmcp": "^3.32.0",
    "zod": "^3.24.0"
  }
}
```

然后在项目根目录运行 `bun install`。

### 2. 实现 MCP 服务器

#### 2.1 server.ts — 创建 FastMCP 实例

```typescript
import { FastMCP } from 'fastmcp';
import { registerAllTools } from './tools';

export interface ServerConfig {
  transport: 'stdio' | 'http';
  host?: string;
  port?: number;
  endpoint?: `/${string}`;
}

export const DEFAULT_CONFIG: ServerConfig = {
  transport: 'stdio',
  host: '0.0.0.0',
  port: 3001,
  endpoint: '/mcp',
};

export function createServer(): FastMCP {
  const mcp = new FastMCP({
    name: 'creator-flow-my-tool',
    version: '0.1.0',
    instructions: '你的 MCP 服务器描述',
  });
  registerAllTools(mcp);
  return mcp;
}

export async function runServer(config: ServerConfig = DEFAULT_CONFIG): Promise<void> {
  const mcp = createServer();

  if (config.transport === 'stdio') {
    await mcp.start({ transportType: 'stdio' });
  } else {
    await mcp.start({
      transportType: 'httpStream',
      httpStream: { port: config.port!, endpoint: config.endpoint! },
    });
  }
}
```

#### 2.2 tools/index.ts — 注册工具

```typescript
import type { FastMCP } from 'fastmcp';
import { registerMyFeatureTools } from './my-feature';

export function registerAllTools(mcp: FastMCP): void {
  // 所有日志必须用 console.error，不能用 console.log
  // stdout 是 MCP 协议通道，任何非 JSON-RPC 输出都会破坏协议
  console.error('[MCP My Tool] Registering tools...');
  registerMyFeatureTools(mcp);
  console.error('[MCP My Tool] All tools registered');
}
```

#### 2.3 tools/my-feature.ts — 实现具体工具

```typescript
import { z } from 'zod';
import type { FastMCP } from 'fastmcp';

const MyInputSchema = z.object({
  query: z.string().describe('搜索关键词'),
});

export function registerMyFeatureTools(mcp: FastMCP): void {
  mcp.addTool({
    name: 'my_tool_search',           // 工具名前缀用下划线风格
    description: '搜索某某内容',
    parameters: MyInputSchema,
    execute: async (input) => {
      // 实现逻辑
      return JSON.stringify({ results: [] });
    },
  });
}
```

> **工具命名规则**: 工具名必须匹配 `/^[a-zA-Z0-9_.-]{1,64}$/`（Anthropic API 要求）。建议用 `{服务名}_{动作}` 格式，如 `my_tool_search`。

#### 2.4 index.ts — 入口文件

```typescript
import { runServer, DEFAULT_CONFIG, type ServerConfig } from './server';
import { applyStdioCompatPatch } from './utils/stdio-compat-patch';

// 必须在 runServer 之前调用
applyStdioCompatPatch();

function parseArgs(): ServerConfig | null {
  const args = process.argv.slice(2);
  const config: ServerConfig = { ...DEFAULT_CONFIG };
  for (let i = 0; i < args.length; i++) {
    switch (args[i]) {
      case '--transport': config.transport = args[++i] as 'stdio' | 'http'; break;
      case '--port': config.port = parseInt(args[++i]!, 10); break;
      case '--help': return null;
    }
  }
  return config;
}

async function main(): Promise<void> {
  const config = parseArgs();
  if (!config) { process.exit(0); }
  try {
    await runServer(config);
  } catch (error) {
    console.error('[MCP My Tool] Failed to start:', error);
    process.exit(1);
  }
}

main();
```

#### 2.5 复制 stdio-compat-patch.ts

从 `packages/video/src/mcp-server/utils/stdio-compat-patch.ts` 复制到你的 `utils/` 目录。这个补丁让服务器同时支持 Content-Length 头协议和换行分隔 JSON 协议。

### 3. 创建源配置预设

在包根目录创建 `my-tool-mcp-source-config.json`：

```json
{
  "id": "my-tool-mcp_preset",
  "name": "我的工具",
  "slug": "my-tool-mcp",
  "enabled": true,
  "provider": "sprouty-my-tool",
  "type": "mcp",
  "mcp": {
    "command": "node",
    "args": ["index.cjs"],
    "cwd": "app:resources/my-tool-mcp-server",
    "transport": "stdio",
    "env": {}
  },
  "icon": "🔧",
  "tagline": "一句话描述",
  "description": "详细描述",
  "categories": ["tools"],
  "tags": ["my-tool"]
}
```

关键字段说明：
- `slug` — 唯一标识，agent 调用工具时的前缀（`mcp__{slug}__{tool_name}`）
- `cwd: "app:resources/..."` — `app:` 前缀表示相对于 Electron 应用资源目录
- `args: ["index.cjs"]` — 构建产物的文件名（会被 server-builder 解析为绝对路径）

### 4. 添加 esbuild 构建配置

编辑 `scripts/electron-build-main.ts`：

```typescript
// 1. 在文件顶部添加路径常量
const MY_TOOL_MCP_ENTRY = join(ROOT_DIR, "packages/my-tool/src/mcp-server/index.ts");
const MY_TOOL_MCP_OUTPUT = join(ROOT_DIR, "apps/electron/resources/my-tool-mcp-server/index.cjs");

// 2. 添加构建函数
async function buildMyToolMcpServer(): Promise<void> {
  console.log("🔧 Building My Tool MCP Server...");

  const outDir = join(MY_TOOL_MCP_OUTPUT, "..");
  if (!existsSync(outDir)) {
    mkdirSync(outDir, { recursive: true });
  }

  const proc = spawn({
    cmd: [
      "bun", "run", "esbuild",
      MY_TOOL_MCP_ENTRY,
      "--bundle",
      "--platform=node",
      "--format=cjs",
      "--outfile=" + MY_TOOL_MCP_OUTPUT,
      // import.meta.url polyfill（CJS 环境必需）
      "--banner:js=#!/usr/bin/env node\nif(typeof globalThis.__import_meta_url__==='undefined'){const{pathToFileURL}=require('url');globalThis.__import_meta_url__=pathToFileURL(__filename).href;}",
      "--define:import.meta.url=globalThis.__import_meta_url__",
      // fastmcp 的可选依赖
      "--external:@valibot/to-json-schema",
      "--external:effect",
      "--external:sury",
      // 如果有其他不需要打包的依赖，也加 --external
    ],
    cwd: ROOT_DIR,
    stdout: "inherit",
    stderr: "inherit",
  });

  const exitCode = await proc.exited;
  if (exitCode !== 0) {
    console.error("❌ My Tool MCP server build failed");
    process.exit(exitCode);
  }
  console.log("✅ My Tool MCP server built successfully");
}

// 3. 在 main() 函数中调用
await buildMyToolMcpServer();
```

### 5. 注册源配置到应用

将源配置预设文件复制或链接到工作区的 sources 目录。源配置最终存储在：

```
~/.sprouty-ai/workspaces/{workspaceId}/sources/{slug}/config.json
```

应用启动时通过 `loadWorkspaceSources()` 加载所有源，然后 `SourceServerBuilder.buildAll()` 构建 MCP 服务器配置传递给 Agent SDK。

## 注意事项

### stdout 纯净性（最重要）

stdio 模式下，stdout 是 MCP JSON-RPC 协议通道。**任何非 JSON-RPC 输出都会导致握手失败，工具列表为空。**

- 所有日志必须用 `console.error`（输出到 stderr）
- 不要用 `console.log`
- 第三方库如果往 stdout 写内容，需要拦截或禁用

### CLI 不支持 cwd

Claude CLI 的 Zod schema 不包含 `cwd` 字段，传递过去会被 strip。`server-builder.ts` 已有 workaround：自动将 `args` 中的相对路径解析为绝对路径。所以源配置中的 `cwd` + 相对路径 `args` 组合是可以正常工作的。

### 工具名命名

- 格式：`/^[a-zA-Z0-9_.-]{1,64}$/`
- 建议：`{服务前缀}_{动作}`，如 `video_list_templates`
- Agent 调用时的完整名称：`mcp__{slug}__{tool_name}`

### 构建产物位置

构建输出到 `apps/electron/resources/{server-name}/index.cjs`，与源配置中的 `cwd: "app:resources/{server-name}"` 对应。

### 测试方法

```bash
# 1. 构建
bun run scripts/electron-build-main.ts

# 2. 独立测试 stdio 协议
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0.0"}}}' | node apps/electron/resources/my-tool-mcp-server/index.cjs

# 3. 验证输出是纯 JSON（无 banner、无日志）
# 期望只看到 {"jsonrpc":"2.0","id":1,"result":{...}}

# 4. 测试 tools/list
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0.0"}}}
{"jsonrpc":"2.0","method":"notifications/initialized"}
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}' | node apps/electron/resources/my-tool-mcp-server/index.cjs
```

## 快速检查清单

- [ ] `package.json` 中声明了 `fastmcp` 依赖
- [ ] 所有日志用 `console.error`，不用 `console.log`
- [ ] 入口文件调用了 `applyStdioCompatPatch()`
- [ ] 源配置 JSON 中 `slug`、`type: "mcp"`、`transport: "stdio"` 正确
- [ ] `electron-build-main.ts` 中添加了构建函数和调用
- [ ] 构建产物路径与源配置 `cwd` + `args` 匹配
- [ ] 独立测试 stdout 输出纯 JSON-RPC
- [ ] 工具名符合 `[a-zA-Z0-9_.-]{1,64}` 规则
