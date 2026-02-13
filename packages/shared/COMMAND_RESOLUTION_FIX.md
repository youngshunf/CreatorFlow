# Command Resolution Fix for MCP Stdio Validation

## 问题描述

在桌面应用中测试 video-mcp source 时，出现错误：
```
Command not found: "bun". Install the required dependency and try again.
```

这是因为 `validateStdioMcpConnection()` 函数直接使用配置中的命令名（如 "bun"），而没有通过 `resolveCommand()` 解析为绝对路径。

## 根本原因

1. `server-builder.ts` 中的 `buildMcpServer()` 函数已经实现了 `resolveCommand()` 来自动探测命令路径
2. 但是 `validation.ts` 中的 `validateStdioMcpConnection()` 函数没有使用这个机制
3. 导致在 Electron 子进程中（PATH 环境变量受限）无法找到 bun 命令

## 修复方案

### 1. 导出 `resolveCommand` 函数

**文件**: `packages/shared/src/sources/server-builder.ts`

将 `resolveCommand` 从私有函数改为导出函数：

```typescript
// 修改前
function resolveCommand(command: string): string {

// 修改后
export function resolveCommand(command: string): string {
```

### 2. 在验证函数中使用 `resolveCommand`

**文件**: `packages/shared/src/mcp/validation.ts`

#### 2.1 导入 `resolveCommand`

```typescript
import { resolveCommand } from '../sources/server-builder.ts';
```

#### 2.2 在函数开始时解析命令

```typescript
export async function validateStdioMcpConnection(
  config: StdioValidationConfig
): Promise<McpValidationResult> {
  const { command, args = [], env = {}, timeout = 30000 } = config;

  // 解析命令路径（自动探测 bun、node 等）
  const resolvedCommand = resolveCommand(command);
  debug(`[stdio-validation] Resolved command: ${command} → ${resolvedCommand}`);
  debug(`[stdio-validation] Spawning: ${resolvedCommand} ${args.join(' ')}`);

  // ... 后续代码
}
```

#### 2.3 在 spawn 和 transport 中使用解析后的命令

```typescript
// spawn 进程时使用 resolvedCommand
childProcess = spawn(resolvedCommand, args, {
  env: { ...process.env, ...env },
  stdio: ['pipe', 'pipe', 'pipe'],
});

// 创建 transport 时使用 resolvedCommand
const transport = new StdioClientTransport({
  command: resolvedCommand,
  args,
  env: { ...processEnv, ...env },
});
```

## 验证测试

### 测试 1: `resolveCommand` 函数

```bash
bun run test-resolve-command.ts
```

结果：
```
✓ resolveCommand('bun') → /Users/mac/.bun/bin/bun
✓ resolveCommand('node') → /Users/mac/.nvm/versions/node/v22.22.0/bin/node
✓ resolveCommand('/usr/bin/env') → /usr/bin/env
✓ resolveCommand('nonexistent-command-xyz') → nonexistent-command-xyz
```

### 测试 2: `validateStdioMcpConnection` 函数

```bash
bun run test-stdio-validation.ts
```

结果：
- ✅ 命令解析成功（不再报 "Command not found: bun"）
- ❌ video-mcp 服务器启动失败（代码错误，与本次修复无关）

## 影响范围

### 修改的文件

1. `packages/shared/src/sources/server-builder.ts` - 导出 `resolveCommand` 函数
2. `packages/shared/src/mcp/validation.ts` - 使用 `resolveCommand` 解析命令

### 受益的功能

1. **MCP Source 连接测试** - 在桌面应用的 Sources 页面测试 stdio MCP 连接
2. **跨平台兼容性** - 自动探测命令路径，支持不同的安装位置
3. **开发体验** - 不需要手动配置绝对路径

## 命令探测策略

`resolveCommand()` 函数使用以下策略：

1. **绝对路径检查** - 如果已经是绝对路径，直接返回
2. **which 命令查找** - 扩展 PATH 包括常见目录：
   - `~/.bun/bin`
   - `~/.deno/bin`
   - `/usr/local/bin`
   - `/opt/homebrew/bin`
3. **Well-known 位置检查** - 检查常见安装位置：
   - bun: `~/.bun/bin/bun`, `/opt/homebrew/bin/bun`, `/usr/local/bin/bun`
   - node: `/usr/local/bin/node`, `/opt/homebrew/bin/node`, `~/.nvm/versions/node/current/bin/node`
   - deno: `~/.deno/bin/deno`, `/usr/local/bin/deno`, `/opt/homebrew/bin/deno`

## 下一步

1. ✅ 修复 video-mcp 服务器的代码错误（`mcp.tool is not a function`）
2. ✅ 在桌面应用中测试 video-mcp source 连接
3. ✅ 验证所有 19 个工具是否可用
4. ✅ 更新 video-mcp 文档

## 相关文件

- `packages/shared/src/sources/server-builder.ts` - 命令解析逻辑
- `packages/shared/src/mcp/validation.ts` - MCP 连接验证
- `packages/shared/test-resolve-command.ts` - 命令解析测试脚本
- `packages/shared/test-stdio-validation.ts` - stdio 验证测试脚本
- `apps/electron/dist/main.cjs` - 已重新构建（包含修复）

## 构建状态

- ✅ 主进程已重新构建（`bun run build:main`）
- ✅ 修改已包含在 `dist/main.cjs` 中
- ✅ 桌面应用已启动（`bun run dev`）
- ⏳ 等待用户在应用中测试 video-mcp 连接
