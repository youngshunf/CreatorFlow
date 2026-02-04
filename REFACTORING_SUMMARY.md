# 视频模块重构总结

## 执行时间
2026-02-04

## 重构内容

### 1. 目录结构调整

**之前的结构：**
```
apps/mcp-video/              # ❌ 位置不合理
├── src/
│   ├── index.ts
│   ├── server.ts
│   ├── tools/
│   ├── services/
│   └── ...
└── package.json

packages/video/              # ✅ 核心包
├── src/
│   ├── compositions/
│   ├── templates/
│   └── ...
└── package.json
```

**重构后的结构：**
```
packages/video/              # ✅ 统一的视频包
├── src/
│   ├── compositions/        # Remotion 组件
│   ├── templates/           # 视频模板
│   ├── components/          # 可复用 UI 组件
│   ├── hooks/               # React Hooks
│   ├── renderer/            # 渲染引擎
│   ├── skills/              # AI Agent Skills
│   ├── utils/               # 工具函数
│   └── mcp-server/          # ✅ MCP 服务器（新位置）
│       ├── index.ts         # MCP 入口
│       ├── server.ts        # FastMCP 配置
│       ├── tools/           # MCP 工具实现
│       ├── services/        # 核心服务
│       ├── types/           # 类型定义
│       ├── utils/           # MCP 工具函数
│       ├── README.md        # MCP 文档
│       └── tsconfig.json    # MCP TypeScript 配置
└── package.json             # 合并后的配置
```

### 2. package.json 更新

**新增内容：**

1. **bin 入口点**：
   ```json
   "bin": {
     "creator-flow-video-mcp": "./src/mcp-server/index.ts"
   }
   ```

2. **MCP 相关脚本**：
   ```json
   "scripts": {
     "mcp:dev": "bun run --hot src/mcp-server/index.ts",
     "mcp:start": "bun run src/mcp-server/index.ts",
     "mcp:start:http": "bun run src/mcp-server/index.ts --transport http",
     "typecheck": "tsc --noEmit"
   }
   ```

3. **MCP 导出**：
   ```json
   "exports": {
     "./mcp-server": "./src/mcp-server/index.ts"
   }
   ```

4. **新增依赖**：
   ```json
   "dependencies": {
     "fastmcp": "^1.0.0",
     "nanoid": "^5.1.5",
     "zod": "^4.0.0"
   },
   "devDependencies": {
     "fast-check": "^3.23.0"
   }
   ```

### 3. 删除的内容

- ✅ 删除 `apps/mcp-video/` 目录
- ✅ 删除 `apps/mcp-video/package.json`（依赖已合并到 packages/video）

## 重构原因

### 1. 功能内聚性
- MCP Server 是视频功能的**对外接口层**
- `packages/video` 是视频功能的**核心实现层**
- 它们是同一个功能模块的不同层次，应该在一起

### 2. 依赖关系清晰
```typescript
// MCP Server 依赖 video 核心包
import { VideoProject, templates } from '@creator-flow/video';
```
这是典型的"接口层依赖实现层"关系，应该在同一个包内。

### 3. 部署单元统一
- MCP Server 不是独立的应用（app）
- 它是视频包的一个**入口点**（entry point）
- 类似于 `packages/video/src/renderer/` 和 `packages/video/src/skills/`

### 4. 架构清晰
现在 `packages/video` 包含所有视频相关代码：
- **核心层**：compositions, templates, components
- **渲染层**：renderer
- **接口层**：skills (AI Agent), mcp-server (MCP Protocol)

## 使用方式

### 1. 启动 MCP Server

**Stdio 模式（用于 Claude Desktop、Kiro）：**
```bash
cd packages/video
bun run mcp:start
```

**HTTP 模式（用于远程部署）：**
```bash
cd packages/video
bun run mcp:start:http --port 3000
```

**开发模式（热重载）：**
```bash
cd packages/video
bun run mcp:dev
```

### 2. 在其他代码中导入

**导入 MCP Server：**
```typescript
import { runServer, createServer } from '@creator-flow/video/mcp-server';
```

**导入视频核心功能：**
```typescript
import { VideoProject, templates } from '@creator-flow/video';
import { TitleAnimation } from '@creator-flow/video/compositions';
import { renderVideo } from '@creator-flow/video/renderer';
```

### 3. Claude Desktop 配置

```json
{
  "mcpServers": {
    "creator-flow-video": {
      "command": "bun",
      "args": ["run", "/path/to/packages/video/src/mcp-server/index.ts"]
    }
  }
}
```

## 后续工作

### 1. 更新文档引用
需要更新以下文档中的路径引用：
- `docs/specs/mcp-video-server/design.md`
- `docs/specs/bun-video-integration/design.md`
- `docs/AI视频创作/remotion-integration-plan-v2.md`

### 2. 更新 Electron 集成代码
如果 `apps/electron` 中有引用 `apps/mcp-video` 的代码，需要更新为：
```typescript
// 旧的
import { ... } from '@creator-flow/mcp-video';

// 新的
import { ... } from '@creator-flow/video/mcp-server';
```

### 3. 更新 Bun 启动脚本
在 `apps/electron/src/main/video/service-manager.ts` 中启动 MCP Server：
```typescript
const mcpServerEntry = path.join(
  __dirname,
  '../../../packages/video/src/mcp-server/index.ts'
);
```

### 4. 安装依赖
```bash
# 在项目根目录
bun install
```

## 验证清单

- [x] MCP Server 代码已移动到 `packages/video/src/mcp-server/`
- [x] `apps/mcp-video/` 目录已删除
- [x] `packages/video/package.json` 已更新
- [x] 目录结构验证完成
- [ ] 依赖安装成功（需要解决证书问题）
- [ ] MCP Server 启动测试
- [ ] 文档路径更新
- [ ] Electron 集成代码更新

## 注意事项

1. **依赖安装问题**：当前遇到 `CERT_HAS_EXPIRED` 错误，需要：
   - 更新系统证书
   - 或使用 `bun install --ignore-scripts`
   - 或配置 npm registry 镜像

2. **导入路径**：所有从 `@creator-flow/mcp-video` 的导入需要改为 `@creator-flow/video/mcp-server`

3. **TypeScript 配置**：`packages/video/src/mcp-server/tsconfig.json` 保留独立配置，用于 MCP Server 的类型检查

## 架构优势

### 之前的问题
```
apps/mcp-video/          # 看起来像独立应用
  └── 依赖 @creator-flow/video  # 但实际上依赖 video 包
```
这种结构让人困惑：MCP Server 到底是独立应用还是 video 的一部分？

### 现在的架构
```
packages/video/
  ├── src/
  │   ├── compositions/      # 核心：视频组件
  │   ├── renderer/          # 核心：渲染引擎
  │   ├── skills/            # 接口：AI Agent 工具
  │   └── mcp-server/        # 接口：MCP 协议服务器
```
清晰表达：MCP Server 是 video 包的一个**接口层**，与 skills 同级。

## 总结

这次重构将 `apps/mcp-video` 移动到 `packages/video/src/mcp-server/`，使得：

1. ✅ **架构更清晰**：视频功能统一在一个包内
2. ✅ **依赖更合理**：接口层和实现层在同一个包
3. ✅ **维护更简单**：相关代码集中管理
4. ✅ **语义更准确**：MCP Server 是入口点，不是独立应用

重构完成后，`packages/video` 成为一个完整的视频创作包，包含：
- 核心实现（compositions, renderer, templates）
- 对外接口（skills, mcp-server）
- 工具函数（utils, hooks）
