# 视频 MCP 服务器内置配置

## 概述

视频 MCP 服务器现已作为内置数据源集成到"自媒体创作"应用中。当用户创建新的工作区时，视频 MCP 服务器会自动配置并可用。

## 实现方案

### 1. 应用级配置

视频 MCP 数据源配置位于：
```
packages/shared/src/apps/bundled-apps/app-creator-media/sources/video-mcp/
├── config.json  # 数据源配置（包含路径占位符）
└── guide.md     # 使用指南
```

### 2. 路径占位符系统

为了解决开发环境和打包后应用的路径差异，我们使用占位符系统：

**config.json 中的占位符**:
```json
{
  "mcp": {
    "command": "{{BUN_PATH}}",
    "args": [
      "run",
      "{{VIDEO_MCP_SERVER_PATH}}"
    ]
  }
}
```

**占位符说明**:
- `{{BUN_PATH}}` - Bun 可执行文件路径
- `{{VIDEO_MCP_SERVER_PATH}}` - MCP 服务器入口文件路径

### 3. 路径解析模块

路径解析逻辑位于 `packages/shared/src/apps/video-mcp-paths.ts`：

```typescript
// 获取 Bun 路径
getBunPath(isPackaged: boolean): string

// 获取 MCP 服务器路径
getVideoMcpServerPath(isPackaged: boolean, resourcesPath?: string): string

// 解析配置中的占位符
resolveSourceConfigPaths(sourceConfig: any, isPackaged: boolean, resourcesPath?: string): any
```

**路径解析策略**:

**开发环境** (`isPackaged = false`):
- Bun: 使用系统 PATH 中的 `bun`
- MCP Server: 使用相对路径推导
  ```
  packages/shared/src/apps/video-mcp-paths.ts
  -> packages/video/src/mcp-server/index.ts
  ```

**打包后** (`isPackaged = true`):
- Bun: 使用系统 PATH 中的 `bun`
- MCP Server: 尝试多个可能的路径
  ```
  1. {resourcesPath}/app/packages/video/src/mcp-server/index.ts
  2. {resourcesPath}/app.asar/packages/video/src/mcp-server/index.ts
  3. {resourcesPath}/packages/video/src/mcp-server/index.ts
  ```

### 4. 工作区初始化流程

当用户创建新工作区时：

1. **调用初始化函数**
   ```typescript
   initializeWorkspaceFromApp({
     name: "我的工作区",
     appId: "app.creator-media",
     isPackaged: app.isPackaged,
     resourcesPath: app.getPath('resources')
   })
   ```

2. **复制数据源配置**
   - 从 `app-creator-media/sources/video-mcp/` 复制到工作区
   - 目标路径: `{workspace}/.creator-flow/sources/video-mcp/`

3. **解析路径占位符**
   - 读取 `config.json`
   - 替换 `{{BUN_PATH}}` 和 `{{VIDEO_MCP_SERVER_PATH}}`
   - 写入解析后的配置

4. **最终配置示例**
   ```json
   {
     "mcp": {
       "command": "bun",
       "args": [
         "run",
         "/path/to/packages/video/src/mcp-server/index.ts"
       ]
     }
   }
   ```

## 使用方式

### 对于用户

1. **创建工作区**
   - 打开 Sprouty AI
   - 选择"自媒体创作"应用
   - 创建新工作区

2. **自动配置**
   - 视频 MCP 服务器自动配置
   - 无需手动设置

3. **开始使用**
   - 在 Agent 会话中直接请求视频创作任务
   - 例如："帮我创建一个 30 秒的抖音视频"

### 对于开发者

#### 在 Electron 主进程中调用

```typescript
import { initializeWorkspaceFromApp } from '@sprouty-ai/shared/apps';
import { app } from 'electron';

// 创建工作区时
const result = initializeWorkspaceFromApp({
  name: workspaceName,
  appId: 'app.creator-media',
  rootPath: workspaceRootPath,
  isPackaged: app.isPackaged,
  resourcesPath: app.getPath('resources')
});

if (result.success) {
  console.log('工作区初始化成功');
  console.log('视频 MCP 服务器已自动配置');
}
```

#### 添加新的内置 MCP 数据源

1. **创建数据源目录**
   ```bash
   mkdir -p packages/shared/src/apps/bundled-apps/app-xxx/sources/my-mcp
   ```

2. **创建 config.json**
   ```json
   {
     "id": "my-mcp_preset",
     "name": "我的 MCP 服务",
     "slug": "my-mcp",
     "enabled": true,
     "type": "mcp",
     "mcp": {
       "transport": "stdio",
       "command": "{{MY_COMMAND}}",
       "args": ["{{MY_ARG}}"]
     }
   }
   ```

3. **创建 guide.md**
   - 提供使用指南和工具说明

4. **扩展路径解析模块**
   - 在 `video-mcp-paths.ts` 中添加新的占位符解析逻辑
   - 或创建新的路径解析模块

## 打包配置

确保 MCP 服务器文件被正确打包：

**electron-builder.yml**:
```yaml
files:
  - packages/video/src/**/*
  - packages/video/package.json
  - packages/video/remotion.config.ts
```

**注意事项**:
- TypeScript 文件会被打包（Bun 可以直接运行 .ts 文件）
- 确保 `node_modules` 中的依赖也被打包
- 测试打包后的应用，验证路径解析正确

## 测试

### 开发环境测试

```bash
# 启动开发服务
bun run electron:dev

# 创建新工作区
# 检查 ~/.creator-flow/workspaces/{id}/sources/video-mcp/config.json
# 验证路径是否正确解析
```

### 打包后测试

```bash
# 打包应用
bun run electron:dist:mac

# 安装并运行打包后的应用
# 创建新工作区
# 检查数据源配置
# 测试视频工具是否可用
```

## 故障排除

### 问题：MCP 服务器无法启动

**检查步骤**:
1. 查看配置文件
   ```bash
   cat ~/.creator-flow/workspaces/{id}/sources/video-mcp/config.json
   ```

2. 验证路径
   ```bash
   # 检查 Bun
   which bun

   # 检查 MCP 服务器文件
   ls -la {path-from-config}
   ```

3. 查看日志
   ```bash
   tail -f ~/.creator-flow/logs/main.log
   ```

### 问题：路径解析错误

**可能原因**:
- `isPackaged` 参数未正确传递
- `resourcesPath` 为空或错误
- 打包配置遗漏文件

**解决方案**:
1. 检查 `initializeWorkspaceFromApp` 调用
2. 验证 `app.isPackaged` 和 `app.getPath('resources')`
3. 检查 electron-builder 配置

### 问题：占位符未被替换

**可能原因**:
- `resolveSourceConfigPaths` 未被调用
- 配置文件格式错误

**解决方案**:
1. 检查 `copyAppSourcesToWorkspace` 实现
2. 验证 JSON 格式
3. 添加调试日志

## 优势

1. **用户友好**: 无需手动配置，开箱即用
2. **环境适配**: 自动处理开发和生产环境差异
3. **可维护性**: 集中管理路径解析逻辑
4. **可扩展性**: 易于添加新的内置 MCP 数据源
5. **版本控制**: 配置随应用版本更新

## 未来改进

1. **动态路径发现**: 使用更智能的路径查找算法
2. **配置验证**: 在初始化时验证 MCP 服务器可用性
3. **错误恢复**: 自动修复损坏的配置
4. **多 MCP 支持**: 支持配置多个 MCP 服务器
5. **版本兼容**: 处理 MCP 服务器版本升级
