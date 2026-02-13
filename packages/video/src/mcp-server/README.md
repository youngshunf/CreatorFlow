# Sprouty AI Video MCP Server

MCP (Model Context Protocol) 服务器，为 AI Agent 提供视频创作能力。基于 Remotion 实现视频渲染，支持项目管理、素材管理、组合编辑、视频渲染和实时预览。

## 特性

- 🎬 **项目管理** - 创建、查询、更新、删除视频项目
- 🎨 **素材管理** - 添加、移除、列出项目素材（视频、音频、图片、字体）
- 🎥 **组合管理** - 管理视频组合（场景配置）
- 📹 **视频渲染** - 基于 Remotion 的高质量视频渲染
- 👁️ **实时预览** - 启动 Remotion Studio 进行实时预览
- 📋 **模板系统** - 内置视频模板，快速创建专业视频

## 安装

```bash
cd Sprouty AI/apps/mcp-video
bun install
```

## 使用方式

### Stdio 模式（推荐用于桌面客户端）

```bash
# 开发模式（热重载）
bun run dev

# 生产模式
bun run start
```

### HTTP 模式（用于远程/云端部署）

```bash
# 默认端口 3000
bun run start:http

# 自定义端口
bun run start -- --transport http --port 8080

# 自定义主机和端口
bun run start -- --transport http --host 127.0.0.1 --port 8080
```

## 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--transport <mode>` | 传输模式：`stdio` 或 `http` | `stdio` |
| `--host <host>` | HTTP 模式的主机地址 | `0.0.0.0` |
| `--port <port>` | HTTP 模式的端口 | `3000` |
| `--endpoint <path>` | HTTP 模式的端点路径 | `/mcp` |
| `--help` | 显示帮助信息 | - |

## 客户端配置

### Claude Desktop

在 `claude_desktop_config.json` 中添加：

```json
{
  "mcpServers": {
    "sprouty-ai-video": {
      "command": "bun",
      "args": ["run", "/path/to/Sprouty AI/apps/mcp-video/src/index.ts"]
    }
  }
}
```

### Kiro

在 `.kiro/settings/mcp.json` 中添加：

```json
{
  "mcpServers": {
    "sprouty-ai-video": {
      "command": "bun",
      "args": ["run", "/path/to/Sprouty AI/apps/mcp-video/src/index.ts"],
      "disabled": false,
      "autoApprove": []
    }
  }
}
```

## 可用工具

### 📁 项目管理

| 工具名 | 说明 |
|--------|------|
| `video_create_project` | 创建新的视频项目 |
| `video_list_projects` | 列出工作区中的所有视频项目 |
| `video_get_project` | 获取项目详细信息 |
| `video_update_project` | 更新项目元数据 |
| `video_delete_project` | 删除视频项目 |

### 🎨 素材管理

| 工具名 | 说明 |
|--------|------|
| `video_add_asset` | 添加素材到项目 |
| `video_remove_asset` | 从项目移除素材 |
| `video_list_assets` | 列出项目中的所有素材 |

### 🎬 组合管理

| 工具名 | 说明 |
|--------|------|
| `video_add_composition` | 添加新的组合到项目 |
| `video_update_composition` | 更新组合配置 |
| `video_remove_composition` | 从项目移除组合 |

### 🎥 视频渲染

| 工具名 | 说明 |
|--------|------|
| `video_render` | 渲染视频项目 |

### 👁️ 实时预览

| 工具名 | 说明 |
|--------|------|
| `video_preview_start` | 启动预览服务器 |
| `video_preview_stop` | 停止预览服务器 |

### 📋 模板管理

| 工具名 | 说明 |
|--------|------|
| `video_list_templates` | 列出可用的视频模板 |

## 工具使用示例

### 创建项目

```json
{
  "tool": "video_create_project",
  "arguments": {
    "workspacePath": "/path/to/workspace",
    "name": "我的视频项目",
    "description": "一个示例视频项目",
    "templateId": "text-reveal"
  }
}
```

### 添加素材

```json
{
  "tool": "video_add_asset",
  "arguments": {
    "workspacePath": "/path/to/workspace",
    "projectName": "我的视频项目",
    "sourcePath": "/path/to/video.mp4",
    "type": "video"
  }
}
```

### 渲染视频

```json
{
  "tool": "video_render",
  "arguments": {
    "workspacePath": "/path/to/workspace",
    "projectName": "我的视频项目",
    "compositionId": "main",
    "quality": "high"
  }
}
```

### 启动预览

```json
{
  "tool": "video_preview_start",
  "arguments": {
    "workspacePath": "/path/to/workspace",
    "projectName": "我的视频项目"
  }
}
```

## 项目目录结构

创建项目后，会在工作区生成以下目录结构：

```
{workspacePath}/
└── 视频创作/
    └── {projectName}/
        ├── project.json    # 项目配置文件
        ├── 素材/           # 素材目录
        │   ├── 视频/
        │   ├── 音频/
        │   ├── 图片/
        │   └── 字体/
        ├── 组合/           # 组合配置目录
        └── 输出/           # 渲染输出目录
```

## 开发

### 运行测试

```bash
bun test
```

### 监视模式测试

```bash
bun test:watch
```

### 类型检查

```bash
bun run typecheck
```

### 构建

```bash
bun run build
```

## 依赖

- [Bun](https://bun.sh/) - JavaScript 运行时
- [FastMCP](https://github.com/jlowin/fastmcp) - MCP 服务器框架
- [Remotion](https://remotion.dev/) - 视频渲染引擎
- [@sprouty-ai/video](../packages/video) - 视频模板和组件

## 许可证

MIT
