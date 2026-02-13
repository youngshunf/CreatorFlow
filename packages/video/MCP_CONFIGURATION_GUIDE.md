# 视频 MCP 服务器配置指南

## 概述

本指南说明如何在 Sprouty AI 中配置视频 MCP 服务器，使 Agent 能够使用视频创作工具。

## 配置步骤

### 1. 创建数据源目录

数据源存储在工作区的 `sources` 目录下：

```bash
# 假设工作区 ID 为 800d2a01-8dd5-b485-533b-0efbebc52bb0
mkdir -p ~/.creator-flow/workspaces/800d2a01-8dd5-b485-533b-0efbebc52bb0/sources/video-mcp
```

### 2. 创建 config.json

在数据源目录下创建 `config.json` 文件：

**文件路径**: `~/.creator-flow/workspaces/{workspaceId}/sources/video-mcp/config.json`

```json
{
  "id": "video-mcp-001",
  "name": "视频创作服务",
  "slug": "video-mcp",
  "enabled": true,
  "provider": "sprouty-ai-video",
  "type": "mcp",
  "mcp": {
    "transport": "stdio",
    "command": "bun",
    "args": [
      "run",
      "/Users/mac/saas/creator-flow/.zcf/Sprouty AI/video-integration/packages/video/src/mcp-server/index.ts"
    ],
    "env": {
      "NODE_ENV": "development"
    }
  },
  "icon": "🎬",
  "tagline": "视频项目管理、素材处理、Remotion 组合编辑和视频渲染",
  "isAuthenticated": true,
  "connectionStatus": "connected",
  "createdAt": 1738684800000,
  "updatedAt": 1738684800000
}
```

### 3. 创建 guide.md

在同一目录下创建 `guide.md` 文件，提供使用指南：

**文件路径**: `~/.creator-flow/workspaces/{workspaceId}/sources/video-mcp/guide.md`

```markdown
# 视频创作服务使用指南

视频创作服务提供完整的视频制作工作流，从项目管理到最终渲染。

## 功能范围

### 项目管理
- 创建、列出、获取、更新、删除视频项目
- 支持多个项目并行工作
- 项目存储在 `{workspace}/视频创作/{项目名称}/` 目录

### 素材管理
- 添加图片、视频、音频、字体素材
- 自动发现工作区中的可用素材
- 素材存储在项目的 `素材/` 目录下

### 组合编辑
- 使用 Remotion 创建视频组合
- 支持 React + TypeScript 代码
- 代码验证和语法检查
- 组合存储在项目的 `组合/` 目录

### 视频渲染
- 渲染视频到 MP4 文件
- 实时跟踪渲染进度
- 支持多种质量和格式选项
- 输出保存在项目的 `输出/` 目录

### 预览功能
- 启动本地预览服务器
- 实时查看视频效果
- 支持热重载

## 可用工具

### 项目管理 (5 个工具)
- `video_create_project` - 创建新的视频项目
- `video_list_projects` - 列出所有视频项目
- `video_get_project` - 获取项目详情
- `video_update_project` - 更新项目配置
- `video_delete_project` - 删除项目

### 素材管理 (4 个工具)
- `video_add_asset` - 添加素材到项目
- `video_remove_asset` - 移除素材
- `video_list_assets` - 列出项目素材
- `video_list_available_assets` - 发现工作区中的可用素材

### 组合管理 (4 个工具)
- `video_add_composition` - 添加 Remotion 组合
- `video_update_composition` - 更新组合代码
- `video_remove_composition` - 移除组合
- `video_validate_composition` - 验证组合代码

### 渲染与预览 (4 个工具)
- `video_render` - 渲染视频到文件
- `video_get_render_status` - 查看渲染进度
- `video_preview_start` - 启动预览服务器
- `video_preview_stop` - 停止预览服务器

### 模板管理 (2 个工具)
- `video_list_templates` - 列出可用模板
- `video_get_template` - 获取模板详情

## 使用指南

### 创建视频项目

```typescript
// 1. 列出可用模板
const templates = await video_list_templates();

// 2. 创建项目
const project = await video_create_project({
  workspacePath: "/path/to/workspace",
  name: "我的视频项目",
  template: "social-media-vertical"
});
```

### 添加素材

```typescript
// 1. 发现可用素材
const assets = await video_list_available_assets({
  workspacePath: "/path/to/workspace",
  assetType: "image"
});

// 2. 添加素材到项目
await video_add_asset({
  workspacePath: "/path/to/workspace",
  projectName: "我的视频项目",
  assetPath: assets[0].path,
  assetType: "image"
});
```

### 创建视频组合

```typescript
// 1. 验证代码
const validation = await video_validate_composition({
  workspacePath: "/path/to/workspace",
  code: remotionCode
});

// 2. 如果验证通过，添加组合
if (validation.valid) {
  await video_add_composition({
    workspacePath: "/path/to/workspace",
    projectName: "我的视频项目",
    compositionId: "MainVideo",
    code: remotionCode,
    props: {
      title: "视频标题",
      duration: 30
    }
  });
}
```

### 预览和渲染

```typescript
// 1. 启动预览
const preview = await video_preview_start({
  workspacePath: "/path/to/workspace",
  projectName: "我的视频项目"
});
// 在浏览器中打开 preview.url

// 2. 渲染视频
const render = await video_render({
  workspacePath: "/path/to/workspace",
  projectName: "我的视频项目",
  compositionId: "MainVideo",
  outputFormat: "mp4",
  quality: "high"
});

// 3. 跟踪渲染进度
const status = await video_get_render_status({
  workspacePath: "/path/to/workspace",
  projectName: "我的视频项目"
});
```

## 最佳实践

1. **始终先预览再渲染** - 渲染需要时间，预览可以快速验证效果
2. **使用代码验证** - 在添加组合前验证代码，避免渲染失败
3. **发现素材** - 使用 `video_list_available_assets` 自动发现工作区中的素材
4. **使用模板** - 从模板开始可以节省时间，避免从零编写代码
5. **检查渲染状态** - 使用 `video_get_render_status` 跟踪长时间渲染任务

## 注意事项

- 项目名称会作为目录名，避免使用特殊字符
- 素材文件会被复制到项目目录，不会修改原文件
- 渲染过程中不要关闭应用
- 预览服务器会占用端口，同时只能预览一个项目
- 代码验证只做基本检查，复杂错误可能在渲染时才发现

## 技术细节

- **框架**: Remotion (React-based video framework)
- **运行时**: Bun
- **存储**: 文件系统 (JSON + TypeScript)
- **渲染**: FFmpeg (通过 Remotion)
- **预览**: Vite + React

## 故障排除

### 问题：MCP 服务器无法启动
- 检查 Bun 是否已安装：`bun --version`
- 检查路径是否正确
- 查看日志：`~/.creator-flow/logs/`

### 问题：渲染失败
- 使用 `video_validate_composition` 验证代码
- 检查素材文件是否存在
- 查看渲染错误信息

### 问题：预览无法打开
- 检查端口是否被占用
- 确认防火墙设置
- 尝试重启预览服务器
```

### 4. 配置权限（可选）

如果需要限制工具使用，可以创建权限配置：

**文件路径**: `~/.creator-flow/workspaces/{workspaceId}/sources/video-mcp/permissions.json`

```json
{
  "allowedMcpPatterns": [
    "^video_.*"
  ],
  "blockedTools": []
}
```

## 快速配置脚本

创建一个脚本来自动配置：

```bash
#!/bin/bash

# 配置变量
WORKSPACE_ID="800d2a01-8dd5-b485-533b-0efbebc52bb0"
SOURCE_DIR="$HOME/.creator-flow/workspaces/$WORKSPACE_ID/sources/video-mcp"
VIDEO_SERVER_PATH="/Users/mac/saas/creator-flow/.zcf/Sprouty AI/video-integration/packages/video/src/mcp-server/index.ts"

# 创建目录
mkdir -p "$SOURCE_DIR"

# 创建 config.json
cat > "$SOURCE_DIR/config.json" << 'EOF'
{
  "id": "video-mcp-001",
  "name": "视频创作服务",
  "slug": "video-mcp",
  "enabled": true,
  "provider": "sprouty-ai-video",
  "type": "mcp",
  "mcp": {
    "transport": "stdio",
    "command": "bun",
    "args": [
      "run",
      "VIDEO_SERVER_PATH_PLACEHOLDER"
    ],
    "env": {
      "NODE_ENV": "development"
    }
  },
  "icon": "🎬",
  "tagline": "视频项目管理、素材处理、Remotion 组合编辑和视频渲染",
  "isAuthenticated": true,
  "connectionStatus": "connected",
  "createdAt": $(date +%s)000,
  "updatedAt": $(date +%s)000
}
EOF

# 替换路径
sed -i '' "s|VIDEO_SERVER_PATH_PLACEHOLDER|$VIDEO_SERVER_PATH|g" "$SOURCE_DIR/config.json"

echo "✅ 视频 MCP 服务器配置已创建"
echo "📁 配置目录: $SOURCE_DIR"
echo ""
echo "下一步："
echo "1. 重启 Sprouty AI 应用"
echo "2. 在数据源列表中查看 '视频创作服务'"
echo "3. 测试连接"
```

## 验证配置

配置完成后，在 Sprouty AI 中：

1. 打开工作区设置
2. 进入"数据源"页面
3. 应该看到"视频创作服务"数据源
4. 点击"测试连接"验证配置
5. 如果连接成功，状态应显示为"已连接"

## 在 Agent 会话中使用

配置完成后，Agent 可以直接使用视频工具：

```
用户: 帮我创建一个 30 秒的抖音视频，介绍我们的 AI 产品

Agent: 我将使用视频创作服务来帮你创建视频...
[调用 video_list_templates]
[调用 video_create_project]
[调用 video_add_composition]
[调用 video_preview_start]
...
```

## 注意事项

1. **路径配置**: 确保 `command` 和 `args` 中的路径正确
2. **Bun 运行时**: 确保系统已安装 Bun
3. **工作区 ID**: 替换为实际的工作区 ID
4. **重启应用**: 配置后需要重启 Sprouty AI 才能生效
5. **日志查看**: 如有问题，查看 `~/.creator-flow/logs/` 目录下的日志
