#!/bin/bash

# 视频 MCP 服务器自动配置脚本
# 用于在 CreatorFlow 中配置视频创作服务数据源

set -e

echo "🎬 视频 MCP 服务器配置脚本"
echo "================================"
echo ""

# 获取工作区 ID
echo "请输入工作区 ID（留空使用默认值）:"
read -r WORKSPACE_ID

if [ -z "$WORKSPACE_ID" ]; then
  # 尝试从配置文件中读取第一个工作区
  CONFIG_FILE="$HOME/.sprouty-ai/config.json"
  if [ -f "$CONFIG_FILE" ]; then
    WORKSPACE_ID=$(cat "$CONFIG_FILE" | grep -o '"id":"[^"]*"' | head -1 | cut -d'"' -f4)
    echo "使用配置文件中的工作区: $WORKSPACE_ID"
  else
    echo "❌ 错误: 未找到配置文件，请手动输入工作区 ID"
    exit 1
  fi
fi

# 配置变量
SOURCE_DIR="$HOME/.sprouty-ai/workspaces/$WORKSPACE_ID/sources/video-mcp"
VIDEO_SERVER_PATH="$(pwd)/packages/video/src/mcp-server/index.ts"

# 检查视频服务器路径是否存在
if [ ! -f "$VIDEO_SERVER_PATH" ]; then
  echo "❌ 错误: 视频服务器文件不存在: $VIDEO_SERVER_PATH"
  echo "请确保在项目根目录运行此脚本"
  exit 1
fi

# 检查 Bun 是否已安装
if ! command -v bun &> /dev/null; then
  echo "❌ 错误: 未找到 Bun 运行时"
  echo "请先安装 Bun: https://bun.sh"
  exit 1
fi

echo ""
echo "配置信息:"
echo "  工作区 ID: $WORKSPACE_ID"
echo "  数据源目录: $SOURCE_DIR"
echo "  服务器路径: $VIDEO_SERVER_PATH"
echo ""

# 创建目录
echo "📁 创建数据源目录..."
mkdir -p "$SOURCE_DIR"

# 创建 config.json
echo "📝 创建 config.json..."
cat > "$SOURCE_DIR/config.json" << EOF
{
  "id": "video-mcp-$(date +%s)",
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
      "$VIDEO_SERVER_PATH"
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

# 创建 guide.md
echo "📝 创建 guide.md..."
cat > "$SOURCE_DIR/guide.md" << 'EOF'
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

### 视频渲染
- 渲染视频到 MP4 文件
- 实时跟踪渲染进度
- 支持多种质量和格式选项

## 可用工具 (19 个)

### 项目管理
- video_create_project, video_list_projects, video_get_project
- video_update_project, video_delete_project

### 素材管理
- video_add_asset, video_remove_asset, video_list_assets
- video_list_available_assets

### 组合管理
- video_add_composition, video_update_composition, video_remove_composition
- video_validate_composition

### 渲染与预览
- video_render, video_get_render_status
- video_preview_start, video_preview_stop

### 模板管理
- video_list_templates, video_get_template

## 使用示例

创建视频项目并渲染：

1. 列出模板: video_list_templates()
2. 创建项目: video_create_project(name, template)
3. 发现素材: video_list_available_assets()
4. 添加素材: video_add_asset()
5. 添加组合: video_add_composition(code, props)
6. 预览: video_preview_start()
7. 渲染: video_render()

## 最佳实践

- 始终先预览再渲染
- 使用代码验证避免错误
- 从模板开始节省时间
- 跟踪渲染进度
EOF

echo ""
echo "✅ 配置完成！"
echo ""
echo "📋 配置摘要:"
echo "  - 数据源名称: 视频创作服务"
echo "  - 数据源 slug: video-mcp"
echo "  - 工具数量: 19 个"
echo "  - 配置文件: $SOURCE_DIR/config.json"
echo "  - 使用指南: $SOURCE_DIR/guide.md"
echo ""
echo "🚀 下一步:"
echo "  1. 重启 CreatorFlow 应用"
echo "  2. 打开工作区设置 → 数据源"
echo "  3. 查看 '视频创作服务' 数据源"
echo "  4. 点击 '测试连接' 验证配置"
echo ""
echo "💡 提示: 如果连接失败，请检查:"
echo "  - Bun 是否已安装: bun --version"
echo "  - 服务器路径是否正确"
echo "  - 查看日志: ~/.sprouty-ai/logs/"
echo ""
