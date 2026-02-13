# Video MCP Source 使用指南

## 概述

Video MCP 是基于 Remotion 的视频创作引擎 MCP 服务，提供完整的视频项目管理、模板系统、素材管理和渲染能力。

## 配置说明

### 开发环境配置

在开发环境中，使用相对路径指向 monorepo 中的 video 包：

```json
{
  "mcp": {
    "command": "bun",
    "args": ["run", "src/mcp-server/index.ts"],
    "cwd": "../../packages/video",
    "transport": "stdio"
  }
}
```

### 生产环境配置（打包后）

在打包后的应用中，使用 `app:` 前缀指向应用资源目录：

```json
{
  "mcp": {
    "command": "bun",
    "args": ["run", "src/mcp-server/index.ts"],
    "cwd": "app:resources/video",
    "transport": "stdio"
  }
}
```

### cwd 路径解析规则

- `app:` 前缀：相对于应用安装目录（`process.resourcesPath`）
  - 示例：`app:resources/video` → `/Applications/Sprouty AI.app/Contents/Resources/resources/video`
- 绝对路径：直接使用
  - 示例：`/usr/local/lib/video-mcp`
- 相对路径：相对于工作区根目录
  - 示例：`packages/video` → `~/.sprouty-ai/workspaces/{id}/packages/video`

## 可用工具（18 个）

### 项目管理（5 个）
- `video_create_project` - 创建新的视频项目
- `video_list_projects` - 列出所有视频项目
- `video_get_project` - 获取项目详情
- `video_update_project` - 更新项目配置
- `video_delete_project` - 删除视频项目

### 素材管理（4 个）
- `video_add_asset` - 添加素材到项目
- `video_list_assets` - 列出项目素材
- `video_get_asset` - 获取素材详情
- `video_remove_asset` - 从项目移除素材

### 组合管理（3 个）
- `video_create_composition` - 创建视频组合
- `video_list_compositions` - 列出所有组合
- `video_get_composition` - 获取组合详情

### 模板系统（2 个）
- `video_list_templates` - 列出可用模板
- `video_get_template` - 获取模板详情

### 渲染与预览（3 个）
- `video_render` - 渲染视频
- `video_preview` - 生成预览
- `video_get_render_status` - 查询渲染状态

### 代码校验（1 个）
- `video_validate_code` - 校验 Remotion 代码

## 支持的模板类型

1. **Slideshow** - 幻灯片展示
   - 支持多张图片轮播
   - 可自定义标题、描述、过渡效果

2. **TitleAnimation** - 标题动画
   - 支持主标题和副标题
   - 多种动画效果（淡入、滑入、缩放等）

3. **DataVisualization** - 数据可视化
   - 柱状图（Bar Chart）
   - 折线图（Line Chart）
   - 饼图（Pie Chart）
   - 支持动态数据和动画效果

4. **ProductShowcase** - 产品展示
   - 产品特性列表
   - 产品图片展示
   - 多种样式变体（Default、Minimal、Premium）

## 使用示例

### 创建项目

```typescript
await video_create_project({
  name: "我的第一个视频",
  description: "使用 Remotion 创建的视频项目",
  fps: 30,
  width: 1920,
  height: 1080
});
```

### 创建组合

```typescript
await video_create_composition({
  projectId: "project-123",
  template: "Slideshow",
  props: {
    title: "产品介绍",
    items: [
      { image: "slide1.jpg", title: "特性 1", description: "描述..." },
      { image: "slide2.jpg", title: "特性 2", description: "描述..." }
    ],
    logo: "logo.png"
  }
});
```

### 渲染视频

```typescript
await video_render({
  compositionId: "comp-456",
  outputPath: "/path/to/output.mp4",
  codec: "h264",
  quality: "high"
});
```

## 故障排除

### 问题：bun 命令找不到

**原因**：Electron 子进程的 PATH 环境变量不包含 bun 的安装路径。

**解决方案**：系统会自动探测 bun 路径，支持以下位置：
- `~/.bun/bin/bun`
- `/opt/homebrew/bin/bun`
- `/usr/local/bin/bun`
- `~/.deno/bin/deno`（如果使用 deno）

如果自动探测失败，可以在配置中使用绝对路径：

```json
{
  "mcp": {
    "command": "/Users/username/.bun/bin/bun",
    "args": ["run", "src/mcp-server/index.ts"]
  }
}
```

### 问题：打包后找不到 video 包

**原因**：打包时没有正确复制 video 包到应用资源目录。

**解决方案**：
1. 确认 `apps/electron/scripts/copy-assets.ts` 包含 video 包复制逻辑
2. 使用 `app:resources/video` 作为 cwd 路径
3. 重新构建应用

### 问题：渲染失败

**常见原因**：
- 素材路径不正确（使用 `staticFile()` 包装路径）
- Props 类型不匹配（参考模板文档）
- 视频时长计算错误（使用 `calculateMetadata`）

**调试步骤**：
1. 使用 `video_validate_code` 校验代码
2. 检查 MCP 服务器日志
3. 在 Remotion Studio 中预览组合

## 最佳实践

1. **资源管理**：所有图片路径使用 `staticFile()` 包装
2. **类型安全**：使用 TypeScript 定义 Props 类型
3. **动态时长**：使用 `calculateMetadata` 根据内容计算视频时长
4. **组织结构**：使用 `<Folder>` 在 Remotion Studio 中分组组合
5. **性能优化**：避免在渲染循环中进行重计算，使用 `useMemo`

## 相关文档

- [Remotion 官方文档](https://www.remotion.dev/docs)
- [Video 包 README](../../packages/video/README.md)
- [优化总结](../../packages/video/OPTIMIZATION_SUMMARY.md)
- [变更日志](../../packages/video/CHANGELOG.md)
