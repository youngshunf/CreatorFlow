# 视频创作服务

视频创作服务提供完整的视频制作工作流，从项目管理到最终渲染。基于 Remotion 框架，支持使用 React + TypeScript 创建专业视频。

## 功能范围

- **项目管理**: 创建、管理多个视频项目
- **素材管理**: 自动发现和管理图片、视频、音频、字体素材
- **组合编辑**: 使用 Remotion 创建视频组合，支持代码验证
- **视频渲染**: 渲染视频到 MP4 文件，实时跟踪进度
- **预览功能**: 本地预览服务器，实时查看效果

## 可用工具

### 项目管理 (5 个)
- `video_create_project` - 创建新项目
- `video_list_projects` - 列出所有项目
- `video_get_project` - 获取项目详情
- `video_update_project` - 更新项目
- `video_delete_project` - 删除项目

### 素材管理 (4 个)
- `video_add_asset` - 添加素材
- `video_remove_asset` - 移除素材
- `video_list_assets` - 列出素材
- `video_list_available_assets` - 发现可用素材

### 组合管理 (4 个)
- `video_add_composition` - 添加组合
- `video_update_composition` - 更新组合
- `video_remove_composition` - 移除组合
- `video_validate_composition` - 验证代码

### 渲染预览 (4 个)
- `video_render` - 渲染视频
- `video_get_render_status` - 查看进度
- `video_preview_start` - 启动预览
- `video_preview_stop` - 停止预览

### 模板管理 (2 个)
- `video_list_templates` - 列出模板
- `video_get_template` - 获取模板

## 使用指南

### 典型工作流

1. 列出可用模板
2. 创建视频项目
3. 发现工作区中的素材
4. 添加素材到项目
5. 验证 Remotion 代码
6. 添加视频组合
7. 启动预览查看效果
8. 渲染最终视频
9. 跟踪渲染进度

### 最佳实践

- 始终先预览再渲染（渲染耗时较长）
- 使用代码验证避免渲染失败
- 从模板开始节省时间
- 使用素材发现工具自动查找资源

## 技术细节

- **框架**: Remotion (React-based)
- **运行时**: Bun
- **渲染**: FFmpeg
- **存储**: `{workspace}/视频创作/{项目名称}/`
