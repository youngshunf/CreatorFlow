# Sprouty AI 视频创作功能使用指南

本文档介绍 Sprouty AI 桌面应用中集成的 Remotion 视频创作功能。

## 功能概述

Sprouty AI 集成了 Remotion 视频创作引擎，提供以下核心能力：

- **程序化视频创作** - 使用 React 组件定义视频内容和动画
- **预设模板系统** - 快速创建社交媒体、营销、教程等类型视频
- **实时预览** - 在编辑过程中即时预览视频效果
- **多格式导出** - 支持 MP4、WebM、GIF 格式输出
- **AI Agent 集成** - 通过自然语言描述生成视频内容

## 快速开始

### 访问视频编辑器

在 Sprouty AI 应用中，通过以下方式访问视频创作功能：

```typescript
import { routes } from '@/shared/routes';
import { useNavigate } from '@/contexts/NavigationContext';

// 打开视频编辑器
const { navigate } = useNavigate();
navigate(routes.view.video());

// 打开特定项目
navigate(routes.view.video('project-id'));
```

### 创建视频项目

```typescript
// 通过 Electron API 创建项目
const project = await window.electronAPI.video.createProject({
  name: '我的视频项目',
  workspaceId: 'workspace-id',
  template: 'social-media',  // 可选模板
  config: {
    width: 1920,
    height: 1080,
    fps: 30,
    durationInFrames: 300,  // 10秒 @ 30fps
  },
});
```

## API 参考

### VideoAPI 接口

通过 `window.electronAPI.video` 访问以下方法：

#### 项目管理

| 方法 | 说明 |
|------|------|
| `createProject(options)` | 创建新视频项目 |
| `listProjects(workspaceId)` | 列出工作区内所有项目 |
| `getProject(projectId)` | 获取单个项目详情 |
| `updateProject(projectId, updates)` | 更新项目配置 |
| `deleteProject(projectId)` | 删除项目 |

#### 素材管理

| 方法 | 说明 |
|------|------|
| `addAsset(projectId, assetPath, assetType)` | 添加素材到项目 |
| `removeAsset(projectId, assetId)` | 从项目移除素材 |

**支持的素材类型：**
- `image` - 图片 (png, jpg, gif, webp, svg)
- `video` - 视频 (mp4, webm, mov)
- `audio` - 音频 (mp3, wav, ogg, m4a)
- `font` - 字体 (ttf, otf, woff, woff2)

#### 渲染

| 方法 | 说明 |
|------|------|
| `render(options)` | 渲染视频到文件 |
| `cancelRender()` | 取消当前渲染 |
| `onRenderProgress(callback)` | 监听渲染进度 |

**渲染选项：**
```typescript
interface RenderOptions {
  projectId: string;
  compositionId: string;
  outputFormat?: 'mp4' | 'webm' | 'gif';
  quality?: 'draft' | 'standard' | 'high';
}
```

#### 预览

| 方法 | 说明 |
|------|------|
| `startPreview(projectId)` | 启动预览服务器 |
| `stopPreview(projectId)` | 停止预览服务器 |

## 预设模板

### 可用模板

| 模板 ID | 名称 | 分辨率 | 适用场景 |
|---------|------|--------|----------|
| `social-media` | 社交媒体横屏 | 1920×1080 | YouTube、B站 |
| `social-media-vertical` | 社交媒体竖屏 | 1080×1920 | 抖音、小红书 |
| `social-media-square` | 社交媒体方形 | 1080×1080 | Instagram |
| `marketing` | 营销视频 | 1920×1080 | 产品宣传 |
| `tutorial` | 教程视频 | 1920×1080 | 操作演示 |

### 使用模板

```typescript
const project = await window.electronAPI.video.createProject({
  name: '产品介绍',
  workspaceId: 'ws-123',
  template: 'marketing',
});
```

## 内置 Composition 组件

### TitleAnimation
标题动画组件，支持主标题和副标题的淡入效果。

### Slideshow
图片轮播组件，支持多图片自动切换展示。

### DataVisualization
数据可视化组件，支持柱状图、折线图等动画展示。

### ProductShowcase
产品展示组件，支持产品特性的动画呈现。

## AI Agent 视频创作

Sprouty AI 的 AI Agent 可以通过自然语言指令创建视频。

### 可用工具

| 工具名称 | 功能 |
|----------|------|
| `video_create_project` | 创建视频项目 |
| `video_generate_composition` | 根据描述生成视频组件 |
| `video_update_composition` | 更新现有组件 |
| `video_preview` | 启动预览 |
| `video_render` | 渲染导出 |
| `video_add_asset` | 添加素材 |

### 示例对话

```
用户: 帮我创建一个 15 秒的产品介绍视频，竖屏格式

Agent: 好的，我来为您创建一个竖屏产品介绍视频。

[调用 video_create_project]
{
  "name": "产品介绍视频",
  "template": "social-media-vertical",
  "durationInSeconds": 15
}

项目已创建！您可以：
1. 添加产品图片素材
2. 让我生成具体的视频内容
3. 直接预览当前效果
```

详细的 Agent 使用指南请参考：`packages/video/src/skills/AGENTS.md`

## UI 组件

视频编辑器包含以下 UI 组件：

| 组件 | 功能 |
|------|------|
| `VideoEditor` | 主编辑器，三栏布局 |
| `VideoPreview` | 视频预览播放器 |
| `VideoTimeline` | 时间轴控制 |
| `VideoProperties` | 属性配置面板 |
| `VideoProjectList` | 项目列表 |
| `VideoTemplates` | 模板选择器 |
| `VideoExport` | 导出配置 |

## 文件结构

```
Sprouty AI/
├── packages/video/              # @sprouty-ai/video 包
│   └── src/
│       ├── components/          # 可复用动画组件
│       ├── compositions/        # 预设 Composition
│       ├── templates/           # 视频模板
│       ├── skills/              # AI Agent 工具
│       ├── hooks/               # React Hooks
│       ├── renderer/            # 渲染工具
│       ├── utils/               # 工具函数
│       ├── types.ts             # 类型定义
│       ├── Root.tsx             # Remotion 根组件
│       └── index.ts             # 包入口
│
└── apps/electron/
    └── src/
        ├── main/video/          # 主进程视频服务
        │   ├── project-manager.ts
        │   ├── render-worker.ts
        │   ├── preview-server.ts
        │   └── ipc-handlers.ts
        ├── preload/
        │   └── video-api.ts     # Preload API
        └── renderer/components/video/  # UI 组件
            ├── VideoEditor.tsx
            ├── VideoPreview.tsx
            ├── VideoTimeline.tsx
            └── ...
```

## 注意事项

1. **依赖安装** - 首次使用前需运行 `bun install` 安装 Remotion 相关依赖
2. **内存使用** - 视频渲染可能消耗较多内存，建议关闭不必要的应用
3. **渲染时间** - 高质量渲染需要较长时间，可先用 `draft` 质量预览
4. **GIF 限制** - GIF 格式不支持音频，且文件较大

## 相关文档

- [Remotion 官方文档](https://www.remotion.dev/docs)
- [AI Agent 视频创作指南](../packages/video/src/skills/AGENTS.md)
- [设计文档](../../.kiro/specs/remotion-video-integration/design.md)
