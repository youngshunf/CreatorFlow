# Remotion 视频创作集成方案 v2

## 概述

本方案将 [Remotion](https://github.com/remotion-dev/remotion) 集成到 Sprouty AI 桌面端，实现 AI 驱动的程序化视频创作能力。

## 一、核心价值

### 为什么选择 Remotion？

1. **React 生态兼容** - Sprouty AI 已使用 React，无缝集成
2. **AI 友好** - Remotion 官方提供 Claude Code Agent Skills，专为 AI 生成视频优化
3. **程序化创作** - 代码即视频，完美契合 AI 生成场景
4. **灵活渲染** - 支持本地渲染、服务端渲染、客户端渲染（实验性）

### 用户场景

```
创作者: "帮我创建一个产品介绍视频，包含标题动画、3个产品特性展示、结尾 CTA"

AI Agent:
1. 分析需求，规划视频结构
2. 生成 Remotion 组件代码
3. 实时预览，支持迭代修改
4. 渲染输出 MP4
```

## 二、技术架构

### 整体架构

```
┌─────────────────────────────────────────────────────────────────┐
│                    Sprouty AI Desktop                          │
├─────────────────────────────────────────────────────────────────┤
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐ │
│  │   AI Agent      │  │  Video Editor   │  │  Render Engine  │ │
│  │  (Claude SDK)   │◄─┤    (React UI)   │◄─┤   (Remotion)    │ │
│  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘ │
│           │                    │                    │          │
│           ▼                    ▼                    ▼          │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │              @sprouty-ai/video (新包)                     ││
│  │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────────┐   ││
│  │  │Templates │ │Compositions│ │ Renderer │ │ AI Skills    │   ││
│  │  └──────────┘ └──────────┘ └──────────┘ └──────────────┘   ││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────┘
```

### 目录结构

```
Sprouty AI/
├── packages/
│   └── video/                          # 新增：视频处理核心包
│       ├── package.json
│       ├── tsconfig.json
│       ├── remotion.config.ts          # Remotion 配置
│       ├── src/
│       │   ├── index.ts                # 导出入口
│       │   ├── Root.tsx                # Remotion 根组件
│       │   ├── compositions/           # 视频组件
│       │   │   ├── index.ts
│       │   │   ├── TitleAnimation.tsx
│       │   │   ├── Slideshow.tsx
│       │   │   ├── DataVisualization.tsx
│       │   │   └── ProductShowcase.tsx
│       │   ├── templates/              # 预设模板
│       │   │   ├── index.ts
│       │   │   ├── social-media.ts     # 社交媒体模板
│       │   │   ├── marketing.ts        # 营销视频模板
│       │   │   └── tutorial.ts         # 教程视频模板
│       │   ├── components/             # 可复用组件
│       │   │   ├── AnimatedText.tsx
│       │   │   ├── Transition.tsx
│       │   │   ├── Background.tsx
│       │   │   └── Logo.tsx
│       │   ├── hooks/                  # Remotion hooks
│       │   │   ├── useAnimation.ts
│       │   │   └── useAssets.ts
│       │   ├── renderer/               # 渲染器
│       │   │   ├── index.ts
│       │   │   ├── local-renderer.ts
│       │   │   └── preview-server.ts
│       │   ├── skills/                 # AI Agent Skills
│       │   │   ├── index.ts
│       │   │   ├── video-creation.ts
│       │   │   └── AGENTS.md           # Remotion 专用指南
│       │   └── utils/
│       │       ├── interpolation.ts
│       │       └── color.ts
│       └── public/                     # 静态资源
│           └── fonts/
│
├── apps/electron/
│   └── src/
│       ├── main/
│       │   └── video/                  # 主进程视频服务
│       │       ├── index.ts
│       │       ├── ipc-handlers.ts     # IPC 处理
│       │       ├── project-manager.ts  # 项目管理
│       │       └── render-worker.ts    # 渲染工作进程
│       │
│       ├── preload/
│       │   └── video-api.ts            # 视频 API 暴露
│       │
│       └── renderer/
│           ├── components/video/       # 视频编辑器 UI
│           │   ├── VideoEditor.tsx     # 主编辑器
│           │   ├── VideoPreview.tsx    # 预览组件
│           │   ├── VideoTimeline.tsx   # 时间轴
│           │   ├── VideoProperties.tsx # 属性面板
│           │   ├── VideoProjectList.tsx# 项目列表
│           │   ├── VideoTemplates.tsx  # 模板选择器
│           │   └── VideoExport.tsx     # 导出对话框
│           │
│           └── pages/
│               └── VideoPage.tsx       # 视频创作页面
```

## 三、核心依赖

### package.json (packages/video)

```json
{
  "name": "@sprouty-ai/video",
  "version": "0.1.0",
  "type": "module",
  "main": "dist/index.js",
  "types": "dist/index.d.ts",
  "scripts": {
    "dev": "remotion studio",
    "build": "tsc",
    "render": "remotion render",
    "upgrade": "remotion upgrade"
  },
  "dependencies": {
    "remotion": "^4.0.0",
    "@remotion/player": "^4.0.0",
    "@remotion/cli": "^4.0.0",
    "@remotion/renderer": "^4.0.0",
    "@remotion/bundler": "^4.0.0",
    "@remotion/zod-types": "^4.0.0",
    "zod": "^3.23.0"
  },
  "peerDependencies": {
    "react": "^18.0.0",
    "react-dom": "^18.0.0"
  }
}
```

## 四、实现细节

### 4.1 视频项目数据模型

```typescript
// packages/video/src/types.ts
import { z } from 'zod';

export const VideoProjectSchema = z.object({
  id: z.string(),
  name: z.string(),
  description: z.string().optional(),
  createdAt: z.string(),
  updatedAt: z.string(),
  
  // 视频配置
  config: z.object({
    width: z.number().default(1920),
    height: z.number().default(1080),
    fps: z.number().default(30),
    durationInFrames: z.number(),
  }),
  
  // 组件代码（AI 生成）
  compositions: z.array(z.object({
    id: z.string(),
    name: z.string(),
    code: z.string(),           // React 组件代码
    props: z.record(z.any()),   // 组件 props
  })),
  
  // 素材资源
  assets: z.array(z.object({
    id: z.string(),
    type: z.enum(['image', 'video', 'audio', 'font']),
    name: z.string(),
    path: z.string(),
  })),
  
  // 渲染历史
  renders: z.array(z.object({
    id: z.string(),
    compositionId: z.string(),
    outputPath: z.string(),
    status: z.enum(['pending', 'rendering', 'completed', 'failed']),
    progress: z.number(),
    createdAt: z.string(),
  })),
});

export type VideoProject = z.infer<typeof VideoProjectSchema>;
```

### 4.2 AI Agent 视频创作技能

```typescript
// packages/video/src/skills/video-creation.ts
import type { Tool } from '@sprouty-ai/shared';

export const videoCreationTools: Tool[] = [
  {
    name: 'video_create_project',
    description: '创建新的视频项目',
    inputSchema: {
      type: 'object',
      properties: {
        name: { 
          type: 'string', 
          description: '项目名称' 
        },
        template: { 
          type: 'string', 
          enum: ['blank', 'title-animation', 'slideshow', 'product-showcase', 'social-media'],
          description: '使用的模板' 
        },
        width: { type: 'number', default: 1920 },
        height: { type: 'number', default: 1080 },
        fps: { type: 'number', default: 30 },
        durationInSeconds: { type: 'number', description: '视频时长（秒）' },
      },
      required: ['name', 'durationInSeconds'],
    },
  },
  
  {
    name: 'video_generate_composition',
    description: '使用 AI 生成 Remotion 视频组件代码',
    inputSchema: {
      type: 'object',
      properties: {
        projectId: { type: 'string' },
        description: { 
          type: 'string', 
          description: '视频内容描述，AI 将据此生成代码' 
        },
        style: { 
          type: 'string', 
          enum: ['modern', 'minimal', 'playful', 'corporate', 'cinematic'],
          description: '视觉风格' 
        },
        colorScheme: {
          type: 'object',
          properties: {
            primary: { type: 'string' },
            secondary: { type: 'string' },
            background: { type: 'string' },
          },
        },
      },
      required: ['projectId', 'description'],
    },
  },
  
  {
    name: 'video_update_composition',
    description: '修改现有视频组件',
    inputSchema: {
      type: 'object',
      properties: {
        projectId: { type: 'string' },
        compositionId: { type: 'string' },
        changes: { 
          type: 'string', 
          description: '需要修改的内容描述' 
        },
      },
      required: ['projectId', 'compositionId', 'changes'],
    },
  },
  
  {
    name: 'video_preview',
    description: '启动视频预览',
    inputSchema: {
      type: 'object',
      properties: {
        projectId: { type: 'string' },
        compositionId: { type: 'string' },
        startFrame: { type: 'number', default: 0 },
      },
      required: ['projectId'],
    },
  },
  
  {
    name: 'video_render',
    description: '渲染视频为 MP4 文件',
    inputSchema: {
      type: 'object',
      properties: {
        projectId: { type: 'string' },
        compositionId: { type: 'string' },
        outputFormat: { 
          type: 'string', 
          enum: ['mp4', 'webm', 'gif'],
          default: 'mp4' 
        },
        quality: {
          type: 'string',
          enum: ['draft', 'standard', 'high'],
          default: 'standard',
        },
      },
      required: ['projectId', 'compositionId'],
    },
  },
  
  {
    name: 'video_add_asset',
    description: '添加素材到视频项目',
    inputSchema: {
      type: 'object',
      properties: {
        projectId: { type: 'string' },
        assetPath: { type: 'string', description: '素材文件路径' },
        assetType: { type: 'string', enum: ['image', 'video', 'audio', 'font'] },
      },
      required: ['projectId', 'assetPath', 'assetType'],
    },
  },
];
```

### 4.3 Remotion Agent Skills 文件

```markdown
<!-- packages/video/src/skills/AGENTS.md -->
# Remotion 视频创作指南

## 概述

本技能用于在 Sprouty AI 中使用 Remotion 创建程序化视频。

## Remotion 核心概念

### 1. Composition（组合）
视频的基本单位，定义了视频的尺寸、帧率和时长。

### 2. useCurrentFrame()
获取当前帧数，用于创建动画。

### 3. useVideoConfig()
获取视频配置（fps、width、height、durationInFrames）。

### 4. interpolate()
创建平滑的数值过渡动画。

## 代码生成规则

### 必须遵循的规则

1. **始终使用 TypeScript**
2. **使用函数组件和 hooks**
3. **避免使用 useEffect 和 useState**（Remotion 有自己的状态管理）
4. **使用 interpolate() 创建动画**
5. **使用 AbsoluteFill 作为根容器**

### 动画最佳实践

```typescript
import { useCurrentFrame, interpolate, Easing } from 'remotion';

// 好的做法：使用 interpolate
const opacity = interpolate(
  frame,
  [0, 30],      // 输入范围
  [0, 1],       // 输出范围
  { extrapolateRight: 'clamp' }
);

// 好的做法：使用 spring 动画
import { spring, useVideoConfig } from 'remotion';

const { fps } = useVideoConfig();
const scale = spring({
  frame,
  fps,
  config: { damping: 200 },
});
```

### 常用组件模式

```typescript
// 标题动画
export const TitleAnimation: React.FC<{ title: string }> = ({ title }) => {
  const frame = useCurrentFrame();
  const opacity = interpolate(frame, [0, 30], [0, 1], { extrapolateRight: 'clamp' });
  const y = interpolate(frame, [0, 30], [50, 0], { extrapolateRight: 'clamp' });
  
  return (
    <AbsoluteFill style={{ justifyContent: 'center', alignItems: 'center' }}>
      <h1 style={{ opacity, transform: `translateY(${y}px)` }}>
        {title}
      </h1>
    </AbsoluteFill>
  );
};

// 序列动画
import { Sequence } from 'remotion';

export const VideoWithSequences: React.FC = () => {
  return (
    <AbsoluteFill>
      <Sequence from={0} durationInFrames={60}>
        <Intro />
      </Sequence>
      <Sequence from={60} durationInFrames={120}>
        <MainContent />
      </Sequence>
      <Sequence from={180} durationInFrames={60}>
        <Outro />
      </Sequence>
    </AbsoluteFill>
  );
};
```

## 模板使用

### 可用模板

1. **title-animation** - 标题动画
2. **slideshow** - 图片轮播
3. **product-showcase** - 产品展示
4. **social-media** - 社交媒体短视频
5. **data-visualization** - 数据可视化

### 模板参数

每个模板都接受标准化的 props：

```typescript
interface TemplateProps {
  title?: string;
  subtitle?: string;
  items?: Array<{ title: string; description?: string; image?: string }>;
  colors?: {
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;
  cta?: { text: string; url?: string };
}
```

## 渲染配置

### 质量预设

- **draft**: 快速预览，720p，15fps
- **standard**: 标准质量，1080p，30fps
- **high**: 高质量，1080p/4K，60fps

### 输出格式

- **mp4**: H.264 编码，最通用
- **webm**: VP9 编码，适合 Web
- **gif**: 动图，适合短循环动画
```

### 4.4 视频编辑器 UI 组件

```typescript
// apps/electron/src/renderer/components/video/VideoEditor.tsx
import React, { useState, useCallback } from 'react';
import { Player, PlayerRef } from '@remotion/player';
import { VideoPreview } from './VideoPreview';
import { VideoTimeline } from './VideoTimeline';
import { VideoProperties } from './VideoProperties';
import { VideoProjectList } from './VideoProjectList';
import type { VideoProject } from '@sprouty-ai/video';

interface VideoEditorProps {
  workspaceId: string;
}

export const VideoEditor: React.FC<VideoEditorProps> = ({ workspaceId }) => {
  const [currentProject, setCurrentProject] = useState<VideoProject | null>(null);
  const [currentFrame, setCurrentFrame] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const playerRef = React.useRef<PlayerRef>(null);

  const handleFrameChange = useCallback((frame: number) => {
    setCurrentFrame(frame);
    playerRef.current?.seekTo(frame);
  }, []);

  const handleRender = useCallback(async () => {
    if (!currentProject) return;
    
    // 调用主进程渲染
    await window.electronAPI.video.render({
      projectId: currentProject.id,
      compositionId: currentProject.compositions[0]?.id,
      quality: 'standard',
    });
  }, [currentProject]);

  return (
    <div className="flex h-full">
      {/* 左侧：项目列表 */}
      <div className="w-64 border-r border-border">
        <VideoProjectList
          workspaceId={workspaceId}
          onSelect={setCurrentProject}
          selected={currentProject?.id}
        />
      </div>

      {/* 中间：预览区域 */}
      <div className="flex-1 flex flex-col">
        <div className="flex-1 p-4">
          <VideoPreview
            project={currentProject}
            playerRef={playerRef}
            onFrameChange={setCurrentFrame}
          />
        </div>
        
        {/* 时间轴 */}
        <div className="h-32 border-t border-border">
          <VideoTimeline
            project={currentProject}
            currentFrame={currentFrame}
            onFrameChange={handleFrameChange}
            isPlaying={isPlaying}
            onPlayPause={() => setIsPlaying(!isPlaying)}
          />
        </div>
      </div>

      {/* 右侧：属性面板 */}
      <div className="w-80 border-l border-border">
        <VideoProperties
          project={currentProject}
          onUpdate={(updates) => {
            // 更新项目属性
          }}
          onRender={handleRender}
        />
      </div>
    </div>
  );
};
```

### 4.5 主进程渲染服务

```typescript
// apps/electron/src/main/video/render-worker.ts
import { bundle } from '@remotion/bundler';
import { renderMedia, selectComposition, getCompositions } from '@remotion/renderer';
import { BrowserWindow } from 'electron';
import path from 'path';

interface RenderOptions {
  projectPath: string;
  compositionId: string;
  outputPath: string;
  quality: 'draft' | 'standard' | 'high';
  onProgress?: (progress: number) => void;
}

const QUALITY_PRESETS = {
  draft: { crf: 28, scale: 0.5 },
  standard: { crf: 18, scale: 1 },
  high: { crf: 12, scale: 1 },
};

export class VideoRenderWorker {
  private mainWindow: BrowserWindow;

  constructor(mainWindow: BrowserWindow) {
    this.mainWindow = mainWindow;
  }

  async render(options: RenderOptions): Promise<string> {
    const { projectPath, compositionId, outputPath, quality, onProgress } = options;
    const preset = QUALITY_PRESETS[quality];

    try {
      // 1. Bundle 项目
      this.sendProgress('bundling', 0);
      const bundled = await bundle({
        entryPoint: path.join(projectPath, 'src/index.ts'),
        webpackOverride: (config) => config,
      });

      // 2. 获取 composition
      this.sendProgress('preparing', 10);
      const compositions = await getCompositions(bundled);
      const composition = compositions.find(c => c.id === compositionId);
      
      if (!composition) {
        throw new Error(`Composition "${compositionId}" not found`);
      }

      // 3. 渲染视频
      this.sendProgress('rendering', 20);
      await renderMedia({
        composition,
        serveUrl: bundled,
        codec: 'h264',
        outputLocation: outputPath,
        crf: preset.crf,
        scale: preset.scale,
        onProgress: ({ progress }) => {
          const totalProgress = 20 + progress * 80;
          this.sendProgress('rendering', totalProgress);
          onProgress?.(totalProgress);
        },
      });

      this.sendProgress('completed', 100);
      return outputPath;
    } catch (error) {
      this.sendProgress('failed', 0, error.message);
      throw error;
    }
  }

  private sendProgress(status: string, progress: number, error?: string) {
    this.mainWindow.webContents.send('video:render-progress', {
      status,
      progress,
      error,
    });
  }
}
```

### 4.6 IPC 通道定义

```typescript
// apps/electron/src/main/video/ipc-handlers.ts
import { ipcMain, dialog } from 'electron';
import { VideoRenderWorker } from './render-worker';
import { VideoProjectManager } from './project-manager';

export function registerVideoIpcHandlers(
  renderWorker: VideoRenderWorker,
  projectManager: VideoProjectManager
) {
  // 创建项目
  ipcMain.handle('video:create-project', async (_, options) => {
    return projectManager.createProject(options);
  });

  // 获取项目列表
  ipcMain.handle('video:list-projects', async (_, workspaceId) => {
    return projectManager.listProjects(workspaceId);
  });

  // 获取项目详情
  ipcMain.handle('video:get-project', async (_, projectId) => {
    return projectManager.getProject(projectId);
  });

  // 更新项目
  ipcMain.handle('video:update-project', async (_, projectId, updates) => {
    return projectManager.updateProject(projectId, updates);
  });

  // 删除项目
  ipcMain.handle('video:delete-project', async (_, projectId) => {
    return projectManager.deleteProject(projectId);
  });

  // 添加素材
  ipcMain.handle('video:add-asset', async (_, projectId, assetPath, assetType) => {
    return projectManager.addAsset(projectId, assetPath, assetType);
  });

  // 渲染视频
  ipcMain.handle('video:render', async (_, options) => {
    const { projectId, compositionId, quality } = options;
    const project = await projectManager.getProject(projectId);
    
    // 选择输出路径
    const { filePath } = await dialog.showSaveDialog({
      title: '保存视频',
      defaultPath: `${project.name}.mp4`,
      filters: [
        { name: 'MP4 Video', extensions: ['mp4'] },
        { name: 'WebM Video', extensions: ['webm'] },
      ],
    });

    if (!filePath) return null;

    return renderWorker.render({
      projectPath: project.path,
      compositionId,
      outputPath: filePath,
      quality,
    });
  });

  // 启动预览服务器
  ipcMain.handle('video:start-preview', async (_, projectId) => {
    return projectManager.startPreview(projectId);
  });

  // 停止预览服务器
  ipcMain.handle('video:stop-preview', async (_, projectId) => {
    return projectManager.stopPreview(projectId);
  });
}
```

## 五、实施计划

### Phase 1: 基础架构（2 周）

- [ ] 创建 `@sprouty-ai/video` 包
- [ ] 安装 Remotion 依赖
- [ ] 实现基础 Composition 组件
- [ ] 配置 Remotion Studio 开发环境

### Phase 2: 渲染引擎（2 周）

- [ ] 实现本地渲染器
- [ ] 实现预览服务器
- [ ] 添加 IPC 通道
- [ ] 实现渲染进度追踪

### Phase 3: UI 集成（3 周）

- [ ] 实现 VideoEditor 组件
- [ ] 实现 VideoPreview 组件
- [ ] 实现 VideoTimeline 组件
- [ ] 实现 VideoProperties 组件
- [ ] 添加视频创作页面路由

### Phase 4: AI 集成（3 周）

- [ ] 实现视频创作 Agent Skills
- [ ] 添加 AGENTS.md 指南
- [ ] 实现代码生成逻辑
- [ ] 实现模板系统
- [ ] 测试 AI 生成视频流程

### Phase 5: 优化与完善（2 周）

- [ ] 性能优化
- [ ] 错误处理
- [ ] i18n 支持
- [ ] 文档完善
- [ ] 用户测试

## 六、许可证说明

### Remotion 许可证

根据 [Remotion 许可证](https://www.remotion.dev/docs/license)：

- ✅ **免费使用**：个人用户、3人以下公司、非营利组织
- 💰 **需要付费**：超过 3 人的营利性公司

### 建议

1. 在 Sprouty AI 设置中添加许可证配置选项
2. 首次使用时提示用户确认许可证状态
3. 提供购买许可证的链接

## 七、风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| 渲染性能慢 | 用户体验差 | 后台渲染 + 进度显示 + 低质量预览 |
| AI 生成代码错误 | 无法渲染 | 代码验证 + 沙箱测试 + 错误提示 |
| 素材管理复杂 | 使用门槛高 | 拖放上传 + 素材库 + 模板预设 |
| Electron 打包体积增大 | 安装包变大 | 按需加载 + 延迟初始化 |

## 八、参考资源

- [Remotion 官方文档](https://www.remotion.dev/docs)
- [Remotion AI 集成](https://www.remotion.dev/docs/ai)
- [Remotion Agent Skills](https://www.remotion.dev/docs/ai/skills)
- [Remotion GitHub](https://github.com/remotion-dev/remotion)
- [Remotion 示例](https://www.remotion.dev/showcase)


## 九、快速启动指南

### 9.1 安装依赖

```bash
cd Sprouty AI

# 安装 Remotion 核心依赖
bun add remotion @remotion/cli @remotion/player @remotion/renderer @remotion/bundler

# 安装可选依赖
bun add @remotion/zod-types @remotion/media-utils
```

### 9.2 创建视频包

```bash
# 创建包目录
mkdir -p packages/video/src/{compositions,templates,components,hooks,renderer,skills,utils}

# 初始化 package.json
cd packages/video
bun init
```

### 9.3 基础配置文件

```typescript
// packages/video/remotion.config.ts
import { Config } from '@remotion/cli/config';

Config.setVideoImageFormat('jpeg');
Config.setOverwriteOutput(true);
```

```typescript
// packages/video/src/index.ts
export * from './compositions';
export * from './templates';
export * from './renderer';
export * from './skills';
export type * from './types';
```

```typescript
// packages/video/src/Root.tsx
import React from 'react';
import { Composition } from 'remotion';
import { HelloWorld } from './compositions/HelloWorld';
import { TitleAnimation } from './compositions/TitleAnimation';

export const RemotionRoot: React.FC = () => {
  return (
    <>
      <Composition
        id="HelloWorld"
        component={HelloWorld}
        durationInFrames={150}
        fps={30}
        width={1920}
        height={1080}
      />
      <Composition
        id="TitleAnimation"
        component={TitleAnimation}
        durationInFrames={90}
        fps={30}
        width={1920}
        height={1080}
        defaultProps={{
          title: 'Sprouty AI',
          subtitle: 'AI 创作平台',
        }}
      />
    </>
  );
};
```

### 9.4 示例 Composition

```typescript
// packages/video/src/compositions/HelloWorld.tsx
import React from 'react';
import { AbsoluteFill, useCurrentFrame, useVideoConfig, interpolate } from 'remotion';

export const HelloWorld: React.FC = () => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();

  const opacity = interpolate(frame, [0, fps], [0, 1], {
    extrapolateRight: 'clamp',
  });

  const scale = interpolate(frame, [0, fps], [0.8, 1], {
    extrapolateRight: 'clamp',
  });

  return (
    <AbsoluteFill
      style={{
        backgroundColor: '#0f0f0f',
        justifyContent: 'center',
        alignItems: 'center',
      }}
    >
      <h1
        style={{
          color: '#ffffff',
          fontSize: 120,
          fontWeight: 'bold',
          opacity,
          transform: `scale(${scale})`,
        }}
      >
        Hello Sprouty AI
      </h1>
    </AbsoluteFill>
  );
};
```

```typescript
// packages/video/src/compositions/TitleAnimation.tsx
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Sequence,
} from 'remotion';

interface TitleAnimationProps {
  title: string;
  subtitle?: string;
  backgroundColor?: string;
  textColor?: string;
}

export const TitleAnimation: React.FC<TitleAnimationProps> = ({
  title,
  subtitle,
  backgroundColor = '#1a1a2e',
  textColor = '#ffffff',
}) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();

  // 标题动画
  const titleOpacity = interpolate(frame, [0, 20], [0, 1], {
    extrapolateRight: 'clamp',
  });

  const titleY = spring({
    frame,
    fps,
    config: { damping: 200, stiffness: 100 },
  });

  // 副标题动画（延迟出现）
  const subtitleOpacity = interpolate(frame, [20, 40], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  return (
    <AbsoluteFill
      style={{
        backgroundColor,
        justifyContent: 'center',
        alignItems: 'center',
        flexDirection: 'column',
        gap: 20,
      }}
    >
      <h1
        style={{
          color: textColor,
          fontSize: 100,
          fontWeight: 'bold',
          opacity: titleOpacity,
          transform: `translateY(${(1 - titleY) * 50}px)`,
          margin: 0,
        }}
      >
        {title}
      </h1>

      {subtitle && (
        <p
          style={{
            color: textColor,
            fontSize: 40,
            opacity: subtitleOpacity,
            margin: 0,
          }}
        >
          {subtitle}
        </p>
      )}
    </AbsoluteFill>
  );
};
```

### 9.5 验证安装

```bash
# 在 packages/video 目录下
cd packages/video

# 启动 Remotion Studio 预览
bun run dev

# 渲染测试视频
bun run render HelloWorld out/test.mp4
```

## 十、与 AI Agent 集成示例

### 10.1 对话示例

```
用户: 帮我创建一个 10 秒的产品介绍视频，标题是"智能创作助手"，
     包含 3 个特性展示：AI 写作、视频生成、多平台发布

AI Agent:
1. 创建视频项目
   → video_create_project({ name: "智能创作助手介绍", durationInSeconds: 10 })

2. 生成视频组件
   → video_generate_composition({
       projectId: "proj_xxx",
       description: "产品介绍视频，标题'智能创作助手'，3个特性：AI写作、视频生成、多平台发布",
       style: "modern"
     })

3. 启动预览
   → video_preview({ projectId: "proj_xxx" })

4. 用户确认后渲染
   → video_render({ projectId: "proj_xxx", quality: "high" })
```

### 10.2 生成的代码示例

```typescript
// AI 自动生成的 Composition
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Sequence,
} from 'remotion';

const features = [
  { icon: '✍️', title: 'AI 写作', desc: '智能生成高质量内容' },
  { icon: '🎬', title: '视频生成', desc: '一键创建专业视频' },
  { icon: '📱', title: '多平台发布', desc: '同步发布到各大平台' },
];

export const ProductIntro: React.FC = () => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  return (
    <AbsoluteFill style={{ backgroundColor: '#0f172a' }}>
      {/* 标题序列 */}
      <Sequence from={0} durationInFrames={fps * 3}>
        <TitleSection />
      </Sequence>

      {/* 特性展示序列 */}
      {features.map((feature, index) => (
        <Sequence
          key={index}
          from={fps * 3 + index * fps * 2}
          durationInFrames={fps * 2}
        >
          <FeatureCard {...feature} />
        </Sequence>
      ))}

      {/* 结尾 CTA */}
      <Sequence from={fps * 9} durationInFrames={fps}>
        <CTASection />
      </Sequence>
    </AbsoluteFill>
  );
};

const TitleSection: React.FC = () => {
  const frame = useCurrentFrame();
  const opacity = interpolate(frame, [0, 30], [0, 1]);
  const scale = interpolate(frame, [0, 30], [0.9, 1]);

  return (
    <AbsoluteFill style={{ justifyContent: 'center', alignItems: 'center' }}>
      <h1
        style={{
          color: '#fff',
          fontSize: 80,
          fontWeight: 'bold',
          opacity,
          transform: `scale(${scale})`,
        }}
      >
        智能创作助手
      </h1>
    </AbsoluteFill>
  );
};

const FeatureCard: React.FC<{ icon: string; title: string; desc: string }> = ({
  icon,
  title,
  desc,
}) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();

  const slideIn = spring({ frame, fps, config: { damping: 200 } });

  return (
    <AbsoluteFill style={{ justifyContent: 'center', alignItems: 'center' }}>
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          gap: 20,
          transform: `translateX(${(1 - slideIn) * 100}px)`,
          opacity: slideIn,
        }}
      >
        <span style={{ fontSize: 80 }}>{icon}</span>
        <h2 style={{ color: '#fff', fontSize: 48, margin: 0 }}>{title}</h2>
        <p style={{ color: '#94a3b8', fontSize: 24, margin: 0 }}>{desc}</p>
      </div>
    </AbsoluteFill>
  );
};

const CTASection: React.FC = () => {
  const frame = useCurrentFrame();
  const opacity = interpolate(frame, [0, 15], [0, 1]);

  return (
    <AbsoluteFill style={{ justifyContent: 'center', alignItems: 'center' }}>
      <div
        style={{
          padding: '20px 40px',
          backgroundColor: '#3b82f6',
          borderRadius: 12,
          opacity,
        }}
      >
        <span style={{ color: '#fff', fontSize: 32, fontWeight: 'bold' }}>
          立即体验 →
        </span>
      </div>
    </AbsoluteFill>
  );
};
```

## 十一、总结

本方案提供了一个完整的 Remotion 集成路径，核心优势：

1. **深度集成** - 作为 Sprouty AI 的原生功能，而非外部工具
2. **AI 驱动** - 利用 Claude Agent SDK 实现自然语言视频创作
3. **模板系统** - 预设模板 + 自定义模板，降低使用门槛
4. **渐进式实施** - 分阶段实施，快速验证价值

预计 MVP 开发周期：**8-10 周**，可实现基础的 AI 视频创作能力。
