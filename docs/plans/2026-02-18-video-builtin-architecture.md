# 视频创作功能重构设计方案 — 完全内置架构

**状态：** 设计评审中
**日期：** 2026-02-18
**作者：** Claude (Brainstorming Session)
**前置文档：** scene-composer-design-v3.rst

---

## 一、问题分析

### 1.1 现有架构的根本问题

当前视频功能被设计成一个**独立的 MCP Server 子进程**，导致一系列问题：

1. **MCP Server 需要 Bun 运行时** → 用户机器上不一定有，打包分发复杂
2. **stdio 通信** → 启动慢、容易崩溃、调试困难
3. **项目数据存在 JSON 文件** → 和主应用的 SQLite (creator.db) 割裂
4. **前端需要先"启动视频服务"** → 多了一个不必要的步骤
5. **依赖安装困难** → Remotion 环境需要在用户项目目录安装，Windows 兼容性差

### 1.2 用户体验断裂

- 用户需要理解"启动视频服务"的概念
- 服务启动失败时没有清晰的错误提示
- 项目管理、素材管理、编辑、导出功能在 UI 上不完整
- Agent 对话和 UI 操作之间数据不同步

### 1.3 设计目标

**核心原则：** 视频功能应该是 Electron 应用的内置能力，而不是外部服务。

**用户体验目标：**
- 启动应用即可使用，无需额外配置
- Agent 驱动创作，UI 提供预览、微调和导出
- 数据统一存储在 creator.db
- 跨平台兼容（Windows/Mac/Linux）

---

## 二、新架构设计

### 2.1 架构总览

```
旧架构:
  Electron App → (spawn) → Bun MCP Video Server → (stdio) → 操作 project.json
  前端 → IPC → main → MCP Server → 返回结果

新架构:
  Electron App (内置一切)
  ├── Agent → session-scoped tools (video_*) → IPC → main 进程 → 操作 SQLite
  ├── 前端 UI → IPC → main 进程 → 操作 SQLite
  └── 渲染: 预览用 @remotion/player (renderer 进程)
            导出用 @remotion/renderer (main 进程, 内嵌 ffmpeg, 复用 Electron Chromium)
```

### 2.2 核心变化

| 模块 | 旧方案 | 新方案 |
|------|--------|--------|
| **数据存储** | project.json 文件 | SQLite (creator.db) |
| **Agent 工具** | MCP Server (独立进程) | Session-scoped tools (进程内) |
| **预览** | Preview Server (子进程) | @remotion/player (renderer 进程) |
| **渲染** | MCP Server → RenderEngine | RenderWorker (main 进程) |
| **依赖** | 用户项目目录需要 Remotion 环境 | 内嵌 ffmpeg，复用 Electron Chromium |

### 2.3 删除的组件

- `packages/video/src/mcp-server/` — 整个 MCP server 目录
- `apps/electron/src/main/video/service-manager.ts` — MCP 进程管理
- `apps/electron/src/main/video/preview-server.ts` — 独立预览服务器
- `apps/electron/src/main/video/project-manager.ts` — 项目管理（移到 Repository）
- 所有 `video:start-service` / `video:stop-service` 相关逻辑

### 2.4 保留的组件

- `packages/video/src/compositions/` — 所有 Remotion 组合组件（SceneComposer 等）
- `packages/video/src/components/` — 可复用的动画组件
- `packages/video/src/templates/` — 模板定义（需修改格式）
- `packages/video/src/Root.tsx` — Remotion 入口（导出渲染用）
- `apps/electron/src/main/video/render-worker.ts` — 渲染 worker（需重构）
- `apps/electron/src/renderer/components/video/` — 前端组件（需完善）

---

## 三、数据层设计

### 3.1 SQLite Schema

新增 3 张表到 creator.db：

```sql
-- 视频项目（和 contents 表关联）
CREATE TABLE video_projects (
  id TEXT PRIMARY KEY,
  content_id TEXT NOT NULL REFERENCES contents(id) ON DELETE CASCADE,
  name TEXT NOT NULL,
  description TEXT,
  width INTEGER NOT NULL DEFAULT 1080,
  height INTEGER NOT NULL DEFAULT 1920,
  fps INTEGER NOT NULL DEFAULT 30,
  metadata TEXT,          -- JSON: 扩展字段
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_video_projects_content_id ON video_projects(content_id);

-- 视频场景（有序列表）
CREATE TABLE video_scenes (
  id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL REFERENCES video_projects(id) ON DELETE CASCADE,
  composition_id TEXT NOT NULL,   -- 引用内置组合 ID，如 "TitleAnimation"
  name TEXT,
  sort_order INTEGER NOT NULL DEFAULT 0,
  duration_in_frames INTEGER NOT NULL DEFAULT 90,
  props TEXT NOT NULL DEFAULT '{}',  -- JSON: 传给组合组件的参数
  transition_type TEXT DEFAULT 'none' CHECK(transition_type IN ('none','fade','slide','wipe','flip','clock-wipe')),
  transition_duration INTEGER DEFAULT 0,
  transition_direction TEXT CHECK(transition_direction IN ('from-left','from-right','from-top','from-bottom') OR transition_direction IS NULL),
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_video_scenes_project_id ON video_scenes(project_id);
CREATE INDEX idx_video_scenes_sort_order ON video_scenes(project_id, sort_order);

-- 视频素材
CREATE TABLE video_assets (
  id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL REFERENCES video_projects(id) ON DELETE CASCADE,
  type TEXT NOT NULL CHECK(type IN ('image','video','audio','font')),
  name TEXT NOT NULL,
  file_path TEXT NOT NULL,   -- 相对于项目目录的路径
  file_size INTEGER,
  metadata TEXT,             -- JSON: 宽高、时长等
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_video_assets_project_id ON video_assets(project_id);
CREATE INDEX idx_video_assets_type ON video_assets(project_id, type);
```

### 3.2 设计要点

1. **video_projects.content_id** 关联到 contents 表 — 一个内容条目对应一个视频项目
2. **过渡效果存在 video_scenes 表** — 每个场景存它"前面"的过渡，而不是单独的 transitions 数组
3. **sort_order** 控制场景顺序，方便拖拽排序
4. **props 是 JSON 字段** — 存储每个组合的参数（文字、颜色、图片路径等）
5. **级联删除** — 删除 content 时自动删除关联的 video_project、scenes、assets

### 3.3 TypeScript 类型

```typescript
export interface VideoProject {
  id: string;
  content_id: string;
  name: string;
  description?: string;
  width: number;
  height: number;
  fps: number;
  metadata?: string;
  created_at: string;
  updated_at: string;

  // 关联数据（查询时填充）
  scenes?: VideoScene[];
  assets?: VideoAsset[];
}

export interface VideoScene {
  id: string;
  project_id: string;
  composition_id: string;
  name?: string;
  sort_order: number;
  duration_in_frames: number;
  props: Record<string, any>; // JSON parsed
  transition_type: 'none' | 'fade' | 'slide' | 'wipe' | 'flip' | 'clock-wipe';
  transition_duration: number;
  transition_direction?: 'from-left' | 'from-right' | 'from-top' | 'from-bottom';
  created_at: string;
  updated_at: string;
}

export interface VideoAsset {
  id: string;
  project_id: string;
  type: 'image' | 'video' | 'audio' | 'font';
  name: string;
  file_path: string;
  file_size?: number;
  metadata?: string;
  created_at: string;
}
```

---

## 四、Agent 工具层设计

### 4.1 Session-Scoped Tools

按照现有的 session-scoped tools 模式，在 `packages/shared/src/agent/session-scoped-tools.ts` 中新增 10 个视频工具：

| 工具名 | 功能 | 参数 |
|--------|------|------|
| `video_list_compositions` | 列出所有可用的内置组合 | 无 |
| `video_create_project` | 创建视频项目 | contentId, name, width, height, fps, description |
| `video_add_scene` | 添加场景 | projectId, compositionId, name, durationInFrames, props, insertAt, transition* |
| `video_update_scene` | 更新场景 | projectId, sceneId, props, durationInFrames, compositionId, transition* |
| `video_remove_scene` | 删除场景 | projectId, sceneId |
| `video_reorder_scenes` | 重新排序场景 | projectId, sceneIds[] |
| `video_get_project` | 获取项目详情 | projectId |
| `video_add_asset` | 添加素材 | projectId, assetPath, assetType |
| `video_list_templates` | 列出模板 | category? |
| `video_create_from_template` | 从模板创建 | contentId, templateId, name |

### 4.2 工具实现模式

```typescript
export function createVideoAddSceneTool(sessionId: string, workspaceRootPath: string) {
  return tool(
    'video_add_scene',
    'Add a scene to a video project',
    {
      projectId: z.string(),
      compositionId: z.string().describe('Built-in composition ID (e.g., "TitleAnimation")'),
      name: z.string().optional(),
      durationInFrames: z.number().describe('Scene duration in frames'),
      props: z.record(z.any()).describe('Props to pass to the composition component'),
      insertAt: z.number().optional().describe('Insert position (default: append to end)'),
      transitionType: z.enum(['none','fade','slide','wipe','flip','clock-wipe']).optional(),
      transitionDuration: z.number().optional(),
      transitionDirection: z.enum(['from-left','from-right','from-top','from-bottom']).optional(),
    },
    async (args) => {
      const result = await window.electronAPI.video.addScene({
        projectId: args.projectId,
        compositionId: args.compositionId,
        name: args.name,
        durationInFrames: args.durationInFrames,
        props: args.props,
        insertAt: args.insertAt,
        transitionType: args.transitionType || 'none',
        transitionDuration: args.transitionDuration || 0,
        transitionDirection: args.transitionDirection,
      });
      return { content: [{ type: 'text', text: `Scene added: ${result.sceneId}` }] };
    }
  );
}
```

### 4.3 注册到 Agent

```typescript
// 在 getSessionScopedTools() 函数中添加
export function getSessionScopedTools(sessionId: string, workspaceRootPath: string) {
  const tools = [
    createSubmitPlanTool(sessionId),
    // ... 其他工具

    // 视频工具
    createVideoListCompositionsTool(sessionId, workspaceRootPath),
    createVideoCreateProjectTool(sessionId, workspaceRootPath),
    createVideoAddSceneTool(sessionId, workspaceRootPath),
    createVideoUpdateSceneTool(sessionId, workspaceRootPath),
    createVideoRemoveSceneTool(sessionId, workspaceRootPath),
    createVideoReorderScenesTool(sessionId, workspaceRootPath),
    createVideoGetProjectTool(sessionId, workspaceRootPath),
    createVideoAddAssetTool(sessionId, workspaceRootPath),
    createVideoListTemplatesTool(sessionId, workspaceRootPath),
    createVideoCreateFromTemplateTool(sessionId, workspaceRootPath),
  ];

  return createSdkMcpServer({ name: 'session', version: '1.0.0', tools });
}
```

---

## 五、IPC 层设计

### 5.1 IPC 通道定义

```typescript
export const VIDEO_IPC_CHANNELS = {
  // 项目管理
  CREATE_PROJECT: 'video:create-project',
  GET_PROJECT: 'video:get-project',
  LIST_PROJECTS: 'video:list-projects',
  UPDATE_PROJECT: 'video:update-project',
  DELETE_PROJECT: 'video:delete-project',

  // 场景管理
  ADD_SCENE: 'video:add-scene',
  UPDATE_SCENE: 'video:update-scene',
  REMOVE_SCENE: 'video:remove-scene',
  REORDER_SCENES: 'video:reorder-scenes',
  GET_SCENES: 'video:get-scenes',

  // 素材管理
  SELECT_ASSET_FILES: 'video:select-asset-files',
  ADD_ASSET: 'video:add-asset',
  REMOVE_ASSET: 'video:remove-asset',
  LIST_ASSETS: 'video:list-assets',

  // 组合和模板
  LIST_COMPOSITIONS: 'video:list-compositions',
  LIST_TEMPLATES: 'video:list-templates',
  GET_TEMPLATE: 'video:get-template',
  CREATE_FROM_TEMPLATE: 'video:create-from-template',

  // 渲染
  RENDER: 'video:render',
  CANCEL_RENDER: 'video:cancel-render',
  RENDER_PROGRESS: 'video:render-progress',
  RENDER_COMPLETED: 'video:render-completed',
  RENDER_FAILED: 'video:render-failed',
} as const;
```

### 5.2 核心 Handler 实现

所有 handler 直接调用 `VideoRepository`，不再通过 MCP Server：

```typescript
// 创建项目
ipcMain.handle(VIDEO_IPC_CHANNELS.CREATE_PROJECT, async (_event, request: CreateProjectRequest) => {
  const videoRepo = new VideoRepository(getCreatorDb());
  return videoRepo.createProject({
    contentId: request.contentId,
    name: request.name,
    description: request.description,
    width: request.config?.width || 1080,
    height: request.config?.height || 1920,
    fps: request.config?.fps || 30,
  });
});

// 添加场景
ipcMain.handle(VIDEO_IPC_CHANNELS.ADD_SCENE, async (_event, request: AddSceneRequest) => {
  const videoRepo = new VideoRepository(getCreatorDb());
  const sceneId = videoRepo.addScene(request);
  return { sceneId };
});

// 获取项目（包含场景列表）
ipcMain.handle(VIDEO_IPC_CHANNELS.GET_PROJECT, async (_event, projectId: string) => {
  const videoRepo = new VideoRepository(getCreatorDb());
  return videoRepo.getProject(projectId);
});
```

---

## 六、渲染引擎设计

### 6.1 RenderWorker 重构

```typescript
export class RenderWorker {
  async render(options: RenderOptions): Promise<string> {
    const { projectId, outputPath, outputFormat, quality, onProgress } = options;

    // 1. 从数据库加载项目数据
    const videoRepo = new VideoRepository(getCreatorDb());
    const project = videoRepo.getProject(projectId);
    if (!project || !project.scenes || project.scenes.length === 0) {
      throw new Error('Project has no scenes');
    }

    // 2. 构造 SceneComposer 的 inputProps
    const inputProps = {
      scenes: project.scenes.map(s => ({
        id: s.id,
        name: s.name || s.composition_id,
        compositionId: s.composition_id,
        durationInFrames: s.duration_in_frames,
        props: s.props,
      })),
      transitions: project.scenes.map(s => ({
        type: s.transition_type || 'none',
        durationInFrames: s.transition_duration || 0,
        direction: s.transition_direction,
      })),
    };

    // 3. 获取 Remotion 入口文件路径
    const entryPoint = this.getEntryPoint();

    // 4. Bundle Remotion 项目
    const bundled = await bundle({ entryPoint });

    // 5. 选择 SceneComposer 组合
    const composition = await selectComposition({
      serveUrl: bundled,
      id: 'SceneComposer',
      inputProps,
    });

    // 6. 渲染视频
    await renderMedia({
      composition,
      serveUrl: bundled,
      codec: outputFormat === 'mp4' ? 'h264' : outputFormat,
      outputLocation: outputPath,
      browserExecutable: this.getElectronChromiumPath(),
      ffmpegExecutable: this.getFfmpegPath(),
      onProgress: ({ progress }) => {
        onProgress?.({ progress: Math.round(progress * 100) });
      },
      cancelSignal: this.abortController.signal,
    });

    return outputPath;
  }

  private getElectronChromiumPath(): string | undefined {
    return process.execPath; // 尝试复用 Electron Chromium
  }

  private getFfmpegPath(): string | undefined {
    if (app.isPackaged) {
      return path.join(process.resourcesPath, 'bin', 'ffmpeg');
    } else {
      return require('ffmpeg-static');
    }
  }
}
```

### 6.2 依赖处理

- **ffmpeg**: 开发环境用 `ffmpeg-static`，打包时捆绑到 `resources/bin/`
- **Chromium**: 优先复用 Electron 内置 Chromium，失败时 Remotion 自动下载

---

## 七、前端 UI 设计

### 7.1 VideoStudio 页面结构

三列布局：

```
┌─────────────────────────────────────────────────────────────┐
│ 左侧面板 (300px)    │ 中心面板 (flex-1)  │ 右侧面板 (350px) │
├─────────────────────┼────────────────────┼──────────────────┤
│ [Tab: 内容 | 模板]  │ 视频预览           │ 属性编辑         │
│                     │ (@remotion/player) │ - 场景参数       │
│ 视频内容列表        │                    │ - 过渡效果       │
│ - 项目1 (creating)  │                    │ - 时长           │
│ - 项目2 (completed) │                    │                  │
│                     ├────────────────────┤ 素材管理         │
│ 或                  │ 场景列表 (拖拽)    │ - 图片           │
│                     │ ┌────────────────┐ │ - 视频           │
│ 模板库              │ │ 场景1: 标题    │ │ - 音频           │
│ - 社交媒体          │ │ 90帧 | fade    │ │ - 字体           │
│ - 营销推广          │ └────────────────┘ │                  │
│ - 教程讲解          │ ┌────────────────┐ │ 导出             │
│                     │ │ 场景2: 内容    │ │ [导出视频]       │
│                     │ │ 180帧 | slide  │ │ 格式: MP4        │
│                     │ └────────────────┘ │ 质量: 标准       │
└─────────────────────┴────────────────────┴──────────────────┘
```

### 7.2 核心功能

1. **项目管理**
   - 创建项目（关联到 content）
   - 切换项目
   - 删除项目

2. **场景编辑**
   - 添加场景（选择 composition）
   - 拖拽排序（@dnd-kit）
   - 编辑参数（动态表单）
   - 删除场景

3. **实时预览**
   - @remotion/player 渲染 SceneComposer
   - 参数变化实时更新
   - 播放控制（播放/暂停/跳转）

4. **素材管理**
   - 添加素材（文件选择对话框）
   - 素材列表（缩略图）
   - 删除素材

5. **导出**
   - 选择格式（mp4/webm/gif）
   - 选择质量（draft/standard/high）
   - 进度显示
   - 取消渲染

---

## 八、实施计划

### Phase 1: 数据层基础（2-3 天）

**任务：**
- 添加 3 张表到 schema.ts
- 创建迁移脚本 011_to_012
- 创建 VideoRepository
- 编写单元测试

**验证：**
- creator.db 自动创建新表
- Repository CRUD 操作正常
- 外键约束生效

### Phase 2: IPC 层重构（2-3 天）

**任务：**
- 重写 ipc-handlers.ts
- 删除 service-manager.ts, preview-server.ts
- 更新 preload/video-api.ts

**验证：**
- IPC 调用能正确读写 SQLite
- 前端调用不报错

### Phase 3: Agent 工具层（1-2 天）

**任务：**
- 添加 10 个 session-scoped video tools
- 注册到 getSessionScopedTools()

**验证：**
- Agent 工具列表包含 video_* 工具
- 对话测试能成功调用

### Phase 4: 渲染引擎重构（2-3 天）

**任务：**
- 重构 RenderWorker
- 内嵌 ffmpeg
- 配置打包

**验证：**
- 能成功导出 MP4
- 进度显示正常
- 取消功能正常

### Phase 5: 前端 UI 完善（3-4 天）

**任务：**
- 完善 VideoStudio.tsx
- 实现场景拖拽排序
- 集成 @remotion/player
- 完善所有子组件

**验证：**
- 能创建项目、添加场景
- 预览实时更新
- 能导出视频

### Phase 6: 模板系统集成（1-2 天）

**任务：**
- 修改模板定义格式
- 实现 CREATE_FROM_TEMPLATE

**验证：**
- 从模板创建项目正常

### Phase 7: 清理与优化（1 天）

**任务：**
- 删除废弃代码
- 性能优化
- 更新文档

### Phase 8: 测试与文档（1-2 天）

**任务：**
- 功能测试
- 跨平台测试
- 编写用户文档

**总工作量：** 13-20 天

---

## 九、风险与缓解

| 风险 | 级别 | 缓解措施 |
|------|------|----------|
| Electron Chromium 不兼容 Remotion | 中 | 首次渲染失败时自动下载 Remotion 的 Chromium |
| ffmpeg 打包路径问题 | 中 | 提供 fallback 到系统 ffmpeg |
| 场景过多导致预览卡顿 | 低 | 添加虚拟滚动 + 预览降帧 |
| 数据库迁移失败 | 低 | 迁移前自动备份 creator.db |

---

## 十、成功标准

1. **用户体验**
   - 启动应用即可使用视频功能，无需额外配置
   - Agent 对话能完整创建视频项目
   - UI 能完整管理项目、场景、素材
   - 导出视频成功率 > 95%

2. **性能**
   - 应用启动时间不增加
   - 预览延迟 < 500ms
   - 渲染速度与 Remotion CLI 相当

3. **稳定性**
   - 无内存泄漏
   - 跨平台兼容
   - 错误处理完善

---

## 十一、后续增强（V2）

- 音频轨道：背景音乐、配音、音效
- 字幕系统：@remotion/captions 集成
- 视频嵌入：在场景中嵌入用户上传的视频片段
- AI 配音：ElevenLabs TTS 集成
- 更多过渡效果：自定义过渡动画
- 运行时编译：esbuild 编译 agent 写的 TSX（高级用户）
