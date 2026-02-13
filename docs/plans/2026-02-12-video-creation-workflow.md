# 视频创作流程集成实施计划

## Context

### 问题背景
当前 Sprouty AI 的内容创作功能只有一个统一的"创建内容"按钮，创建后所有内容都进入相同的流程（`researching → scripting → creating → adapting → scheduled → published`）。但实际上：

1. **一个内容可以同时创作为图文和视频**：同一个选题可以产出图文版本和视频版本
2. **图文和视频的产出分别记录在 `content_stages` 表中**：通过不同的 `stage` 值区分（如 `draft_image_text` 和 `draft_video`）
3. **目录结构按内容类型分离**：
   - 图文：`工作区/项目名/图文/日期/序号_标题/`
   - 视频：`工作区/项目名/视频/日期/序号_标题/`

### 设计目标
1. 将"创建内容"按钮拆分为"图文创作"和"视频创作"两个独立入口
2. 支持一个内容同时创作图文和视频版本
3. 通过 `content_stages` 表的 `stage` 字段区分不同类型的产出
4. 视频工作台与内容表串联，支持从内容列表直接跳转到视频编辑
5. 完善视频创作的状态流转和数据管理

### 关键约束
- **删除 `contents.content_type` 字段**：内容类型信息完全由 `content_stages` 管理
- **目录结构遵循现有规范**：`工作区/项目名/视频/日期/序号_标题/`
- **stage 字段区分类型**：使用 `draft_image_text`、`draft_video` 等明确的 stage 值
- 视频脚本技能已存在（`/skill-marketplace/apps/app-creator-media/skills/video-script-create`）
- Video MCP 服务器已实现（`packages/video`），提供 18 个视频创作工具
- 视频工作台已实现（`VideoStudio.tsx`），需要集成到内容流程中

---

## 实施方案

### 核心设计决策

#### 1. 数据库架构调整

**决策 1：删除 `contents.content_type` 字段**

理由：
- 一个内容可以同时创作为图文和视频
- 内容类型信息完全由 `content_stages` 表管理
- 避免数据冗余和不一致

**决策 2：通过 `stage` 字段区分内容类型**

`content_stages` 表的 `stage` 枚举值：

| stage 值 | 说明 | 文件路径示例 |
|---------|------|------------|
| `topic_recommend` | 选题推荐 | `项目名/选题推荐/2026-02-12/01_标题_85.md` |
| `research` | 灵感调研 | `项目名/选题推荐/2026-02-12/01_标题_调研.md` |
| `script_image_text` | 图文脚本 | `项目名/创作脚本/2026-02-12/01_标题_图文脚本.md` |
| `script_video` | 视频脚本 | `项目名/创作脚本/2026-02-12/01_标题_视频脚本.md` |
| `draft_image_text` | 图文原稿 | `项目名/图文/2026-02-12/01_标题/原稿.md` |
| `draft_video` | 视频项目 | `项目名/视频/2026-02-12/01_标题/project.json` |
| `platform_adapt_image_text` | 图文平台适配 | `项目名/图文/2026-02-12/01_标题/平台适配/小红书.md` |
| `platform_adapt_video` | 视频平台适配 | `项目名/视频/2026-02-12/01_标题/平台适配/抖音.md` |

**数据库迁移 SQL：**

```sql
-- 迁移版本 12：删除 content_type 字段
ALTER TABLE contents DROP COLUMN content_type;
```

#### 2. 目录结构规范

**图文创作目录：**
```
工作区/
└── 项目名/
    └── 图文/
        └── 2026-02-12/
            └── 01_标题/
                ├── 原稿.md
                ├── assets/
                │   ├── cover.jpg
                │   └── img_01.jpg
                └── 平台适配/
                    ├── 小红书.md
                    ├── 抖音.md
                    └── B站.md
```

**视频创作目录：**
```
工作区/
└── 项目名/
    └── 视频/
        └── 2026-02-12/
            └── 01_标题/
                ├── 视频脚本.md
                ├── project.json      # VideoProject 配置
                ├── assets/           # 素材文件
                │   ├── images/
                │   ├── videos/
                │   └── audio/
                ├── output/           # 渲染输出
                │   └── final.mp4
                └── 平台适配/
                    ├── 抖音.md
                    └── B站.md
```

**共享目录：**
```
工作区/
└── 项目名/
    ├── 选题推荐/
    │   └── 2026-02-12/
    │       └── 01_标题_85.md
    └── 创作脚本/
        └── 2026-02-12/
            ├── 01_标题_图文脚本.md
            └── 01_标题_视频脚本.md
```

#### 3. UI 改造
**决策：拆分为两个并列按钮**

```
[图文创作] [视频创作]
```

优势：
- 视频创作是核心功能，应该一级入口
- 减少点击次数
- 视觉上更清晰

#### 4. 内容状态流转

**状态保持不变**：`researching → scripting → creating → adapting → scheduled → published → archived`

**关键点**：
- 状态是内容级别的，不区分图文或视频
- 一个内容可以同时处于多个创作阶段（图文 creating + 视频 scripting）
- 通过查询 `content_stages` 表判断具体哪些类型的产出已完成

**判断逻辑示例：**
```typescript
// 判断是否有图文产出
const hasImageText = await db.query(
  'SELECT 1 FROM content_stages WHERE content_id = ? AND stage LIKE "%_image_text" LIMIT 1',
  [contentId]
)

// 判断是否有视频产出
const hasVideo = await db.query(
  'SELECT 1 FROM content_stages WHERE content_id = ? AND stage LIKE "%_video" LIMIT 1',
  [contentId]
)
```

---

## 详细实施步骤

### 阶段 1：数据库迁移（优先级：最高）

#### 1.1 创建数据库迁移脚本

**文件：** `packages/shared/src/db/migrations.ts`

**改动点：**
添加版本 12 迁移，删除 `contents.content_type` 字段

**具体实现：**

```typescript
// 在 MIGRATIONS 数组中添加
{
  version: 12,
  description: '删除 contents.content_type 字段，类型信息由 content_stages 管理',
  up: (db: Database) => {
    // SQLite 不支持 DROP COLUMN，需要重建表
    db.exec(`
      -- 1. 创建新表（无 content_type 字段）
      CREATE TABLE contents_new (
        id                TEXT PRIMARY KEY,
        project_id        TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
        title             TEXT,
        status            TEXT NOT NULL DEFAULT 'idea',
        target_platforms  TEXT,
        pipeline_mode     TEXT DEFAULT 'semi-auto',
        content_dir_path  TEXT,
        viral_pattern_id  TEXT,
        metadata          TEXT,
        created_at        DATETIME DEFAULT CURRENT_TIMESTAMP,
        updated_at        DATETIME DEFAULT CURRENT_TIMESTAMP
      );

      -- 2. 复制数据（排除 content_type）
      INSERT INTO contents_new
      SELECT id, project_id, title, status, target_platforms,
             pipeline_mode, content_dir_path, viral_pattern_id,
             metadata, created_at, updated_at
      FROM contents;

      -- 3. 删除旧表
      DROP TABLE contents;

      -- 4. 重命名新表
      ALTER TABLE contents_new RENAME TO contents;

      -- 5. 重建索引
      CREATE INDEX IF NOT EXISTS idx_content_project ON contents(project_id);
      CREATE INDEX IF NOT EXISTS idx_content_status ON contents(status);
    `)
  }
}
```

#### 1.2 更新类型定义

**文件：** `packages/shared/src/db/types.ts`

**改动点：**
1. 删除 `ContentType` 类型导出
2. 从 `Content` 接口中删除 `content_type` 字段
3. 添加 `ContentStageType` 枚举

**具体实现：**

```typescript
// 删除这些行：
export type ContentType =
  | 'video'
  | 'image-text'
  | 'article'
  | 'short-video'
  | 'live';

// 从 Content 接口中删除：
interface Content {
  // content_type: ContentType  // 删除此行
  ...
}

// 新增 stage 类型枚举
export type ContentStageType =
  | 'topic_recommend'
  | 'research'
  | 'script_image_text'
  | 'script_video'
  | 'draft_image_text'
  | 'draft_video'
  | 'platform_adapt_image_text'
  | 'platform_adapt_video';

// 更新 ContentStageRecord 接口
export interface ContentStageRecord {
  id: string
  content_id: string
  stage: ContentStageType  // 使用新的类型
  file_path: string
  status: ContentStageStatus
  version: number
  source_type: ContentStageSourceType
  metadata: string | null
  created_at: string
  updated_at: string
}
```

### 阶段 2：UI 改造（优先级：高）

#### 2.1 简化创建内容对话框

**文件：** `apps/electron/src/renderer/pages/creator-media/components/CreateContentDialog.tsx`

**改动点：**
1. 删除内容类型选择器（不再需要）
2. 删除视频相关字段（视频模板、分辨率）
3. 简化为纯粹的内容创建表单

**具体实现：**

```typescript
// 删除这些状态：
const [contentType, setContentType] = useState<ContentType>('image-text')  // 删除
const [videoTemplateId, setVideoTemplateId] = useState<string | undefined>()  // 删除
const [videoAspectRatio, setVideoAspectRatio] = useState<string>('16:9')  // 删除
const isVideoType = contentType === 'video' || contentType === 'short-video'  // 删除

// 删除视频相关的 UI 部分（第 184-221 行）

// 简化 handleCreate 函数：
const handleCreate = async () => {
  if (!title.trim()) return
  setSaving(true)
  try {
    await onCreateContent({
      title: title.trim(),
      topic: topic.trim() || null,
      status: 'researching',
      target_platforms: selectedPlatforms.length > 0 ? JSON.stringify(selectedPlatforms) : null,
      pipeline_mode: 'manual',
      scheduled_at: scheduledAt || null,
      metadata: null,  // 不再需要视频元数据
    })
    handleClose(false)
  } finally {
    setSaving(false)
  }
}
```

#### 2.2 拆分创建按钮

**文件：** `apps/electron/src/renderer/pages/creator-media/ProjectDashboard.tsx`

**改动点：**
将单个"新建内容"按钮拆分为"图文创作"和"视频创作"两个按钮

**具体实现：**

```typescript
// 原代码（第 401-410 行）：
<button onClick={() => setShowCreateContent(true)}>
  新建内容
</button>

// 改为：
const [showCreateContent, setShowCreateContent] = useState(false)
const [creationType, setCreationType] = useState<'image-text' | 'video'>('image-text')

<div className="flex items-center gap-2">
  <button
    onClick={() => {
      setCreationType('image-text')
      setShowCreateContent(true)
    }}
    className="..."
  >
    <FileText className="w-3.5 h-3.5" />
    图文创作
  </button>

  <button
    onClick={() => {
      setCreationType('video')
      setShowCreateContent(true)
    }}
    className="..."
  >
    <Video className="w-3.5 h-3.5" />
    视频创作
  </button>
</div>

// 创建内容后，根据 creationType 自动触发对应的技能
const handleCreateContent = async (data) => {
  const content = await window.electronAPI.creatorMedia.contents.create(workspace.id, data)

  // 根据创建类型自动触发对应的脚本创作技能
  if (creationType === 'image-text') {
    // 触发图文脚本创作
    navigate(routes.view.appView('creator-media', 'chat', {
      input: `[skill:${workspace.id}:script-create-image-text] 为内容「${data.title}」创建图文脚本。\n\n内容 ID: ${content.id}`,
      send: true,
    }))
  } else if (creationType === 'video') {
    // 触发视频脚本创作
    navigate(routes.view.appView('creator-media', 'chat', {
      input: `[skill:${workspace.id}:video-script-create] 为内容「${data.title}」创建视频脚本。\n\n内容 ID: ${content.id}`,
      send: true,
    }))
  }
}
```

#### 2.3 增强 ContentTable 显示逻辑

**文件：** `apps/electron/src/renderer/pages/creator-media/components/ContentTable.tsx`

**改动点：**
1. 删除基于 `content_type` 的判断逻辑
2. 通过查询 `content_stages` 判断内容有哪些类型的产出
3. 根据产出类型显示对应的操作按钮

**具体实现：**

```typescript
// 新增辅助函数：查询内容的产出类型
const getContentOutputTypes = async (contentId: string): Promise<{
  hasImageText: boolean
  hasVideo: boolean
}> => {
  const stages = await window.electronAPI.creatorMedia.contentStages.list(workspace.id, contentId)

  return {
    hasImageText: stages.some(s => s.stage.includes('_image_text')),
    hasVideo: stages.some(s => s.stage.includes('_video')),
  }
}

// 在渲染时使用
{item.status === 'creating' && (
  <>
    {/* 如果有图文产出，显示图文编辑按钮 */}
    {outputTypes[item.id]?.hasImageText && (
      <button onClick={() => onEditImageText?.(item)}>
        <FileText className="h-3.5 w-3.5" />
        {t('编辑图文')}
      </button>
    )}

    {/* 如果有视频产出，显示视频工作台按钮 */}
    {outputTypes[item.id]?.hasVideo && (
      <button onClick={() => onOpenVideoStudio?.(item)}>
        <Video className="h-3.5 w-3.5" />
        {t('视频工作台')}
      </button>
    )}
  </>
)}
```

#### 2.4 视频工作台路由集成

**文件：** `apps/electron/src/renderer/pages/creator-media/VideoStudio.tsx`

**改动点：**
1. 支持通过 URL 参数 `contentId` 加载关联内容
2. 修改"完成制作"按钮的状态更新逻辑（从 `reviewing` 改为 `adapting`）

**具体实现：**

```typescript
// 1. 接收 contentId 参数
const searchParams = new URLSearchParams(window.location.search)
const contentId = searchParams.get('contentId')

useEffect(() => {
  if (contentId) {
    const content = videoContents.find(c => c.id === contentId)
    if (content) {
      handleContentSelect(content)
    }
  }
}, [contentId, videoContents])

// 2. 修改完成制作逻辑（第 287-298 行）
const handleFinishCreation = useCallback(async () => {
  if (!activeContentId) return
  try {
    // 改为 adapting 状态，而不是 reviewing
    await window.electronAPI.creatorMedia.contents.updateStatus(
      workspace.id,
      activeContentId,
      'adapting'
    )
    toast.success(t('视频制作完成，进入平台适配阶段'))
    navigate(routes.view.appView('creator-media', 'dashboard'))
  } catch {
    toast.error(t('状态更新失败'))
  }
}, [activeContentId, workspace.id, navigate, t])
```

在 `ProjectDashboard.tsx` 中添加跳转逻辑：
```typescript
const handleOpenVideoStudio = useCallback((content: Content) => {
  navigate(routes.view.appView('creator-media', 'video-studio', {
    contentId: content.id
  }))
}, [navigate])

// 传递给 ContentTable
<ContentTable
  items={contents}
  onOpenVideoStudio={handleOpenVideoStudio}
  // ... 其他 props
/>
```

---

### 阶段 2：技能集成（优先级：高）

#### 2.1 视频脚本创作技能调用

**文件：** `apps/electron/src/renderer/pages/creator-media/ProjectDashboard.tsx`

**改动点：**
添加"生成脚本"按钮的处理函数，触发视频脚本创作技能

**具体实现：**

```typescript
const handleGenerateScript = useCallback((content: Content) => {
  // 构建技能调用提示词
  const prompt = `为内容「${content.title}」生成视频脚本。

内容 ID: ${content.id}
内容类型: ${content.content_type}
选题描述: ${content.topic || '无'}
目标平台: ${content.target_platforms ? JSON.parse(content.target_platforms).join('、') : '未指定'}

请根据选题生成结构化的视频脚本，包含分镜描述、旁白、字幕和时长规划。`

  // 导航到新 session 并激活技能
  navigate(routes.view.appView('creator-media', 'chat', {
    input: `[skill:${workspace.id}:video-script-create] ${prompt}`,
    send: true,  // 自动发送
  }))
}, [workspace.id, navigate])

// 传递给 ContentTable
<ContentTable
  items={contents}
  onGenerateScript={handleGenerateScript}
  // ... 其他 props
/>
```

#### 2.2 创建视频制作技能

**位置：** 需要在 `skill-marketplace/apps/app-creator-media/skills/` 下创建新技能

**技能 ID：** `video-production`

**技能文件结构：**
```
skill-marketplace/apps/app-creator-media/skills/video-production/
├── SKILL.md           # 技能提示词和工作流程
├── skill.json         # 技能元数据配置
└── examples/          # 示例脚本和输出
```

**SKILL.md 核心内容：**

```markdown
---
name: 视频制作
description: 基于视频脚本调用 video-mcp 完成素材准备、项目创建和视频渲染
icon: 🎬
---

# 视频制作

你是一个专业的视频制作助手。你的职责是根据视频脚本，调用 video-mcp 工具完成完整的视频制作流程。

## 工作流程

### 步骤一：读取视频脚本

从 `content_stages` 表读取最新的视频脚本：

```bash
sqlite3 .sprouty-ai/db/creator.db "
  SELECT file_path FROM content_stages
  WHERE content_id = '[content_id]' AND stage = 'script'
  ORDER BY version DESC LIMIT 1
"
```

读取脚本文件内容（Markdown 格式）。

### 步骤二：解析脚本结构

从脚本中提取：
- 视频基本信息（标题、时长、分辨率）
- 分镜列表（每个分镜的画面描述、旁白、字幕、时长）
- 素材需求（图片、视频、音频）
- 转场效果

### 步骤三：创建视频项目

调用 video-mcp 的 `video_create_project` 工具：

```typescript
const project = await mcp.call('video_create_project', {
  name: content.title,
  width: metadata.width,
  height: metadata.height,
  fps: 30,
  output_dir: `{content_dir_path}/视频创作/output/`
})
```

### 步骤四：素材准备

按优先级匹配素材：
1. 用户上传素材（检查 `{content_dir_path}/assets/` 目录）
2. AI 生成（调用 `baoyu-image-gen` 技能）
3. 免费素材库（Pexels/Unsplash API）
4. 标记需要用户录制

调用 `video_add_asset` 添加素材到项目：

```typescript
await mcp.call('video_add_asset', {
  project_id: project.id,
  asset_path: assetPath,
  asset_type: 'image' | 'video' | 'audio'
})
```

### 步骤五：创建 Remotion 组合

根据分镜脚本生成 Remotion composition.json：

```json
{
  "id": "video-{content_id}",
  "fps": 30,
  "width": 1080,
  "height": 1920,
  "durationInFrames": 1800,
  "scenes": [
    {
      "sceneIndex": 0,
      "from": 0,
      "durationInFrames": 90,
      "componentType": "TextOverlay",
      "props": { "text": "...", "fontSize": 48 }
    }
  ]
}
```

调用 `video_add_composition` 添加到项目：

```typescript
await mcp.call('video_add_composition', {
  project_id: project.id,
  composition_id: `video-${content.id}`,
  composition_data: compositionJson
})
```

### 步骤六：渲染视频

调用 `video_render` 开始渲染：

```typescript
const renderResult = await mcp.call('video_render', {
  project_id: project.id,
  composition_id: `video-${content.id}`,
  quality: 'high',
  format: 'mp4'
})
```

### 步骤七：更新内容元数据

渲染完成后，更新 `contents.metadata`：

```bash
sqlite3 .sprouty-ai/db/creator.db "
  UPDATE contents
  SET metadata = json_set(
    COALESCE(metadata, '{}'),
    '$.video_project_id', '${project.id}',
    '$.video_project_name', '${project.name}',
    '$.video_render_status', 'completed',
    '$.video_output_path', '${renderResult.output_path}',
    '$.video_duration', ${renderResult.duration}
  ),
  status = 'adapting',
  updated_at = datetime('now')
  WHERE id = '[content_id]'
"
```

## 注意事项

- 渲染过程可能需要数分钟，需要显示进度提示
- 如果素材不足，标记缺失素材并提示用户补充
- 渲染失败时，更新 `video_render_status` 为 'failed' 并记录错误信息
```

**skill.json 配置：**

```json
{
  "id": "video-production",
  "name": "视频制作",
  "description": "基于视频脚本调用 video-mcp 完成素材准备、项目创建和视频渲染",
  "version": "1.0.0",
  "author": "Sprouty AI",
  "category": "video-production",
  "tags": ["video", "production", "remotion", "mcp"],
  "mcp_servers": ["video"],
  "required_tools": [
    "video_create_project",
    "video_add_asset",
    "video_add_composition",
    "video_render"
  ]
}
```

---

### 阶段 3：数据流优化（优先级：中）

#### 3.1 确保 IPC 通道完整

**文件：** `apps/electron/src/main/creator-media-ipc.ts`

**验证点：**
- ✅ `contents.list` - 已实现
- ✅ `contents.get` - 已实现
- ✅ `contents.create` - 已实现
- ✅ `contents.update` - 已实现
- ✅ `contents.updateStatus` - 已实现
- ✅ `contents.delete` - 已实现

**无需改动**，现有 IPC 通道已完整支持所有操作。

#### 3.2 元数据辅助函数

**文件：** `packages/shared/src/db/repositories/contents.ts`

**验证点：**
- ✅ `parseContentMetadata()` - 已实现（第 23-35 行）
- ✅ `updateContentVideoMetadata()` - 已实现（第 37-58 行）

**无需改动**，现有辅助函数已支持视频元数据的解析和更新。

---

### 阶段 4：文件系统组织（优先级：中）

#### 4.1 内容目录结构

**标准结构：**
```
{Workspace}/
└── 创作媒体/
    └── {项目名}/
        └── 内容/
            └── {内容ID}/
                ├── 选题分析.md          # researching 阶段
                ├── 视频脚本.md          # scripting 阶段（视频专属）
                ├── 视频创作/            # creating 阶段（视频专属）
                │   ├── project.json    # VideoProject 配置
                │   ├── assets/         # 素材文件
                │   │   ├── images/
                │   │   ├── videos/
                │   │   └── audio/
                │   └── output/         # 渲染输出
                │       └── final.mp4
                └── 平台适配/            # adapting 阶段
                    ├── 小红书.md
                    ├── 抖音.md
                    └── B站.md
```

**实现方式：**
- 在视频脚本创作技能中创建 `视频脚本.md`
- 在视频制作技能中创建 `视频创作/` 目录和子目录
- 通过 `content_dir_path` 字段关联内容和文件系统

#### 4.2 content_stages 记录

**视频脚本阶段：**
```typescript
{
  content_id: contentId,
  stage: 'script',
  file_path: '{content_dir_path}/视频脚本.md',
  status: 'completed',
  version: 1,
  source_type: 'agent',
  metadata: JSON.stringify({
    word_count: 1200,
    scene_count: 8,
    duration_total: 60
  })
}
```

**视频项目阶段：**
```typescript
{
  content_id: contentId,
  stage: 'draft',  // 复用 draft 阶段
  file_path: '{content_dir_path}/视频创作/project.json',
  status: 'completed',
  version: 1,
  source_type: 'agent',
  metadata: JSON.stringify({
    video_project_id: projectId,
    render_status: 'completed'
  })
}
```

---

## 关键文件清单

### 必须修改的文件（优先级：最高）

1. **`packages/shared/src/db/migrations.ts`**
   - 添加版本 12 迁移：删除 `contents.content_type` 字段
   - 重建 contents 表结构

2. **`packages/shared/src/db/types.ts`**
   - 删除 `ContentType` 类型导出
   - 从 `Content` 接口删除 `content_type` 字段
   - 添加 `ContentStageType` 枚举（包含 `script_video`, `draft_video` 等）

3. **`apps/electron/src/renderer/pages/creator-media/components/CreateContentDialog.tsx`**
   - 删除内容类型选择器
   - 删除视频相关字段（模板、分辨率）
   - 简化为纯粹的内容创建表单

4. **`apps/electron/src/renderer/pages/creator-media/ProjectDashboard.tsx`**
   - 拆分"新建内容"按钮为"图文创作"和"视频创作"
   - 创建后自动触发对应的脚本创作技能

5. **`apps/electron/src/renderer/pages/creator-media/components/ContentTable.tsx`**
   - 删除基于 `content_type` 的判断逻辑
   - 通过查询 `content_stages` 判断产出类型
   - 根据产出类型显示对应操作按钮

### 需要增强的文件（优先级：高）

6. **`apps/electron/src/main/creator-media-ipc.ts`**
   - 添加 `content_stages` 相关 IPC 通道
   - `list`, `getLatest`, `create` 方法

7. **`apps/electron/src/shared/types.ts`**
   - 添加 `content_stages` IPC 通道常量
   - 扩展 `ElectronAPI` 接口

8. **`apps/electron/src/preload/index.ts`**
   - 绑定 `content_stages` IPC 通道

9. **`apps/electron/src/renderer/pages/creator-media/VideoStudio.tsx`**
   - 支持 `contentId` URL 参数
   - 修改"完成制作"状态更新为 `adapting`

### 需要创建的文件（优先级：高）

10. **`skill-marketplace/apps/app-creator-media/skills/video-production/SKILL.md`**
    - 视频制作技能提示词
    - 完整工作流程定义

11. **`skill-marketplace/apps/app-creator-media/skills/video-production/skill.json`**
    - 技能元数据配置

### 需要适配的文件（优先级：中）

12. **`skill-marketplace/apps/app-creator-media/skills/video-script-create/SKILL.md`**
    - 更新脚本保存路径为新目录结构
    - 使用 `stage = 'script_video'`
    - 更新状态为 `creating`

---

### 必须修改的文件（优先级：高）

1. **`apps/electron/src/renderer/pages/creator-media/ProjectDashboard.tsx`**
   - 拆分创建按钮为"图文创作"和"视频创作"
   - 添加 `handleGenerateScript` 和 `handleOpenVideoStudio` 函数
   - 传递回调给 `ContentTable`

2. **`apps/electron/src/renderer/pages/creator-media/components/CreateContentDialog.tsx`**
   - 添加 `presetContentType` prop
   - 根据预设类型初始化表单

3. **`apps/electron/src/renderer/pages/creator-media/components/ContentTable.tsx`**
   - 添加"生成脚本"按钮（`scripting` 状态 + 视频类型）
   - 添加 `onGenerateScript` 回调 prop

4. **`apps/electron/src/renderer/pages/creator-media/VideoStudio.tsx`**
   - 支持 `contentId` URL 参数
   - 修改"完成制作"按钮状态更新为 `adapting`

### 需要创建的文件（优先级：高）

5. **`skill-marketplace/apps/app-creator-media/skills/video-production/SKILL.md`**
   - 视频制作技能提示词
   - 完整的工作流程定义

6. **`skill-marketplace/apps/app-creator-media/skills/video-production/skill.json`**
   - 技能元数据配置

### 可选优化的文件（优先级：低）

7. **`packages/shared/src/db/repositories/content-stages.ts`**
   - 添加视频脚本和视频项目的辅助查询函数（可选）

---

## 测试验证步骤

### 1. 图文创作流程测试

```
1. 点击"图文创作"按钮
2. 验证对话框打开，content_type 预设为 'image-text'
3. 填写标题、选题、平台
4. 创建后验证：
   - content_type = 'image-text'
   - status = 'researching'
   - metadata = null
5. 继续现有流程...
```

### 2. 视频创作流程测试

```
1. 点击"视频创作"按钮
2. 验证对话框打开，content_type 预设为 'video'
3. 填写标题、选题、选择分辨率、模板
4. 创建后验证：
   - content_type = 'video'
   - status = 'researching'
   - metadata 包含 videoTemplateId, aspectRatio, width, height
5. 执行灵感调研，状态转为 'scripting'
6. 点击"生成脚本"按钮
7. 验证：
   - 跳转到聊天界面
   - 技能自动激活：[skill:ws:video-script-create]
   - 提示词包含内容 ID 和选题信息
8. 脚本生成完成后，验证：
   - 脚本文件已创建：{content_dir_path}/视频脚本.md
   - content_stages 记录已创建
   - status 转为 'creating'
9. 点击"视频工作台"按钮
10. 验证：
    - 正确跳转到 VideoStudio 视图
    - contentId 参数正确传递
    - 关联的内容自动加载
11. 在视频工作台中编辑和渲染
12. 点击"完成制作"
13. 验证：
    - status 转为 'adapting'
    - metadata.video_render_status = 'completed'
    - 返回内容列表
```

### 3. 视频工作台集成测试

```
1. 从内容表点击"视频工作台"按钮（creating 状态）
2. 验证：
   - 正确跳转到 VideoStudio 视图
   - URL 包含 contentId 参数
   - 关联的内容自动加载到左侧面板
3. 编辑视频项目（添加素材、调整属性）
4. 导出渲染
5. 点击"完成制作"
6. 验证：
   - 返回内容列表
   - 状态已更新为 'adapting'
   - 可以继续平台适配流程
```

---

## 风险与注意事项

### 技术风险

1. **video-mcp 工具可用性**
   - 风险：video-mcp 服务器可能未完全实现或存在 bug
   - 缓解：先在视频工作台手动测试 video-mcp 工具，确保基本功能可用

2. **技能系统集成**
   - 风险：技能调用格式或注册机制可能与预期不符
   - 缓解：参考现有技能（如 `video-script-create`）的实现方式

3. **文件路径管理**
   - 风险：跨平台路径兼容性问题
   - 缓解：使用 Node.js `path` 模块处理路径，确保跨平台兼容

### 用户体验风险

1. **视频渲染时间**
   - 风险：长视频渲染可能需要数分钟，用户体验差
   - 缓解：显示渲染进度条，支持后台渲染

2. **素材管理**
   - 风险：大量素材文件占用磁盘空间
   - 缓解：提供素材清理功能，定期清理未使用的素材

3. **错误处理**
   - 风险：脚本生成失败、素材下载失败、渲染失败
   - 缓解：友好的错误提示和重试机制

### 建议

1. **分阶段实施**
   - 先实现 UI 改造和基本流程（阶段 1）
   - 再实现技能集成（阶段 2）
   - 最后优化数据流和文件系统（阶段 3-4）

2. **保持向后兼容**
   - 现有图文创作流程不受影响
   - 旧数据可以正常访问

3. **文档更新**
   - 更新 CLAUDE.md 记录新的视频创作流程
   - 添加技能使用文档

---

## 实施时间估算

- **阶段 1（UI 改造）**：1-2 天
- **阶段 2（技能集成）**：2-3 天
- **阶段 3（数据流优化）**：1 天
- **阶段 4（集成测试）**：1-2 天

**总计：5-8 天**

---

## 成功标准

1. ✅ 用户可以通过"图文创作"和"视频创作"两个独立按钮创建不同类型的内容
2. ✅ 视频创作流程包含视频脚本生成和视频制作环节
3. ✅ 视频工作台可以从内容列表直接跳转，并正确加载关联内容
4. ✅ 视频渲染完成后，状态正确更新为 `adapting`
5. ✅ 所有数据正确存储在数据库和文件系统中
6. ✅ 图文创作流程不受影响，保持向后兼容
