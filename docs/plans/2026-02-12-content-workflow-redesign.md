# 内容创作工作流重构设计

> 日期：2026-02-12
> 状态：待实施

## 1. 设计目标

- 每一步创作产出保存为 Markdown 文件，Agent 可直接读写
- 数据库只记录文件路径、状态和元数据，不存储内容本身
- 工作区支持多项目，内容目录结构清晰规范
- 精简 `contents` 表，新增 `content_stages` 表记录每个阶段的产出
- 废弃 `content_versions` 表，由 `content_stages` 的 version 字段替代
- 全流程 Skill 驱动，每个阶段由独立 Skill 创建 Agent Session 完成

---

## 2. 内容目录结构

内容文件存放在工作区根目录下（非隐藏目录），按项目组织：

```
工作区根目录/
├── .sprouty-ai/                              # 系统数据（隐藏）
│   ├── db/creator.db
│   └── ...
├── 项目A/
│   ├── 选题推荐/
│   │   └── 2026-02-12/
│   │       ├── 01_AI短视频趋势_85.md
│   │       └── 02_春季穿搭指南_72.md
│   ├── 创作脚本/
│   │   └── 2026-02-12/
│   │       └── 01_AI短视频趋势_创作脚本.md
│   ├── 图文/
│   │   └── 2026-02-12/
│   │       └── 01_AI短视频趋势/
│   │           ├── 原稿.md
│   │           ├── assets/
│   │           │   ├── cover.jpg
│   │           │   └── img_01.jpg
│   │           └── 平台适配/
│   │               ├── 小红书.md
│   │               └── 公众号.md
│   └── 视频/
│       └── 2026-02-12/
│           └── 01_AI短视频趋势/
│               ├── 原稿.md
│               ├── assets/
│               │   ├── bgm.mp3
│               │   └── footage_01.mp4
│               └── 平台适配/
│                   ├── 抖音.md
│                   └── B站.md
└── 项目B/
    └── ...
```

### 文件命名规范

| 类型 | 格式 | 示例 |
|------|------|------|
| 选题推荐 | `{序号}_{标题}_{推荐指数}.md` | `01_AI短视频趋势_85.md` |
| 创作脚本 | `{序号}_{标题}_创作脚本.md` | `01_AI短视频趋势_创作脚本.md` |
| 内容文件夹 | `{序号}_{内容标题}/` | `01_AI短视频趋势/` |
| 原稿 | `原稿.md` | — |
| 平台适配 | `{平台名称}.md` | `小红书.md` |

- 序号：当日自增，两位数字补零
- 标题：取内容标题，过长时截断（建议 ≤ 20 字符）
- 日期格式：`YYYY-MM-DD`

### 路径约定

- 所有 `file_path` 和 `content_dir_path` 均为相对于工作区根目录的相对路径
- 示例：`项目A/图文/2026-02-12/01_AI短视频趋势/原稿.md`

---

## 3. 数据库变更

### 3.1 精简 contents 表

**保留字段：**

```sql
CREATE TABLE contents (
  id                TEXT PRIMARY KEY,
  project_id        TEXT NOT NULL,
  title             TEXT,
  content_type      TEXT,              -- video | short-video | image-text | article
  status            TEXT NOT NULL,     -- 见状态定义
  target_platforms  TEXT,              -- JSON: ["小红书", "抖音"]
  pipeline_mode     TEXT DEFAULT 'semi-auto',
  content_dir_path  TEXT,              -- 内容文件夹相对路径
  viral_pattern_id  TEXT,
  metadata          TEXT,              -- JSON: 扩展元数据
  created_at        TEXT NOT NULL,
  updated_at        TEXT NOT NULL
);
```

**删除字段：**

| 字段 | 原用途 | 迁移去向 |
|------|--------|---------|
| `topic` | 选题内容 | content_stages(topic_recommend) 文件 |
| `topic_source` | 选题来源 | content_stages metadata |
| `source_topic_id` | 关联选题 ID | content_stages metadata |
| `script_path` | 脚本路径 | content_stages(script).file_path |
| `pipeline_state` | 流水线状态 | 通过查询 content_stages 推导 |
| `files` | 文件列表 | content_dir_path + 文件系统 |
| `review_summary` | 复盘摘要 | metadata |

### 3.2 新增 content_stages 表

```sql
CREATE TABLE content_stages (
  id            TEXT PRIMARY KEY,
  content_id    TEXT NOT NULL,
  stage         TEXT NOT NULL,        -- 阶段类型
  file_path     TEXT NOT NULL,        -- 产出文件相对路径
  status        TEXT DEFAULT 'draft', -- draft | completed | revised
  version       INTEGER DEFAULT 1,   -- 同阶段多版本
  source_type   TEXT,                 -- auto | user_edit | agent
  metadata      TEXT,                 -- JSON: 阶段特有元数据
  created_at    TEXT NOT NULL,
  updated_at    TEXT NOT NULL,
  FOREIGN KEY (content_id) REFERENCES contents(id)
);

CREATE INDEX idx_content_stages_content_id ON content_stages(content_id);
CREATE INDEX idx_content_stages_stage ON content_stages(stage);
```

**stage 枚举值：**

| stage | 触发时机 | file_path 示例 | metadata 示例 |
|-------|---------|---------------|--------------|
| `topic_recommend` | 采纳选题时 | `项目A/选题推荐/2026-02-12/01_标题_85.md` | `{ "heat_index": 85, "source_topic_id": "xxx" }` |
| `research` | 灵感调研完成 | `项目A/选题推荐/2026-02-12/01_标题_调研.md` | `{ "source": "idea" }` |
| `script` | 脚本创作完成 | `项目A/创作脚本/2026-02-12/01_标题_创作脚本.md` | `{ "word_count": 1200 }` |
| `draft` | 原稿创作完成 | `项目A/图文/2026-02-12/01_标题/原稿.md` | `{ "content_type": "image-text" }` |
| `platform_adapt` | 平台适配完成 | `项目A/图文/2026-02-12/01_标题/平台适配/小红书.md` | `{ "platform": "xiaohongshu", "char_count": 800 }` |

### 3.3 废弃 content_versions 表

`content_stages` 的 `version` 字段替代版本管理功能：
- 同一 `content_id` + `stage` 可有多条记录，`version` 自增
- 查询最新版本：`ORDER BY version DESC LIMIT 1`
- 历史版本保留，文件不删除

### 3.4 recommended_topics 表（保留）

作为选题池继续使用，采纳时的关联逻辑调整：

```
采纳选题:
  1. recommended_topics.status → 1（已采纳）
  2. 创建 contents 记录
  3. 创建 content_stages 记录（stage: topic_recommend）
     - file_path = recommended_topics.md_file_path
     - metadata = { "source_topic_id": recommended_topics.id }
```

---

## 4. 内容状态定义

```
idea → researching → scripting → creating → adapting → scheduled → published → archived
```

| 状态 | 含义 | 入口 |
|------|------|------|
| `idea` | 灵感记录，待研究 | 用户手动创建 |
| `researching` | Agent 调研中，生成选题分析 | 从 idea 推进 / Agent 触发 |
| `scripting` | 脚本创作中 | 从 researching 推进 / 选题池采纳 |
| `creating` | 内容创作中（原稿 + 素材） | 脚本完成后推进 |
| `adapting` | 平台适配中 | 原稿完成后推进 |
| `scheduled` | 已排期待发布 | 适配完成，加入发布队列 |
| `published` | 已发布 | 发布成功 |
| `archived` | 已归档 | 用户手动归档 |

### 两条入口路径

**路径 A — 选题池入口：**
```
recommended_topics 采纳
  → contents(status: scripting)
  → content_stages(topic_recommend)
  → scripting → creating → adapting → scheduled → published
```

**路径 B — 灵感入口：**
```
用户记录灵感
  → contents(status: idea)
  → idea → researching
  → content_stages(research)
  → scripting → creating → adapting → scheduled → published
```

---

## 5. Skill 驱动体系

每个创作阶段由独立 Skill 驱动，通过 Agent Session 执行。Skill 定义在 `skill-marketplace/skills/` 目录下。

### 5.1 Skill 总览

| Skill ID | 名称 | 驱动阶段 | 状态变更 | 现有情况 |
|----------|------|---------|---------|---------|
| `topic-generator` | 选题推荐 | 热榜 → 选题池 | — | 已有，需适配 |
| `idea-researcher` | 灵感调研 | idea → researching → scripting | idea → researching | 新建 |
| `script-create` | 脚本创作 | scripting → creating | scripting → creating | 已有，需适配 |
| `content-creator` | 内容创作 | creating → adapting | creating → adapting | 新建 |
| `platform-adapter` | 平台适配 | adapting → scheduled | adapting → scheduled | 新建 |

### 5.2 topic-generator 适配改动

现有 Skill 需要以下调整：

1. **目录结构适配**：文件路径从 `选题推荐/YYYY-MM-DD/` 改为 `{项目名称}/选题推荐/YYYY-MM-DD/`
2. **文件名分隔符**：从 `-` 改为 `_`（`01-标题-85.md` → `01_标题_85.md`）
3. **数据库操作**：`md_file_path` 路径格式更新，包含项目名称前缀

### 5.3 idea-researcher（新建）

**职责**：用户记录灵感后，Agent 调研相关热点和资料，生成选题分析报告。

**工作流程**：
1. 加载账号画像（同 topic-generator）
2. 读取灵感内容（从 contents 表获取 title + metadata）
3. 调研相关信息（Web 搜索 + 热榜匹配）
4. 生成选题分析报告（复用选题推荐的 Markdown 格式）
5. 写入文件：`{项目名称}/选题推荐/YYYY-MM-DD/{序号}_{标题}_{推荐指数}.md`
6. 创建 content_stages 记录（stage: research）
7. 更新 contents.status → researching
8. 输出摘要

**产出格式**：与 topic-generator 的 Markdown 格式一致，metadata 中标记 `"source": "idea"` 区分来源。

### 5.4 script-create 适配改动

现有 Skill 需要以下调整：

1. **目录结构适配**：文件路径从 `创作脚本/YYYY-MM-DD/` 改为 `{项目名称}/创作脚本/YYYY-MM-DD/`
2. **文件名分隔符**：从 `-` 改为 `_`
3. **数据库操作**：
   - 不再更新 `contents.script_path`（字段已删除）
   - 改为创建 `content_stages` 记录（stage: script）
   - 更新 `contents.status` → creating（脚本完成后进入内容创作阶段）

### 5.5 content-creator（新建）

**职责**：基于脚本生成完整原稿，创建内容文件夹和素材目录。内部根据 `content_type` 走不同模板。

**工作流程**：
1. 加载账号画像
2. 读取脚本文件（从 content_stages 查询最新 script 阶段的 file_path）
3. 确定内容类型（从 contents.content_type 获取）
4. 根据内容类型选择模板生成原稿
5. 创建内容文件夹和素材目录
6. 写入文件和更新数据库
7. 输出摘要

**图文模板产出**：
```
{项目名称}/图文/YYYY-MM-DD/{序号}_{内容标题}/
├── 原稿.md          # 完整图文文章
└── assets/          # 素材目录（配图建议清单）
```

原稿.md 内容结构：
```markdown
# [标题]

> 内容形式：图文 | 目标平台：[平台列表] | 创建时间：YYYY-MM-DD HH:MM

## 封面图建议

[封面图描述和风格建议]

## 正文

[完整图文内容，包含配图位置标记]

### [段落标题1]

[内容]

📷 配图建议：[图片描述]

### [段落标题2]

[内容]

📷 配图建议：[图片描述]

...

## 标签建议

[推荐标签列表]

## 创作备注

- **总字数**: [字数]
- **预估阅读时间**: [分钟]
- **配图数量**: [数量]
- **风格调性**: [调性说明]
```

**视频模板产出**：
```
{项目名称}/视频/YYYY-MM-DD/{序号}_{内容标题}/
├── 原稿.md          # 详细分镜脚本
└── assets/          # 素材目录（素材清单）
```

原稿.md 内容结构：
```markdown
# [标题]

> 内容形式：视频 | 目标平台：[平台列表] | 创建时间：YYYY-MM-DD HH:MM

## 视频概要

- **预估时长**: [时长]
- **视频比例**: [16:9 / 9:16 / 1:1]
- **风格**: [口播 / Vlog / 教程 / 混剪]

## 分镜脚本

### 场景 1：开场钩子（0:00 - 0:03）

| 项目 | 内容 |
|------|------|
| **画面** | [画面描述] |
| **台词** | [口播文案] |
| **字幕** | [字幕文案] |
| **BGM** | [音乐建议] |
| **转场** | [转场方式] |

### 场景 2：问题引入（0:03 - 0:13）

...

### 场景 N：行动号召

...

## 素材清单

| 序号 | 类型 | 描述 | 来源建议 |
|------|------|------|---------|
| 1 | 实拍 | [描述] | [建议] |
| 2 | 图片 | [描述] | [建议] |
| 3 | BGM | [描述] | [建议] |

## 创作备注

- **总时长**: [时长]
- **场景数**: [数量]
- **风格调性**: [调性说明]
- **注意事项**: [版权、敏感内容等]
```

**数据库操作**：
- 创建 content_stages 记录（stage: draft）
- 更新 contents.content_dir_path
- 更新 contents.status → adapting

### 5.6 platform-adapter（新建）

**职责**：基于原稿，按目标平台的内容规范生成适配版本。平台规范预置在 `references/` 目录中。

**工作流程**：
1. 读取原稿文件（从 content_stages 查询最新 draft 阶段的 file_path）
2. 获取目标平台列表（从 contents.target_platforms）
3. 加载平台规范（从 references/ 目录）
4. 为每个平台生成适配版本
5. 写入文件和更新数据库
6. 输出摘要

**平台规范 references/ 目录**：
```
platform-adapter/
├── SKILL.md
├── config.yaml
└── references/
    ├── 小红书.md        # 小红书内容规范
    ├── 抖音.md          # 抖音内容规范
    ├── 公众号.md        # 微信公众号内容规范
    ├── B站.md           # B站内容规范
    ├── 微博.md          # 微博内容规范
    ├── 知乎.md          # 知乎内容规范
    └── 快手.md          # 快手内容规范
```

每个平台规范文件包含：
- 内容字数/时长限制
- 标题规范（字数、关键词策略）
- 标签/话题策略
- 排版规范（段落长度、emoji 使用、分隔符）
- 封面图/缩略图规范
- 发布时间建议
- 平台算法偏好
- 违规红线

**适配版本 Markdown 格式**：
```markdown
# [适配后标题]

> 平台：[平台名称] | 适配时间：YYYY-MM-DD HH:MM | 原稿：[原稿路径]

## 发布内容

[按平台规范适配后的完整内容]

## 标签/话题

[平台专属标签列表]

## 封面图建议

[按平台规范的封面图描述]

## 发布建议

- **最佳发布时间**: [时间段]
- **字数**: [字数] / 限制 [上限]
- **注意事项**: [平台特有注意事项]
```

**数据库操作**：
- 每个平台创建一条 content_stages 记录（stage: platform_adapt）
- 所有平台适配完成后，更新 contents.status → scheduled

### 5.7 Skill 间的数据传递

Skill 之间通过数据库 + 文件系统松耦合传递数据，不依赖内存状态：

```
topic-generator
  → recommended_topics 表 + Markdown 文件
  → 用户采纳 → contents 表 + content_stages 表

idea-researcher
  → contents 表（读取灵感）
  → Markdown 文件 + content_stages 表

script-create
  → content_stages 表（读取 topic_recommend/research 阶段的 file_path）
  → 读取选题 Markdown 文件
  → 生成脚本 Markdown 文件 + content_stages 表

content-creator
  → content_stages 表（读取 script 阶段的 file_path）
  → 读取脚本 Markdown 文件
  → 生成原稿 + 内容文件夹 + content_stages 表

platform-adapter
  → content_stages 表（读取 draft 阶段的 file_path）
  → 读取原稿 Markdown 文件
  → 生成平台适配文件 + content_stages 表
```

每个 Skill 的输入查询模式：
```sql
-- 获取上一阶段的最新产出
SELECT file_path, metadata FROM content_stages
WHERE content_id = '[content_id]' AND stage = '[上一阶段]'
ORDER BY version DESC LIMIT 1
```

---

## 6. 创作工作流完整串联

```
┌─────────────────────────────────────────────────────────┐
│                     入口层                               │
├──────────────────────┬──────────────────────────────────┤
│   热榜 → 选题池      │        用户灵感                   │
│   (recommended_topics)│     (contents.status=idea)       │
│   [topic-generator]  │     [idea-researcher]             │
│         │            │              │                    │
│     采纳选题          │         researching               │
│         │            │     Agent 调研 → research 文件     │
│         └────────────┴──────────┐                       │
│                                 ▼                       │
│                           scripting                     │
│                       [script-create]                   │
│              Agent 生成脚本 → script 文件                 │
│                    content_stages(script)                │
│                                 │                       │
│                                 ▼                       │
│                            creating                     │
│                       [content-creator]                  │
│           创建内容文件夹 + 原稿 + 素材                     │
│                    content_stages(draft)                 │
│              contents.content_dir_path 更新              │
│                                 │                       │
│                                 ▼                       │
│                            adapting                     │
│                      [platform-adapter]                  │
│            按 target_platforms 生成适配版本               │
│              content_stages(platform_adapt) × N          │
│                                 │                       │
│                                 ▼                       │
│                           scheduled                     │
│                  加入 publish_queue 排期                  │
│                                 │                       │
│                                 ▼                       │
│                           published                     │
│              publish_records 记录发布结果                 │
│              review_tasks 定时采集数据                    │
│                                 │                       │
│                                 ▼                       │
│                            archived                     │
└─────────────────────────────────────────────────────────┘
```

### 每步操作的原子逻辑

每个 Skill 执行时，完成两个原子操作：

1. **写文件** — 将产出内容写入对应路径的 Markdown 文件
2. **更新数据库** — 创建 content_stages 记录 + 更新 contents.status

```sql
-- 1. 创建阶段记录
INSERT INTO content_stages (id, content_id, stage, file_path, status, version, source_type, metadata, created_at, updated_at)
VALUES ('[uuid]', '[content_id]', '[stage]', '[相对路径]', 'completed', [版本号], 'agent', '[JSON]', datetime('now'), datetime('now'));

-- 2. 更新内容状态
UPDATE contents SET status = '[下一状态]', updated_at = datetime('now') WHERE id = '[content_id]';
```

---

## 7. 需要改动的文件清单

### 7.1 数据库层（Sprouty AI）

| 文件 | 改动 |
|------|------|
| `packages/shared/src/db/schema.ts` | 精简 contents 表字段，新增 content_stages 表定义 |
| `packages/shared/src/db/types.ts` | 更新 Content 类型，新增 ContentStage 类型，更新 ContentStatus 枚举 |
| `packages/shared/src/db/migrations.ts` | 新增 v10 迁移：重建 contents 表 + CREATE content_stages（开发阶段无需数据迁移） |
| `packages/shared/src/db/repositories/contents.ts` | 精简字段相关的 CRUD 方法 |
| `packages/shared/src/db/repositories/content-versions.ts` | 删除（废弃） |
| 新增 `packages/shared/src/db/repositories/content-stages.ts` | content_stages 的 CRUD |

### 7.2 服务层（Sprouty AI）

| 文件 | 改动 |
|------|------|
| `packages/shared/src/services/topic-recommend-service.ts` | `adoptTopic()` 逻辑调整：同时创建 content_stages |

### 7.3 APP 配置（skill-marketplace/apps/app-creator-media）

#### statuses/config.json

将 `reviewing`（审核中）改为 `adapting`（适配中）：

```json
{
  "version": 2,
  "statuses": [
    { "id": "idea",        "label": "选题",   "category": "open",   "isDefault": true, "order": 0 },
    { "id": "researching", "label": "研究中", "category": "open",   "order": 1 },
    { "id": "scripting",   "label": "写脚本", "category": "open",   "order": 2 },
    { "id": "creating",    "label": "创作中", "category": "open",   "order": 3 },
    { "id": "adapting",    "label": "适配中", "category": "open",   "order": 4 },
    { "id": "scheduled",   "label": "待发布", "category": "open",   "order": 5 },
    { "id": "published",   "label": "已发布", "category": "closed", "order": 6 },
    { "id": "archived",    "label": "已归档", "category": "closed", "order": 7 }
  ],
  "defaultStatusId": "idea"
}
```

改动：`reviewing` → `adapting`，新增 `adapting.svg` 图标（可复用 `reviewing.svg`）

#### manifest.json

1. **`capabilities.skills`** — 更新为新的 skill 列表：

```json
{
  "capabilities": {
    "skills": [
      "topic-generator",
      "idea-researcher",
      "script-create",
      "content-creator",
      "platform-adapter",
      "hot-topic-scout",
      "competitor-monitor",
      "auto-publisher",
      "data-review",
      "video-script-create"
    ]
  }
}
```

2. **`workspace.directoryStructure`** — 清空，目录由 Skill 运行时动态创建：

```json
{
  "workspace": {
    "directoryStructure": {}
  }
}
```

### 7.4 Skill 文件（skill-marketplace/skills）

| Skill | 改动 |
|-------|------|
| `topic-generator/SKILL.md` | 适配新目录结构（项目名称前缀 + 文件名分隔符 `_`） |
| `script-create/SKILL.md` | 适配新数据库操作（content_stages 替代 contents.script_path） |
| 新建 `idea-researcher/` | 灵感调研 Skill |
| 新建 `content-creator/` | 内容创作 Skill |
| 新建 `platform-adapter/` | 平台适配 Skill（含 references/ 平台规范） |
