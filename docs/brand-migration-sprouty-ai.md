# 智小芽 (Sprouty AI) 品牌迁移设计文档

> 版本: 1.0
> 日期: 2026-02-07
> 状态: 待审核

---

## 1. 背景

项目原名 **Craft Agents**，后更名为 **CreatorFlow**，现品牌统一为 **智小芽 (Sprouty AI)**。

当前代码中存在三代品牌名混用：
- `craft` / `Craft` — 最早期命名（环境变量、OAuth 类型、域名）
- `creator-flow` / `CreatorFlow` — 中期命名（包名、目录、配置路径）
- `sprouty` / `智小芽` — 当前品牌（仅部分插件 manifest 和 Electron 应用名）

本次迁移目标：**将所有内部命名统一到 `sprouty-ai` 命名空间**，同时完成插件目录结构调整。

---

## 2. 命名规范

| 用途 | 旧值 | 新值 |
|------|------|------|
| 配置目录 | `~/.creator-flow/` | `~/.sprouty-ai/` |
| 工作区数据目录 | `{workspace}/.creator-flow/` | `{workspace}/.sprouty-ai/` |
| 环境变量前缀 | `CREATOR_FLOW_` / `CRAFT_` | `SPROUTY_` |
| NPM 包命名空间 | `@creator-flow/*` | `@sprouty-ai/*` |
| Deeplink scheme | `creatorflow://` | `sproutyai://` |
| 插件 manifest name（全局） | `sprouty` | `sprouty` (不变) |
| 插件 manifest name（工作区） | `sprouty-workspace` | `sprouty-workspace` (不变) |
| Electron appId | `com.zhixiaoya.app` | `com.zhixiaoya.app` (不变) |
| Electron productName | `智小芽` | `智小芽` (不变) |
| 代码中的类/函数名 | `CreatorFlowAgent` 等 | `SproutyAgent` 等 |

---

## 3. 目录结构调整

### 3.1 全局配置目录

```
旧: ~/.creator-flow/
新: ~/.sprouty-ai/
    ├── config.json                    # 应用配置
    ├── credentials.enc                # 加密凭证
    ├── drafts.json                    # 会话草稿
    ├── provider-domains.json          # 提供商域名缓存
    ├── theme.json                     # 主题配置
    ├── themes/                        # 预设主题
    ├── tool-icons/                    # 工具图标
    ├── docs/                          # 内置文档
    ├── permissions/                   # 权限配置
    │   └── default.json
    ├── marketplace/                   # 市场缓存
    │   └── cache/
    ├── workspaces/                    # 默认工作区存储
    │   └── {workspace-id}/
    └── .claude-plugin/                # 全局插件 manifest（SDK 兼容）
        └── plugin.json
```

### 3.2 工作区数据目录

```
旧: {workspace}/.creator-flow/
新: {workspace}/.sprouty-ai/
    ├── config.json                    # 工作区配置
    ├── sessions/                      # 会话数据
    ├── sources/                       # 数据源
    ├── skills/                        # 工作区技能（用户自由管理）
    │   └── {skill-name}/
    │       └── SKILL.md
    ├── prompts/                       # 提示词
    ├── guides/                        # 指南
    ├── labels/                        # 标签
    ├── statuses/                      # 状态
    ├── plugins/                       # 插件目录 ← 新位置
    │   └── {plugin-name}/
    │       ├── .claude-plugin/        # Claude SDK 兼容结构
    │       │   ├── plugin.json
    │       │   └── translations.json
    │       ├── commands/              # 插件斜杠命令
    │       │   └── {command}.md
    │       ├── skills/                # 插件自带技能（随插件安装/删除）
    │       │   └── {skill-name}/
    │       │       └── SKILL.md
    │       └── agents/                # 插件子智能体定义（可选）
    │           └── {agent}.json
    ├── agents/                        # 工作区子智能体定义
    │   └── {agent}.json
    └── .claude-plugin/                # 工作区插件 manifest（SDK 兼容）
        └── plugin.json
```

### 3.3 插件目录迁移

**旧结构**（嵌套在 `.claude-plugin/` 下）：
```
{workspace}/.creator-flow/.claude-plugin/productivity/
    ├── plugin.json
    ├── commands/
    └── skills/
```

**新结构**（独立 `plugins/` 目录）：
```
{workspace}/.sprouty-ai/plugins/productivity/
    ├── .claude-plugin/
    │   ├── plugin.json
    │   └── translations.json
    ├── commands/
    ├── skills/
    └── agents/
```

---

## 4. 需要修改的文件清单

### 4.1 核心路径定义

| 文件 | 修改内容 |
|------|---------|
| `packages/shared/src/config/paths.ts` | `CREATOR_FLOW_CONFIG_DIR` → `SPROUTY_CONFIG_DIR`，`.creator-flow` → `.sprouty-ai` |
| `packages/shared/src/workspaces/storage.ts` | `CONFIG_DIR` 引用更新，`WORKSPACE_DATA_DIR` 从 `.creator-flow` → `.sprouty-ai` |

### 4.2 环境变量

| 旧变量 | 新变量 | 所在文件 |
|--------|--------|---------|
| `CREATOR_FLOW_CONFIG_DIR` | `SPROUTY_CONFIG_DIR` | `config/paths.ts` |
| `CRAFT_DEBUG` | `SPROUTY_DEBUG` | `electron/main/index.ts`, 多处 |
| `CRAFT_LOCAL_MCP_ENABLED` | `SPROUTY_LOCAL_MCP_ENABLED` | `workspaces/storage.ts` |
| `CRAFT_APP_NAME` | `SPROUTY_APP_NAME` | `electron/main/index.ts` |
| `CREATORFLOW_DEEPLINK_SCHEME` | `SPROUTY_DEEPLINK_SCHEME` | `electron/main/index.ts`, `sessions.ts` |

### 4.3 NPM 包名

| 文件 | 旧值 | 新值 |
|------|------|------|
| `package.json` (root) | `creator-flow` | `sprouty-ai` |
| `packages/shared/package.json` | `@creator-flow/shared` | `@sprouty-ai/shared` |
| `apps/electron/package.json` | `@creator-flow/electron` | `@sprouty-ai/electron` |
| `apps/server/package.json` | `@creator-flow/server` | `@sprouty-ai/server` |
| `packages/ui/package.json` | `@creator-flow/ui` | `@sprouty-ai/ui` |
| `packages/core/package.json` | `@creator-flow/core` | `@sprouty-ai/core` |
| `packages/mermaid/package.json` | `@creator-flow/mermaid` | `@sprouty-ai/mermaid` |

### 4.4 代码中的 import 路径

所有 `@creator-flow/*` 的 import 需要替换为 `@sprouty-ai/*`，涉及文件约 100+。

### 4.5 TypeScript 类名/常量名

| 旧名 | 新名 | 所在文件 |
|------|------|---------|
| `CreatorFlowAgent` | `SproutyAgent` | `creator-flow-agent.ts` → `sprouty-agent.ts` |
| `CREATOR_FLOW_LOGO` | `SPROUTY_LOGO` | `branding.ts` |
| `CREATOR_FLOW_LOGO_HTML` | `SPROUTY_LOGO_HTML` | `branding.ts` |
| `CRAFT_LOGO` / `CRAFT_LOGO_HTML` | 删除 | `branding.ts` |

### 4.6 Deeplink scheme

| 文件 | 修改内容 |
|------|---------|
| `apps/electron/src/main/index.ts` | `creatorflow` → `sproutyai` |
| `apps/electron/src/main/sessions.ts` | deeplink scheme 引用 |
| `scripts/electron-dev.ts` | deeplink scheme 引用 |

### 4.7 插件系统路径

| 文件 | 修改内容 |
|------|---------|
| `packages/shared/src/workspaces/storage.ts` | 插件目录从 `.claude-plugin/` 迁移到 `plugins/` |
| `packages/shared/src/agent/creator-flow-agent.ts` | SDK plugins 配置路径更新 |
| `packages/shared/src/agent/slash-command-translations.ts` | 命令扫描路径更新 |
| `apps/electron/src/main/sessions.ts` | 命令推送路径更新 |

### 4.8 OAuth / 凭证

| 文件 | 修改内容 |
|------|---------|
| `packages/shared/src/credentials/types.ts` | `craft_oauth` → `claude_oauth`（已完成） |
| `packages/shared/src/credentials/backends/secure-storage.ts` | 路径引用更新 |

### 4.9 文档和注释

| 文件 | 修改内容 |
|------|---------|
| `CLAUDE.md` (root + packages) | 品牌名和路径引用 |
| `README.md` | 品牌名、仓库 URL |
| `CONTRIBUTING.md` | 仓库名引用 |
| `TRADEMARK.md` | 品牌名引用 |
| `CODE_OF_CONDUCT.md` | 联系邮箱 |
| `SECURITY.md` | 联系邮箱 |
| 代码注释中的 `CreatorFlow` / `creator-flow` | 统一替换 |

### 4.10 构建脚本

| 文件 | 修改内容 |
|------|---------|
| `scripts/install-app.sh` | 品牌名引用 |
| `scripts/install-app.ps1` | 品牌名引用 |
| `scripts/electron-dev.ts` | 环境变量名 |

---

## 5. 迁移策略

### 5.1 向后兼容

为避免用户数据丢失，需要在首次启动时自动迁移：

```typescript
// config/paths.ts 中添加迁移逻辑
function migrateConfigDir(): string {
  const newDir = join(homedir(), '.sprouty-ai');
  const oldDir = join(homedir(), '.creator-flow');

  if (!existsSync(newDir) && existsSync(oldDir)) {
    // 旧目录存在但新目录不存在 → 重命名
    renameSync(oldDir, newDir);
    console.log(`Migrated config: ${oldDir} → ${newDir}`);
  }

  return newDir;
}
```

工作区数据目录同理：
```typescript
// workspaces/storage.ts 中添加迁移逻辑
function migrateWorkspaceDataDir(workspaceRoot: string): void {
  const newDir = join(workspaceRoot, '.sprouty-ai');
  const oldDir = join(workspaceRoot, '.creator-flow');

  if (!existsSync(newDir) && existsSync(oldDir)) {
    renameSync(oldDir, newDir);
  }
}
```

### 5.2 环境变量兼容

新旧环境变量同时支持，优先使用新变量：

```typescript
export const CONFIG_DIR =
  process.env.SPROUTY_CONFIG_DIR ||
  process.env.CREATOR_FLOW_CONFIG_DIR ||  // 向后兼容
  migrateConfigDir();
```

### 5.3 插件目录迁移

插件从 `.claude-plugin/{name}/` 迁移到 `plugins/{name}/`：

```typescript
function migratePluginDir(workspaceDataDir: string): void {
  const newPluginsDir = join(workspaceDataDir, 'plugins');
  const oldClaudePluginDir = join(workspaceDataDir, '.claude-plugin');

  if (!existsSync(oldClaudePluginDir)) return;

  // 扫描旧目录中的嵌套插件
  for (const entry of readdirSync(oldClaudePluginDir, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue;
    if (entry.name === 'commands' || entry.name === 'skills') continue; // 跳过顶层目录

    const oldPath = join(oldClaudePluginDir, entry.name);
    const newPath = join(newPluginsDir, entry.name);

    if (existsSync(join(oldPath, 'plugin.json')) || existsSync(join(oldPath, '.claude-plugin', 'plugin.json'))) {
      mkdirSync(newPluginsDir, { recursive: true });
      renameSync(oldPath, newPath);
    }
  }
}
```

---

## 6. 执行计划

### 阶段 1：核心路径和常量（影响面最大，优先处理）

1. 修改 `config/paths.ts` — 配置目录路径 + 迁移逻辑
2. 修改 `workspaces/storage.ts` — 工作区数据目录 + 迁移逻辑
3. 修改环境变量名（保留旧变量兼容）
4. 修改 Deeplink scheme

### 阶段 2：NPM 包名和 import 路径

1. 修改所有 `package.json` 的包名
2. 批量替换 `@creator-flow/` → `@sprouty-ai/` import 路径
3. 更新 `tsconfig.json` 中的路径映射（如有）
4. 更新 monorepo workspace 配置

### 阶段 3：类名和文件名

1. `creator-flow-agent.ts` → `sprouty-agent.ts`（文件重命名 + 类名替换）
2. `branding.ts` 常量名更新
3. 其他 TypeScript 类名/常量名替换

### 阶段 4：插件目录结构

1. 插件路径从 `.claude-plugin/{name}/` → `plugins/{name}/`
2. 更新 SDK plugins 配置扫描逻辑
3. 更新命令扫描路径
4. 添加插件目录迁移逻辑

### 阶段 5：文档和清理

1. 更新所有文档中的品牌名和路径
2. 删除旧的向后兼容导出（`CRAFT_LOGO` 等）
3. 更新 `.gitignore` 中的路径
4. 更新构建脚本

---

## 7. 风险和注意事项

### 7.1 用户数据迁移
- 必须在首次启动时自动迁移 `~/.creator-flow/` → `~/.sprouty-ai/`
- 必须迁移每个工作区的 `.creator-flow/` → `.sprouty-ai/`
- 迁移失败时不能丢失数据，应保留旧目录

### 7.2 多实例开发
- 当前支持 `~/.creator-flow-1/`、`~/.creator-flow-2/` 等多实例
- 迁移后需支持 `~/.sprouty-ai-1/`、`~/.sprouty-ai-2/`
- 环境变量 `SPROUTY_CONFIG_DIR` 覆盖一切

### 7.3 已发布版本兼容
- 已安装的用户升级后，旧配置目录必须自动迁移
- 凭证文件（`credentials.enc`）迁移后必须仍可解密
- 工作区路径存储在 `config.json` 中，迁移时需更新内部引用

### 7.4 构建系统
- NPM 包名变更后，`bun.lock` 需要重新生成
- monorepo workspace 引用需要同步更新
- CI/CD 脚本中的包名引用需要更新

### 7.5 不变更的内容
- `com.zhixiaoya.app` — Electron appId 不变（避免 macOS 签名问题）
- `智小芽` — 用户可见的应用名不变
- `craft.do` 域名引用 — 这是上游服务域名，不属于本项目品牌
- OAuth credential type 名称 — `claude_oauth` 等已是通用名，不含旧品牌

---

## 8. 验证清单

- [ ] `~/.sprouty-ai/` 目录正确创建
- [ ] 旧 `~/.creator-flow/` 自动迁移到新路径
- [ ] 工作区 `.sprouty-ai/` 目录正确创建
- [ ] 旧工作区 `.creator-flow/` 自动迁移
- [ ] 插件目录 `plugins/{name}/` 结构正确
- [ ] SDK 正确加载插件（commands + skills）
- [ ] 斜杠命令菜单正常显示
- [ ] 技能列表正常加载
- [ ] Deeplink `sproutyai://` 正常工作
- [ ] 凭证加解密正常
- [ ] 多实例开发正常
- [ ] 构建打包正常
- [ ] 所有 import 路径正确解析
- [ ] TypeScript 编译无错误
