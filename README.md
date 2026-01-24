# CreatorFlow

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

CreatorFlow 是我们（craft.do 团队）为了更高效地与 AI 智能体协作而构建的工具。它支持直观的多任务处理、无缝连接任何 API 或服务、共享会话，以及更注重文档（而非代码）的工作流程——所有这些都在一个精美流畅的界面中呈现。

它基于 Claude Agent SDK 构建，延续了 Claude Code 的优秀特性，同时改进了我们认为需要提升的领域。

它以 Agent Native 软件原则为核心设计，开箱即用且高度可定制。这是同类产品中的先驱之一。

CreatorFlow 采用 Apache 2.0 许可证开源——您可以自由修改、重组任何内容。我们自己就完全使用 CreatorFlow 来构建 CreatorFlow，不使用任何代码编辑器——因此，任何定制都只需一个提示词即可实现。

我们构建 CreatorFlow 是因为我们希望有一种更好、更有主见的方式（最好是非命令行的方式）来与世界上最强大的智能体协作。我们将基于自身的经验和直觉持续改进它。

<img width="1578" height="894" alt="image" src="https://github.com/user-attachments/assets/3f1f2fe8-7cf6-4487-99ff-76f6c8c0a3fb" />

## 安装

### 一键安装（推荐）

**macOS / Linux:**
```bash
curl -fsSL https://agents.craft.do/install-app.sh | bash
```

**Windows (PowerShell):**
```powershell
irm https://agents.craft.do/install-app.ps1 | iex
```

### 从源码构建

```bash
git clone https://github.com/lukilabs/craft-agents-oss.git
cd craft-agents-oss
bun install
bun run electron:start
```

## 功能特性

- **多会话收件箱**：桌面应用，支持会话管理、状态工作流和标记功能
- **Claude Code 体验**：流式响应、工具可视化、实时更新
- **Craft MCP 集成**：访问 32+ Craft 文档工具（区块、集合、搜索、任务）
- **数据源**：连接 MCP 服务器、REST API（Google、Slack、Microsoft）和本地文件系统
- **权限模式**：三级权限系统（探索、请求编辑、自动），支持自定义规则
- **后台任务**：运行长时间操作并跟踪进度
- **动态状态系统**：可自定义的会话工作流状态（待办、进行中、已完成等）
- **主题系统**：支持应用级和工作区级的级联主题
- **多文件差异对比**：VS Code 风格的窗口，查看一次对话中的所有文件更改
- **技能**：按工作区存储的专用智能体指令
- **文件附件**：拖放图片、PDF、Office 文档并自动转换
- **多语言支持**：内置中英双语，可在设置中切换语言

## 快速开始

1. **启动应用**：安装完成后启动
2. **选择 API 连接**：使用您自己的 Anthropic API 密钥或 Claude Max 订阅
3. **创建工作区**：设置工作区来组织您的会话
4. **连接数据源**（可选）：添加 MCP 服务器、REST API 或本地文件系统
5. **开始对话**：创建会话并与 Claude 交互

## 桌面应用功能

### 会话管理

- **收件箱/归档**：按工作流状态组织会话
- **标记**：标记重要会话以便快速访问
- **状态工作流**：待办 → 进行中 → 待审核 → 已完成
- **会话命名**：AI 生成标题或手动命名
- **会话持久化**：完整对话历史保存到磁盘

### 数据源

连接外部数据源到您的工作区：

| 类型 | 示例 |
|------|------|
| **MCP 服务器** | Craft、Linear、GitHub、Notion、自定义服务器 |
| **REST API** | Google（Gmail、日历、云盘）、Slack、Microsoft |
| **本地文件** | 文件系统、Obsidian 仓库、Git 仓库 |

### 权限模式

| 模式 | 显示名称 | 行为 |
|------|----------|------|
| `safe` | 探索 | 只读，阻止所有写入操作 |
| `ask` | 请求编辑 | 需要批准（默认） |
| `allow-all` | 自动 | 自动批准所有命令 |

使用 **SHIFT+TAB** 在聊天界面中切换模式。

### 键盘快捷键

| 快捷键 | 操作 |
|--------|------|
| `Cmd+N` | 新建聊天 |
| `Cmd+1/2/3` | 聚焦侧边栏/列表/聊天 |
| `Cmd+/` | 快捷键对话框 |
| `SHIFT+TAB` | 切换权限模式 |
| `Enter` | 发送消息 |
| `Shift+Enter` | 换行 |

## 架构

```
craft-agent/
├── apps/
│   └── electron/              # 桌面 GUI（主要）
│       └── src/
│           ├── main/          # Electron 主进程
│           ├── preload/       # 上下文桥接
│           └── renderer/      # React UI (Vite + shadcn)
└── packages/
    ├── core/                  # 共享类型
    └── shared/                # 业务逻辑
        └── src/
            ├── agent/         # CraftAgent，权限
            ├── auth/          # OAuth，令牌
            ├── config/        # 存储，偏好设置，主题
            ├── credentials/   # AES-256-GCM 加密存储
            ├── locale/        # 多语言服务
            ├── sessions/      # 会话持久化
            ├── sources/       # MCP、API、本地数据源
            └── statuses/      # 动态状态系统
```

## 开发

```bash
# 热重载开发
bun run electron:dev

# 构建并运行
bun run electron:start

# 类型检查
bun run typecheck:all

# 多语言翻译
bun run i18n:scan        # 扫描代码中的中文文本
bun run i18n:translate   # 翻译为目标语言

# 调试日志（写入 ~/Library/Logs/CreatorFlow/）
# 开发模式下自动启用日志
```

### 环境变量

OAuth 集成（Google、Slack、Microsoft）需要凭据。创建 `.env` 文件：

```bash
MICROSOFT_OAUTH_CLIENT_ID=your-client-id
GOOGLE_OAUTH_CLIENT_SECRET=your-google-client-secret
GOOGLE_OAUTH_CLIENT_ID=your-client-id.apps.googleusercontent.com
SLACK_OAUTH_CLIENT_ID=your-slack-client-id
SLACK_OAUTH_CLIENT_SECRET=your-slack-client-secret
```

参见 [Google Cloud Console](https://console.cloud.google.com/apis/credentials) 创建 OAuth 凭据。

## 配置

配置存储在 `~/.craft-agent/`：

```
~/.craft-agent/
├── config.json              # 主配置（工作区、认证类型）
├── credentials.enc          # 加密凭据（AES-256-GCM）
├── preferences.json         # 用户偏好设置
├── theme.json               # 应用级主题
└── workspaces/
    └── {id}/
        ├── config.json      # 工作区设置
        ├── theme.json       # 工作区主题覆盖
        ├── sessions/        # 会话数据（JSONL）
        ├── sources/         # 已连接的数据源
        ├── skills/          # 自定义技能
        └── statuses/        # 状态配置
```

## 高级功能

### 大型响应处理

超过约 60KB 的工具响应会自动使用 Claude Haiku 进行意图感知的摘要。`_intent` 字段会注入到 MCP 工具模式中以保持摘要焦点。

### 深度链接

外部应用可以使用 `craftagents://` URL 导航：

```
craftagents://allChats                    # 所有聊天视图
craftagents://allChats/chat/session123    # 特定聊天
craftagents://settings                    # 设置
craftagents://sources/source/github       # 数据源信息
craftagents://action/new-chat             # 创建新聊天
```

### 多语言支持

CreatorFlow 内置多语言支持系统，使用方式：

```typescript
// 在 React 组件中使用
import { useT } from '@/context/LocaleContext'

function MyComponent() {
  const t = useT()
  return <div>{t('设置')}</div>
}

// 在非 React 代码中使用
import { t } from '@craft-agent/shared/locale'
console.log(t('操作成功'))
```

语言包位于 `apps/electron/src/renderer/locales/`，支持自动翻译脚本。

## 技术栈

| 层 | 技术 |
|----|------|
| 运行时 | [Bun](https://bun.sh/) |
| AI | [@anthropic-ai/claude-agent-sdk](https://www.npmjs.com/package/@anthropic-ai/claude-agent-sdk) |
| 桌面 | [Electron](https://www.electronjs.org/) + React |
| UI | [shadcn/ui](https://ui.shadcn.com/) + Tailwind CSS v4 |
| 构建 | esbuild (main) + Vite (renderer) |
| 凭据 | AES-256-GCM 加密文件存储 |

## 许可证

本项目采用 Apache License 2.0 许可证 - 详见 [LICENSE](LICENSE) 文件。

### 第三方许可证

本项目使用 [Claude Agent SDK](https://www.npmjs.com/package/@anthropic-ai/claude-agent-sdk)，受 [Anthropic 商业服务条款](https://www.anthropic.com/legal/commercial-terms) 约束。

### 商标

"Craft" 和 "CreatorFlow" 是 Craft Docs Ltd. 的商标。使用指南请参见 [TRADEMARK.md](TRADEMARK.md)。

## 贡献

我们欢迎贡献！请参阅 [CONTRIBUTING.md](CONTRIBUTING.md) 了解指南。

## 安全

### 本地 MCP 服务器隔离

在启动本地 MCP 服务器（stdio 传输）时，敏感环境变量会被过滤以防止凭据泄露到子进程。被阻止的变量包括：

- `ANTHROPIC_API_KEY`、`CLAUDE_CODE_OAUTH_TOKEN`（应用认证）
- `AWS_ACCESS_KEY_ID`、`AWS_SECRET_ACCESS_KEY`、`AWS_SESSION_TOKEN`
- `GITHUB_TOKEN`、`GH_TOKEN`、`OPENAI_API_KEY`、`GOOGLE_API_KEY`、`STRIPE_SECRET_KEY`、`NPM_TOKEN`

要显式将环境变量传递给特定 MCP 服务器，请在源配置中使用 `env` 字段。

报告安全漏洞，请参见 [SECURITY.md](SECURITY.md)。
