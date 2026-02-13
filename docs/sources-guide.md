# Sources 数据源配置指南

Sources（数据源）是 Sprouty AI 中连接外部数据的核心机制。通过配置数据源，AI 助手可以访问各种外部服务和本地文件，从而提供更强大的辅助能力。

## 概述

### 什么是 Sources？

Sources 是对外部数据连接的抽象，支持三种类型：

| 类型 | 说明 | 典型场景 |
|------|------|----------|
| `api` | REST API 接口 | GitHub、Unsplash、Notion 等云服务 |
| `mcp` | MCP 协议服务 | Linear、Slack 等 MCP 服务器 |
| `local` | 本地文件系统 | 项目目录、笔记库、素材文件夹 |

### 存储位置

Sources 配置存储在工作区目录下：

```
~/.creator-flow/workspaces/{workspaceId}/sources/{sourceSlug}/
├── config.json   # 数据源配置
├── guide.md      # 使用指南（可选）
└── icon.png      # 图标文件（可选）
```

### 工作流程

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  应用预设配置    │ ──▶ │  工作区初始化     │ ──▶ │  用户启用配置    │
│  (bundled-apps) │     │  (复制到工作区)   │     │  (添加认证信息)  │
└─────────────────┘     └──────────────────┘     └─────────────────┘
```

1. **应用预设**：内置应用可包含预设的数据源配置
2. **工作区初始化**：创建工作区时，预设配置自动复制到工作区
3. **用户配置**：用户启用数据源并配置认证信息

## 配置结构

### config.json 完整结构

```typescript
interface FolderSourceConfig {
  // 基础信息
  id: string;              // 唯一标识，格式: {slug}_{random8}
  name: string;            // 显示名称
  slug: string;            // URL 友好的标识符
  enabled: boolean;        // 是否启用
  
  // 提供者和类型
  provider: string;        // 提供者标识 (github, google, local 等)
  type: SourceType;        // 'api' | 'mcp' | 'local'
  
  // 类型特定配置 (三选一)
  api?: ApiSourceConfig;
  mcp?: McpSourceConfig;
  local?: LocalSourceConfig;
  
  // 显示信息
  icon?: string;           // emoji 或图标 URL
  tagline?: string;        // 简短描述
  
  // 状态信息
  isAuthenticated?: boolean;
  connectionStatus?: 'connected' | 'needs_auth' | 'failed' | 'untested';
  connectionError?: string;
  lastTestedAt?: number;
  
  // 时间戳
  createdAt?: number;
  updatedAt?: number;
}
```

## API 类型配置

API 类型用于连接 REST API 服务。

### 配置结构

```typescript
interface ApiSourceConfig {
  baseUrl: string;                    // API 基础 URL
  authType: ApiAuthType;              // 认证类型
  headerName?: string;                // Header 认证时的头名称
  queryParam?: string;                // Query 认证时的参数名
  authScheme?: string;                // Bearer 认证的 scheme (默认 "Bearer")
  defaultHeaders?: Record<string, string>;  // 默认请求头
  testEndpoint?: ApiTestEndpoint;     // 连接测试端点
  
  // OAuth 相关 (特定提供者)
  googleService?: GoogleService;
  googleScopes?: string[];
  slackService?: SlackService;
  slackUserScopes?: string[];
  microsoftService?: MicrosoftService;
  microsoftScopes?: string[];
}

type ApiAuthType = 'bearer' | 'header' | 'query' | 'basic' | 'none';
```

### 认证类型说明

| 认证类型 | 说明 | 配置字段 |
|---------|------|---------|
| `bearer` | Bearer Token | `authScheme` (默认 "Bearer") |
| `header` | 自定义 Header | `headerName` (如 "X-API-Key") |
| `query` | URL 查询参数 | `queryParam` (如 "api_key") |
| `basic` | HTTP Basic Auth | 用户名:密码 |
| `none` | 无认证 | - |

### 示例：GitHub API

```json
{
  "id": "github_a1b2c3d4",
  "name": "GitHub",
  "slug": "github",
  "enabled": true,
  "provider": "github",
  "type": "api",
  "api": {
    "baseUrl": "https://api.github.com",
    "authType": "bearer",
    "authScheme": "Bearer",
    "defaultHeaders": {
      "Accept": "application/vnd.github+json",
      "X-GitHub-Api-Version": "2022-11-28"
    },
    "testEndpoint": {
      "method": "GET",
      "path": "/user"
    }
  },
  "icon": "🐙",
  "tagline": "访问 GitHub 仓库、Issues、Pull Requests"
}
```

### 示例：Unsplash API

```json
{
  "id": "unsplash_e5f6g7h8",
  "name": "Unsplash",
  "slug": "unsplash",
  "enabled": true,
  "provider": "unsplash",
  "type": "api",
  "api": {
    "baseUrl": "https://api.unsplash.com",
    "authType": "header",
    "headerName": "Authorization",
    "authScheme": "Client-ID",
    "defaultHeaders": {
      "Accept-Version": "v1"
    },
    "testEndpoint": {
      "method": "GET",
      "path": "/me"
    }
  },
  "icon": "📷",
  "tagline": "高质量免费图片素材库"
}
```

### 示例：OpenAI API

```json
{
  "id": "openai_i9j0k1l2",
  "name": "OpenAI",
  "slug": "openai",
  "enabled": true,
  "provider": "openai",
  "type": "api",
  "api": {
    "baseUrl": "https://api.openai.com/v1",
    "authType": "bearer",
    "authScheme": "Bearer",
    "defaultHeaders": {
      "Content-Type": "application/json"
    },
    "testEndpoint": {
      "method": "GET",
      "path": "/models"
    }
  },
  "icon": "🤖",
  "tagline": "OpenAI GPT 模型 API"
}
```

### 示例：自定义 API (Header 认证)

```json
{
  "id": "custom-api_m3n4o5p6",
  "name": "内部 API",
  "slug": "custom-api",
  "enabled": true,
  "provider": "custom",
  "type": "api",
  "api": {
    "baseUrl": "https://api.internal.company.com",
    "authType": "header",
    "headerName": "X-API-Key",
    "defaultHeaders": {
      "Content-Type": "application/json"
    },
    "testEndpoint": {
      "method": "GET",
      "path": "/health"
    }
  },
  "icon": "🔧",
  "tagline": "公司内部服务 API"
}
```

### 示例：Google API (OAuth)

```json
{
  "id": "google-calendar_q7r8s9t0",
  "name": "Google 日历",
  "slug": "google-calendar",
  "enabled": true,
  "provider": "google",
  "type": "api",
  "api": {
    "baseUrl": "https://www.googleapis.com/calendar/v3",
    "authType": "bearer",
    "googleService": "calendar",
    "testEndpoint": {
      "method": "GET",
      "path": "/users/me/calendarList"
    }
  },
  "icon": "📅",
  "tagline": "Google 日历事件管理"
}
```

## MCP 类型配置

MCP (Model Context Protocol) 类型用于连接 MCP 服务器。

### 配置结构

```typescript
interface McpSourceConfig {
  transport?: McpTransport;  // 'http' | 'sse' | 'stdio'
  
  // HTTP/SSE 传输
  url?: string;              // 服务器 URL
  authType?: SourceMcpAuthType;  // 'oauth' | 'bearer' | 'none'
  clientId?: string;         // OAuth Client ID
  
  // Stdio 传输 (本地 MCP 服务器)
  command?: string;          // 启动命令
  args?: string[];           // 命令参数
  env?: Record<string, string>;  // 环境变量
}
```

### 示例：Linear MCP (OAuth)

```json
{
  "id": "linear_u1v2w3x4",
  "name": "Linear",
  "slug": "linear",
  "enabled": true,
  "provider": "linear",
  "type": "mcp",
  "mcp": {
    "transport": "http",
    "url": "https://mcp.linear.app",
    "authType": "oauth",
    "clientId": "your-client-id"
  },
  "icon": "📊",
  "tagline": "项目管理和 Issue 追踪"
}
```

### 示例：本地 MCP 服务器 (Stdio)

```json
{
  "id": "sqlite_y5z6a7b8",
  "name": "SQLite 数据库",
  "slug": "sqlite",
  "enabled": true,
  "provider": "sqlite",
  "type": "mcp",
  "mcp": {
    "transport": "stdio",
    "command": "npx",
    "args": ["-y", "@anthropic-ai/mcp-server-sqlite", "/path/to/database.db"],
    "env": {}
  },
  "icon": "🗄️",
  "tagline": "本地 SQLite 数据库访问"
}
```

### 示例：文件系统 MCP 服务器

```json
{
  "id": "filesystem_c9d0e1f2",
  "name": "文件系统",
  "slug": "filesystem",
  "enabled": true,
  "provider": "filesystem",
  "type": "mcp",
  "mcp": {
    "transport": "stdio",
    "command": "npx",
    "args": ["-y", "@anthropic-ai/mcp-server-filesystem", "/Users/me/Documents"],
    "env": {}
  },
  "icon": "📁",
  "tagline": "本地文件系统访问"
}
```

### 示例：Slack MCP

```json
{
  "id": "slack_g3h4i5j6",
  "name": "Slack",
  "slug": "slack",
  "enabled": true,
  "provider": "slack",
  "type": "mcp",
  "mcp": {
    "transport": "http",
    "url": "https://mcp.slack.com",
    "authType": "oauth",
    "clientId": "your-slack-client-id"
  },
  "icon": "💬",
  "tagline": "Slack 消息和频道管理"
}
```

## Local 类型配置

Local 类型用于连接本地文件系统目录。

### 配置结构

```typescript
interface LocalSourceConfig {
  path: string;      // 目录路径
  format?: string;   // 格式提示: 'filesystem' | 'obsidian' | 'git' | 'sqlite'
}
```

### 路径格式支持

| 路径格式 | 说明 | 示例 |
|---------|------|------|
| `./xxx` | 相对于工作区根目录 | `./素材库`, `./项目` |
| `~/xxx` | 相对于用户主目录 | `~/Documents`, `~/Projects` |
| `/xxx` | 绝对路径 | `/Users/me/data` |

> **推荐**: 应用预设配置应使用 `./` 前缀的相对路径，这样路径会相对于工作区根目录展开。同时确保在 `manifest.json` 的 `directoryStructure` 中定义对应目录，这样创建工作区时会自动创建该目录。

### 示例：本地项目目录

```json
{
  "id": "local-projects_k7l8m9n0",
  "name": "本地项目",
  "slug": "local-projects",
  "enabled": true,
  "provider": "local",
  "type": "local",
  "local": {
    "path": "~/Projects",
    "format": "git"
  },
  "icon": "📂",
  "tagline": "管理本地代码项目"
}
```

### 示例：Obsidian 笔记库

```json
{
  "id": "obsidian-vault_o1p2q3r4",
  "name": "Obsidian 笔记",
  "slug": "obsidian-vault",
  "enabled": true,
  "provider": "local",
  "type": "local",
  "local": {
    "path": "~/Documents/ObsidianVault",
    "format": "obsidian"
  },
  "icon": "📝",
  "tagline": "Obsidian 知识库"
}
```

### 示例：素材库

```json
{
  "id": "media-library_s5t6u7v8",
  "name": "本地素材库",
  "slug": "media-library",
  "enabled": true,
  "provider": "local",
  "type": "local",
  "local": {
    "path": "~/CreatorMedia/素材库",
    "format": "filesystem"
  },
  "icon": "📁",
  "tagline": "管理图片、视频、音频素材"
}
```

### 示例：写作目录

```json
{
  "id": "writing_w9x0y1z2",
  "name": "写作目录",
  "slug": "writing",
  "enabled": true,
  "provider": "local",
  "type": "local",
  "local": {
    "path": "~/Documents/Writing",
    "format": "filesystem"
  },
  "icon": "✍️",
  "tagline": "文章草稿和写作素材"
}
```

## guide.md 使用指南

每个数据源可以包含一个 `guide.md` 文件，提供使用指南和上下文信息。AI 助手会读取这些信息来更好地理解如何使用该数据源。

### 推荐结构

```markdown
# 数据源名称

简短描述，说明这个数据源的用途。

## Scope

- 可以做什么
- 数据范围
- 主要功能

## Guidelines

### 使用规范
1. 规则一
2. 规则二

### 认证说明
如何获取和配置认证信息。

### 常用操作
- 操作一
- 操作二

## Context

补充上下文信息，如限制、注意事项等。

## API Notes

API 相关的技术说明（仅 API 类型需要）。
```

### 示例：GitHub guide.md

```markdown
# GitHub

访问 GitHub 代码仓库、Issues、Pull Requests 和项目管理功能。

## Scope

- 浏览和搜索代码仓库
- 管理 Issues 和 Pull Requests
- 查看提交历史和代码变更
- 项目看板和里程碑管理

## Guidelines

### API 认证
1. 访问 GitHub Settings > Developer settings > Personal access tokens
2. 生成 Personal Access Token (classic 或 fine-grained)
3. 选择需要的权限范围 (repo, read:user 等)

### 常用端点
- `GET /user` - 获取当前用户信息
- `GET /user/repos` - 列出用户仓库
- `GET /repos/{owner}/{repo}/issues` - 列出 Issues

## Context

速率限制：
- 认证用户: 5000 次/小时
- 未认证: 60 次/小时

## API Notes

响应格式为 JSON，分页使用 `page` 和 `per_page` 参数。
```

## 应用预设配置

应用可以在 `sources/` 目录下预设数据源配置，这些配置会在创建工作区时自动复制。

### 目录结构

```
bundled-apps/{app-id}/
├── manifest.json
├── AGENTS.md
├── labels/
├── statuses/
└── sources/
    ├── {source-slug-1}/
    │   ├── config.json
    │   └── guide.md
    └── {source-slug-2}/
        ├── config.json
        └── guide.md
```

### 预设配置注意事项

1. **默认禁用**：预设配置应设置 `enabled: false`，让用户主动启用
2. **状态标记**：设置 `connectionStatus: "needs_auth"` 或 `"untested"`
3. **路径占位符**：使用 `~` 作为用户主目录占位符
4. **中文友好**：名称和描述使用中文

### 示例：自媒体创作应用预设

```
bundled-apps/app-creator-media/sources/
├── unsplash/
│   ├── config.json    # Unsplash API 配置
│   └── guide.md       # 使用指南
└── media-library/
    ├── config.json    # 本地素材库配置
    └── guide.md       # 使用指南
```

## 凭证管理

数据源的敏感凭证（API Key、Token 等）不存储在 config.json 中，而是通过 Sprouty AI 的加密凭证系统管理。

### 凭证存储位置

```
~/.creator-flow/credentials.enc  # AES-256-GCM 加密存储
```

### 凭证键名格式

- API 源：`source_api::{workspaceId}::{sourceSlug}`
- MCP 源：`workspace_oauth::{workspaceId}::{sourceSlug}`

### 安全注意事项

1. 凭证文件使用强加密存储
2. 每个工作区的凭证相互隔离
3. 不要在 config.json 中硬编码敏感信息
4. 定期轮换 API 密钥

## 常见问题

### Q: 如何添加新的数据源？

1. 在 UI 中进入「数据源」设置
2. 点击「添加数据源」
3. 选择数据源类型并填写配置
4. 配置认证信息并测试连接

### Q: 预设数据源没有出现在工作区？

检查以下几点：
- 应用的 `sources/` 目录是否存在
- config.json 格式是否正确
- 工作区初始化时是否跳过了预设数据

### Q: 如何更新数据源配置？

直接编辑工作区下的 `sources/{slug}/config.json` 文件，或通过 UI 修改。修改后可能需要重新测试连接。

### Q: API 认证失败怎么办？

1. 检查 API Key 或 Token 是否正确
2. 确认认证类型配置正确
3. 查看 `testEndpoint` 是否可访问
4. 检查 API 服务的速率限制

## 参考资料

- [MCP 协议规范](https://modelcontextprotocol.io/)
- [Sprouty AI 架构文档](./architecture.md)
- [应用开发指南](./app-development.md)
