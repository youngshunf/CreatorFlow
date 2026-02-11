---
name: find-skills
description: 帮助用户发现和安装代理技能 agent skills，当用户询问"如何做X"、"找一个X的技能"、"有没有能...的技能"，或表达对扩展功能的兴趣时使用。此技能应在用户寻找可能作为可安装技能存在的功能时使用。
---

# 查找技能

此技能帮助您从开放代理技能生态系统中发现和安装技能。

## 何时使用此技能

在以下情况下使用此技能：

- 用户询问"如何做X"，其中X可能是一个有现成技能的常见任务
- 用户说"找一个X的技能"或"有没有X的技能"
- 用户询问"你能做X吗"，其中X是一个专业能力
- 用户表达对扩展代理能力的兴趣
- 用户想要搜索工具、模板或工作流
- 用户提到希望在特定领域（设计、测试、部署等）获得帮助

## 什么是技能 CLI？

技能 CLI（`npx skills`）是开放代理技能生态系统的包管理器。技能是模块化包，通过专业知识、工作流和工具扩展代理能力。

**关键命令：**

- `npx skills find [query]` - 交互式搜索技能或按关键词搜索
- `npx skills add <package>` - 从 GitHub 或其他来源安装技能
- `npx skills check` - 检查技能更新
- `npx skills update` - 更新所有已安装的技能

**浏览技能：** https://skills.sh/

## 如何帮助用户查找技能

### 步骤 1：了解他们的需求

当用户寻求帮助时，识别：

1. 领域（例如：React、测试、设计、部署）
2. 具体任务（例如：编写测试、创建动画、审查 PR）
3. 这是否是一个足够常见的任务，可能存在相应技能

### 步骤 2：搜索技能

使用相关查询运行查找命令：

```bash
npx skills find [query]
```

例如：

- 用户询问"如何让我的 React 应用更快？" → `npx skills find react performance`
- 用户询问"你能帮我审查 PR 吗？" → `npx skills find pr review`
- 用户询问"我需要创建变更日志" → `npx skills find changelog`

命令将返回如下结果：

```
Install with npx skills add <owner/repo@skill>

vercel-labs/agent-skills@vercel-react-best-practices
└ https://skills.sh/vercel-labs/agent-skills/vercel-react-best-practices
```

### 步骤 3：向用户展示选项

当您找到相关技能时，向用户展示：

1. 技能名称及其功能
2. 他们可以运行的安装命令
3. 在 skills.sh 上了解更多的链接

示例响应：

```
我找到了一个可能有帮助的技能！"vercel-react-best-practices"技能提供
来自 Vercel 工程团队的 React 和 Next.js 性能优化指南。

安装方法：
npx skills add vercel-labs/agent-skills@vercel-react-best-practices

了解更多：https://skills.sh/vercel-labs/agent-skills/vercel-react-best-practices
```

### 步骤 4：提供安装

如果用户想要继续，您可以为他们安装技能：

```bash
npx skills add <owner/repo@skill> -g -y
```

`-g` 标志表示全局安装（用户级别），`-y` 跳过确认提示。

## 常见技能类别

搜索时，考虑这些常见类别：

| 类别     | 示例查询                                 |
| -------- | ---------------------------------------- |
| Web 开发 | react, nextjs, typescript, css, tailwind |
| 测试     | testing, jest, playwright, e2e           |
| DevOps   | deploy, docker, kubernetes, ci-cd        |
| 文档     | docs, readme, changelog, api-docs        |
| 代码质量 | review, lint, refactor, best-practices   |
| 设计     | ui, ux, design-system, accessibility     |
| 效率工具 | workflow, automation, git                |

## 有效搜索的技巧

1. **使用具体关键词**："react testing"比仅仅"testing"更好
2. **尝试替代术语**：如果"deploy"不起作用，尝试"deployment"或"ci-cd"
3. **查看热门来源**：许多技能来自 `vercel-labs/agent-skills` 或 `ComposioHQ/awesome-claude-skills`

## 当没有找到技能时

如果不存在相关技能：

1. 承认没有找到现有技能
2. 提供使用您的通用能力直接帮助完成任务
3. 建议用户可以使用 `npx skills init` 创建自己的技能

示例：

```
我搜索了与"xyz"相关的技能，但没有找到任何匹配项。
我仍然可以直接帮助您完成此任务！您想让我继续吗？

如果这是您经常做的事情，您可以创建自己的技能：
npx skills init my-xyz-skill
```
