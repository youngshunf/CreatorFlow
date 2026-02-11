/**
 * EditPopover
 *
 * A popover with title, subtitle, and multiline textarea for editing settings.
 * Supports two modes:
 * - Legacy: Opens a new focused window with a chat session
 * - Inline: Executes mini agent inline within the popover using compact ChatDisplay
 */

import * as React from 'react'
import { useState, useRef, useEffect, useCallback, useMemo } from 'react'
import { GripHorizontal } from 'lucide-react'
import { motion, AnimatePresence } from 'motion/react' // motion used for backdrop only
import { Popover, PopoverTrigger, PopoverContent } from './popover'
import { Button } from './button'
import { cn } from '@/lib/utils'
import { useT } from '@/context/LocaleContext'
import { usePlatform } from '@sprouty-ai/ui'
import type { ContentBadge, Session, CreateSessionOptions } from '../../../shared/types'
import { useActiveWorkspace, useAppShellContext, useSession } from '@/context/AppShellContext'
import { useEscapeInterrupt } from '@/context/EscapeInterruptContext'
import { ChatDisplay } from '../app-shell/ChatDisplay'

/**
 * Context passed to the new chat session so the agent knows exactly
 * what is being edited and can execute quickly.
 *
 * Simplified structure: label for display, filePath for the agent to know
 * where to edit, and optional context for additional instructions.
 */
export interface EditContext {
  /** Human-readable label for badge display and agent context (e.g., "Permissions") */
  label: string
  /** Absolute path to the file being edited */
  filePath: string
  /** Optional additional context/instructions for the agent */
  context?: string
}

/* ============================================================================
 * EDIT CONTEXT REGISTRY - SINGLE SOURCE OF TRUTH
 * ============================================================================
 * ALL edit contexts MUST be defined here. This is the canonical location.
 *
 * DO NOT create EditContext objects inline elsewhere in the codebase.
 * Instead, use getEditConfig() exported from this file.
 *
 * To add a new edit context:
 * 1. Add a new key to EditContextKey type
 * 2. Add the config to EDIT_CONFIGS
 * 3. Use via getEditConfig(key, location)
 *
 * This pattern ensures:
 * - All edit prompts and examples are reviewed in one place
 * - Consistent messaging to the agent
 * - Easy updates when context format changes
 * ============================================================================ */

/** Available edit context keys - add new ones here */
export type EditContextKey =
  | 'workspace-permissions'
  | 'default-permissions'
  | 'skill-instructions'
  | 'skill-metadata'
  | 'source-guide'
  | 'source-config'
  | 'source-permissions'
  | 'source-tool-permissions'
  | 'preferences-notes'
  | 'add-source'
  | 'add-source-api'   // Filter-specific: user is viewing APIs
  | 'add-source-mcp'   // Filter-specific: user is viewing MCPs
  | 'add-source-local' // Filter-specific: user is viewing Local Folders
  | 'add-skill'
  | 'edit-statuses'
  | 'edit-labels'
  | 'edit-auto-rules'
  | 'add-label'
  | 'edit-views'
  | 'edit-tool-icons'
  | 'creator-media-create-project'
  | 'creator-media-edit-profile'
  | 'creator-media-create-content'
  | 'creator-media-edit-scheduled-task'

/**
 * Full edit configuration including context for agent and example for UI.
 * Returned by getEditConfig() for use in EditPopover.
 */
export interface EditConfig {
  /** Context passed to the agent */
  context: EditContext
  /** Example text shown in the popover placeholder */
  example: string
  /** Optional custom placeholder text - overrides the default "Describe what you'd like to change" */
  overridePlaceholder?: string
  /** Optional model for mini agent (e.g., 'haiku', 'sonnet') */
  model?: string
  /** Optional system prompt preset for mini agent (e.g., 'mini' for focused edits) */
  systemPromptPreset?: 'default' | 'mini'
  /** When true, executes inline within the popover instead of opening a new window */
  inlineExecution?: boolean
}

/**
 * Registry of all edit configurations.
 * Each entry contains all strings needed for the edit popover and agent context.
 */
const EDIT_CONFIGS: Record<EditContextKey, (location: string) => EditConfig> = {
  'workspace-permissions': (location) => ({
    context: {
      label: '权限设置',
      filePath: `${location}/.sprouty-ai/permissions.json`,
      context:
        '用户在设置页面点击了工作区权限设置的编辑按钮。' +
        '除非另有说明，用户的意图很可能是立即更新设置。' +
        'permissions.json 文件配置探索模式规则，可包含：allowedBashPatterns、' +
        'allowedMcpPatterns、allowedApiEndpoints、blockedTools 和 allowedWritePaths。' +
        '编辑后，调用 config_validate（target 为 "permissions"）验证更改。' +
        '完成后明确确认。',
    },
    example: '在探索模式中允许运行 make build',
    model: 'sonnet',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  'default-permissions': (location) => ({
    context: {
      label: '默认权限',
      filePath: location,
      context:
        '用户正在编辑应用级默认权限（~/.sprouty-ai/permissions/default.json）。' +
        '此文件配置适用于所有工作区的探索模式规则。' +
        '可包含：allowedBashPatterns、allowedMcpPatterns、allowedApiEndpoints、blockedTools 和 allowedWritePaths。' +
        '每个模式可以是字符串或包含 pattern 和 comment 字段的对象。' +
        '请谨慎操作——这些是应用级全局默认值。' +
        '编辑后，调用 config_validate（target 为 "permissions"）验证更改。' +
        '完成后明确确认。',
    },
    example: '允许 git fetch 命令',
    model: 'sonnet',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Skill editing contexts
  'skill-instructions': (location) => ({
    context: {
      label: '技能指令',
      filePath: `${location}/SKILL.md`,
      context:
        '用户正在编辑 SKILL.md 中的技能指令。' +
        '重要：保留文件顶部的 YAML frontmatter（--- 标记之间的内容）。' +
        '专注于编辑 frontmatter 之后的 markdown 内容。' +
        '技能指令指导 AI 如何使用此技能。' +
        '编辑后，调用 skill_validate 并传入技能 slug 验证更改。' +
        '完成后明确确认。',
    },
    example: '添加错误处理指南',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  'skill-metadata': (location) => ({
    context: {
      label: '技能元数据',
      filePath: `${location}/SKILL.md`,
      context:
        '用户正在编辑 SKILL.md 的 YAML frontmatter 中的技能元数据。' +
        'Frontmatter 字段：name（必填）、description（必填）、globs（可选数组）、alwaysAllow（可选数组）。' +
        '除非特别要求，保持 frontmatter 之后的内容不变。' +
        '编辑后，调用 skill_validate 并传入技能 slug 验证更改。' +
        '完成后明确确认。',
    },
    example: '更新技能描述',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Source editing contexts
  'source-guide': (location) => ({
    context: {
      label: '数据源文档',
      filePath: `${location}/guide.md`,
      context:
        '用户正在编辑数据源文档（guide.md）。' +
        '此文件为 AI 提供关于如何使用此数据源的上下文——速率限制、API 模式、最佳实践等。' +
        '保持内容清晰且可操作。' +
        '完成后明确确认。',
    },
    example: '补充限流文档',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  'source-config': (location) => ({
    context: {
      label: '数据源配置',
      filePath: `${location}/config.json`,
      context:
        '用户正在编辑数据源配置（config.json）。' +
        '注意 JSON 语法。字段包括：type、slug、name、tagline、iconUrl，以及传输相关设置（mcp、api、local）。' +
        '除非明确要求，不要修改 slug。' +
        '编辑后，调用 source_test 并传入数据源 slug 验证配置。' +
        '完成后明确确认。',
    },
    example: '更新显示名称',
    model: 'sonnet',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  'source-permissions': (location) => ({
    context: {
      label: '数据源权限',
      filePath: `${location}/permissions.json`,
      context:
        '用户正在编辑数据源级权限（permissions.json）。' +
        '这些规则自动限定在此数据源范围内——编写简单模式即可，无需前缀。' +
        'MCP 类型：使用 allowedMcpPatterns（如 "list"、"get"）。API 类型：使用 allowedApiEndpoints。' +
        '编辑后，调用 config_validate（target 为 "permissions"）并传入数据源 slug 验证更改。' +
        '完成后明确确认。',
    },
    example: '探索模式允许 list 操作',
    model: 'sonnet',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  'source-tool-permissions': (location) => ({
    context: {
      label: '工具权限',
      filePath: `${location}/permissions.json`,
      context:
        '用户正在查看 MCP 数据源的工具列表，想要修改工具权限。' +
        '编辑 permissions.json 文件来控制探索模式下允许使用的工具。' +
        '使用 allowedMcpPatterns 允许特定工具（如 ["list_*", "get_*"] 表示只读）。' +
        '使用 blockedTools 明确阻止特定工具。' +
        '模式自动限定在此数据源范围内。' +
        '编辑后，调用 config_validate（target 为 "permissions"）并传入数据源 slug 验证更改。' +
        '完成后明确确认。',
    },
    example: '仅允许只读操作（list、get、search）',
    model: 'sonnet',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Preferences editing context
  'preferences-notes': (location) => ({
    context: {
      label: '偏好备注',
      filePath: location,
      context:
        '用户正在编辑偏好设置中的备注字段（~/.sprouty-ai/preferences.json）。' +
        '这是一个 JSON 文件。除非明确要求，只修改 "notes" 字段。' +
        'notes 字段是自由格式文本，为 AI 提供关于用户的上下文信息。' +
        '编辑后，调用 config_validate（target 为 "preferences"）验证更改。' +
        '完成后明确确认。',
    },
    example: '添加编码风格偏好',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Add new source/skill contexts - use overridePlaceholder for inspiring, contextual prompts
  'add-source': (location) => ({
    context: {
      label: '添加数据源',
      filePath: `${location}/.sprouty-ai/sources/`,
      context:
        '用户想要向工作区添加新的数据源。' +
        '数据源可以是 MCP 服务器（HTTP/SSE 或 stdio）、REST API 或本地文件系统。' +
        '如有需要可询问澄清问题：什么服务？MCP 还是 API？认证类型？' +
        '在工作区 .sprouty-ai/sources/ 目录中创建数据源文件夹和 config.json。' +
        '参照 ~/.sprouty-ai/docs/sources.md 中的模式。' +
        '创建后，调用 source_test 并传入数据源 slug 验证配置。',
    },
    example: '连接我的 Craft 空间',
    overridePlaceholder: '你想连接什么？',
  }),

  // Filter-specific add-source contexts: user is viewing a filtered list and wants to add that type
  'add-source-api': (location) => ({
    context: {
      label: '添加 API',
      filePath: `${location}/.sprouty-ai/sources/`,
      context:
        '用户正在查看 API 数据源列表，想要添加新的 REST API。' +
        '默认创建 API 数据源（type: "api"），除非用户另有说明。' +
        'API 通过 REST 端点连接，支持认证方式：bearer、header、basic 或 query。' +
        '询问 API 端点 URL 和认证类型。' +
        '在工作区 .sprouty-ai/sources/ 目录中创建数据源文件夹和 config.json。' +
        '参照 ~/.sprouty-ai/docs/sources.md 中的模式。' +
        '创建后，调用 source_test 并传入数据源 slug 验证配置。',
    },
    example: '连接 OpenAI API',
    overridePlaceholder: '你想连接哪个 API？',
  }),

  'add-source-mcp': (location) => ({
    context: {
      label: '添加 MCP 服务器',
      filePath: `${location}/.sprouty-ai/sources/`,
      context:
        '用户正在查看 MCP 数据源列表，想要添加新的 MCP 服务器。' +
        '默认创建 MCP 数据源（type: "mcp"），除非用户另有说明。' +
        'MCP 服务器可使用 HTTP/SSE 传输（远程）或 stdio 传输（本地子进程）。' +
        '询问用户想连接的服务以及是远程 URL 还是本地命令。' +
        '在工作区 .sprouty-ai/sources/ 目录中创建数据源文件夹和 config.json。' +
        '参照 ~/.sprouty-ai/docs/sources.md 中的模式。' +
        '创建后，调用 source_test 并传入数据源 slug 验证配置。',
    },
    example: '连接 Linear',
    overridePlaceholder: '你想连接哪个 MCP 服务器？',
  }),

  'add-source-local': (location) => ({
    context: {
      label: '添加本地文件夹',
      filePath: `${location}/.sprouty-ai/sources/`,
      context:
        '用户想要添加本地文件夹数据源。' +
        '首先查阅指南：mcp__creator-flows-docs__SearchSproutyAgents({ query: "filesystem" })。' +
        '本地文件夹是书签——使用 type: "local" 并设置 local.path 字段。' +
        '它们使用现有的 Read、Write、Glob、Grep 工具——不需要 MCP 服务器。' +
        '如不确定，询问用户想连接的文件夹路径。' +
        '在工作区 .sprouty-ai/sources/ 目录中创建数据源文件夹和 config.json。' +
        '参照 ~/.sprouty-ai/docs/sources.md 中的模式。' +
        '创建后，调用 source_test 并传入数据源 slug 验证配置。',
    },
    example: '连接我的 Obsidian 知识库',
    overridePlaceholder: '你想连接哪个本地文件夹？',
  }),

  'add-skill': (location) => ({
    context: {
      label: '添加技能',
      filePath: `${location}/.sprouty-ai/skills/`,
      context:
        '用户想要向工作区添加新技能。' +
        '技能是包含 SKILL.md 文件的专用指令，文件包含 YAML 前置元数据（name、description）和 Markdown 指令。' +
        '如有需要可询问澄清问题：技能做什么？何时触发？' +
        '在工作区 .sprouty-ai/skills/ 目录中创建技能文件夹和 SKILL.md。' +
        '参照 ~/.sprouty-ai/docs/skills.md 中的模式。' +
        '创建后，调用 skill_validate 并传入技能 slug 验证 SKILL.md 文件。',
    },
    example: '按照代码规范审查 PR',
    overridePlaceholder: '要让我学会做什么？',
  }),

  // Status configuration context
  'edit-statuses': (location) => ({
    context: {
      label: '状态配置',
      filePath: `${location}/.sprouty-ai/statuses/config.json`,
      context:
        '用户想要自定义会话状态（工作流状态）。' +
        '状态存储在 .sprouty-ai/statuses/config.json 中，字段包括：id、label、icon、category（open/closed）、order、isFixed、isDefault。' +
        '固定状态（todo、done、cancelled）不能删除，但可以重新排序或修改标签。' +
        'icon 可以是 { type: "file", value: "name.svg" }（自定义图标在 .sprouty-ai/statuses/icons/ 中）或 { type: "lucide", value: "icon-name" }（Lucide 图标）。' +
        'category 为 "open" 显示在收件箱，"closed" 显示在归档。' +
        '编辑后，调用 config_validate（target 为 "statuses"）验证更改。' +
        '完成后明确确认。',
    },
    example: '添加"阻塞"状态',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Label configuration context
  'edit-labels': (location) => ({
    context: {
      label: '标签配置',
      filePath: `${location}/.sprouty-ai/labels/config.json`,
      context:
        '用户想要自定义会话标签（分类/标记）。' +
        '标签存储在 .sprouty-ai/labels/config.json 中，以层级树结构组织。' +
        '每个标签包含：id（slug，全局唯一）、name（显示名称）、color（可选 EntityColor）、children（子标签数组）。' +
        '颜色使用 EntityColor 格式：字符串简写（如 "blue"）或 { light, dark } 对象用于主题感知颜色。' +
        '标签仅有颜色（无图标）——在 UI 中渲染为彩色圆点。' +
        'children 形成递归树结构——数组位置决定显示顺序。' +
        '参阅 ~/.sprouty-ai/docs/labels.md 获取完整格式参考。' +
        '完成后明确确认。',
    },
    example: '添加"缺陷"标签（红色）',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Auto-label rules context (focused on regex patterns within labels)
  'edit-auto-rules': (location) => ({
    context: {
      label: '自动应用规则',
      filePath: `${location}/.sprouty-ai/labels/config.json`,
      context:
        '用户想要编辑自动应用规则（自动标记会话的正则表达式模式）。' +
        '规则位于 .sprouty-ai/labels/config.json 中各标签的 autoRules 数组内。' +
        '每条规则包含：pattern（带捕获组的正则）、flags（默认 "gi"）、valueTemplate（$1/$2 替换）、description。' +
        '同一标签上的多条规则 = 多种触发方式。"g" 标志始终强制启用。' +
        '避免灾难性回溯模式（如 (a+)+）。' +
        '参阅 ~/.sprouty-ai/docs/labels.md 获取完整格式参考。' +
        '完成后明确确认。',
    },
    example: '添加检测 GitHub Issue 链接的规则',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Add new label context (triggered from the # menu when no labels match)
  'add-label': (location) => ({
    context: {
      label: '添加标签',
      filePath: `${location}/.sprouty-ai/labels/config.json`,
      context:
        '用户想要从 # 内联菜单创建新标签。' +
        '标签存储在 .sprouty-ai/labels/config.json 中，以层级树结构组织。' +
        '每个标签包含：id（slug，全局唯一）、name（显示名称）、color（可选 EntityColor）、children（子标签数组）。' +
        '颜色使用 EntityColor 格式：字符串简写（如 "blue"）或 { light, dark } 对象用于主题感知颜色。' +
        '标签仅有颜色（无图标）——在 UI 中渲染为彩色圆点。' +
        '参阅 ~/.sprouty-ai/docs/labels.md 获取完整格式参考。' +
        '完成后明确确认。',
    },
    example: '一个红色的"缺陷"标签',
    overridePlaceholder: '你想创建什么标签？',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Views configuration context
  'edit-views': (location) => ({
    context: {
      label: '视图配置',
      filePath: `${location}/.sprouty-ai/views.json`,
      context:
        '用户想要编辑视图（基于表达式的动态过滤器）。' +
        '视图存储在工作区 .sprouty-ai/views.json 中的 "views" 数组内。' +
        '每个视图包含：id（唯一 slug）、name（显示文本）、description（可选）、color（可选 EntityColor）、expression（Filtrex 字符串）。' +
        '表达式基于会话上下文字段求值：name、preview、todoState、permissionMode、model、lastMessageRole、' +
        'lastUsedAt、createdAt、messageCount、labelCount、isFlagged、hasUnread、isProcessing、hasPendingPlan、tokenUsage.*、labels。' +
        '可用函数：daysSince(timestamp)、contains(array, value)。' +
        '颜色使用 EntityColor 格式：字符串简写（如 "orange"）或 { light, dark } 对象。' +
        '完成后明确确认。',
    },
    example: '添加"过期"视图（7 天未活跃）',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Tool icons configuration context
  'edit-tool-icons': (location) => ({
    context: {
      label: '工具图标',
      filePath: location,
      context:
        '用户想要编辑 CLI 工具图标映射。' +
        '文件是 ~/.sprouty-ai/tool-icons/ 中的 tool-icons.json。图标图片文件在同一目录中。' +
        'Schema：{ version: 1, tools: [{ id, displayName, icon, commands }] }。' +
        '每个工具包含：id（唯一 slug）、displayName（UI 中显示）、icon（文件名如 "git.ico"）、commands（CLI 命令名数组）。' +
        '支持的图标格式：.png、.ico、.svg、.jpg。图标显示尺寸为 20x20px。' +
        '参阅 ~/.sprouty-ai/docs/tool-icons.md 获取完整格式参考。' +
        '编辑后，调用 config_validate（target 为 "tool-icons"）验证更改。' +
        '完成后明确确认。',
    },
    example: '为我的自定义 CLI 工具 "deploy" 添加图标',
    model: 'haiku',
    systemPromptPreset: 'mini',
    inlineExecution: true,
  }),

  // Creator-media contexts
  'creator-media-create-project': (location) => ({
    context: {
      label: '创建项目',
      filePath: `${location}/.sprouty-ai/db/creator.db`,
      context:
        '用户想要通过 AI 智能创建一个自媒体项目。请根据用户提供的信息（可能是账号链接、账号描述、或领域描述），分析并生成完整的项目信息和账号画像。\n\n' +
        '你需要：\n' +
        '1. 如果用户提供了账号链接，尝试访问并分析该账号的公开信息\n' +
        '2. 根据分析结果，使用 sqlite3 命令直接将项目和画像写入数据库\n' +
        '3. 所有字段都必须填写，不能为 NULL —— 如果用户未提供某些信息，请根据上下文智能推断合理的值\n\n' +
        '数据库路径：.sprouty-ai/db/creator.db（相对于工作目录）\n\n' +
        '操作步骤（必须严格按照以下顺序执行）：\n\n' +
        '步骤1 - 生成 UUID 和头像：\n' +
        'PROJECT_ID=$(uuidgen | tr \'[:upper:]\' \'[:lower:]\')\n' +
        'PROFILE_ID=$(uuidgen | tr \'[:upper:]\' \'[:lower:]\')\n\n' +
        '头像处理（avatar_path）：\n' +
        '- 如果用户提供了图片 URL，直接使用该 URL 作为 avatar_path\n' +
        '- 否则，根据项目名称和领域生成一个 SVG 头像文件，保存到 .sprouty-ai/avatars/<PROJECT_ID>.svg，avatar_path 填写该相对路径\n' +
        '- SVG 头像要求：简洁美观，使用领域相关图标，配合适合该领域的渐变背景色，尺寸 128x128\n\n' +
        '步骤2 - 创建项目并设为活跃：\n' +
        'sqlite3 .sprouty-ai/db/creator.db "UPDATE projects SET is_active = 0, updated_at = CURRENT_TIMESTAMP WHERE is_active = 1;"\n' +
        'sqlite3 .sprouty-ai/db/creator.db "INSERT INTO projects (id, name, description, platform, platforms, avatar_path, is_active) VALUES (...)"\n\n' +
        '步骤3 - 创建账号画像（所有字段必填）：\n' +
        'sqlite3 .sprouty-ai/db/creator.db "INSERT INTO account_profiles (id, project_id, niche, sub_niche, persona, target_audience, tone, keywords, bio, content_pillars, posting_frequency, best_posting_time, style_references, taboo_topics, pillar_weights) VALUES (...)"\n\n' +
        '字段说明：\n' +
        '- platform 取值：xiaohongshu / douyin / bilibili / weibo / wechat / zhihu / toutiao / kuaishou / youtube / tiktok / instagram / twitter / threads / facebook / linkedin / pinterest / medium / substack / other\n' +
        '- platforms：JSON 数组字符串，如 \'["xiaohongshu","douyin"]\'\n' +
        '- posting_frequency 取值：daily / 3_per_week / weekly / biweekly / monthly\n' +
        '- best_posting_time：如 \'20:00\' 或 \'12:00,20:00\'\n' +
        '- pillar_weights：JSON 对象，键为内容支柱名称，值为权重(0-1)，所有权重之和为1\n' +
        '- 文本中的单引号需要转义为两个单引号（SQL 标准）\n\n' +
        '重要：\n' +
        '- 所有字段都必须有值，不允许 NULL，请根据用户信息和领域知识智能推断\n' +
        '- 必须使用 sqlite3 命令操作数据库，不要调用任何 MCP 工具\n' +
        '- 创建 SVG 头像前先确保目录存在：mkdir -p .sprouty-ai/avatars\n' +
        '- 创建完成后告知用户项目已创建成功',
    },
    example: '我是一个小红书美妆博主，主要做护肤品测评',
    overridePlaceholder: '描述你的账号或粘贴主页链接',
    model: 'sonnet',
    systemPromptPreset: 'default',
    inlineExecution: true,
  }),

  'creator-media-edit-profile': (location) => ({
    context: {
      label: '编辑画像',
      filePath: `${location}/.sprouty-ai/db/creator.db`,
      context:
        '用户想要修改当前活跃项目的账号画像。\n\n' +
        '操作步骤：\n' +
        '1. 先读取当前活跃项目和画像：\n' +
        '   sqlite3 -header -column .sprouty-ai/db/creator.db "SELECT p.id, p.name, ap.* FROM projects p LEFT JOIN account_profiles ap ON ap.project_id = p.id WHERE p.is_active = 1;"\n' +
        '2. 根据用户需求，使用 sqlite3 UPDATE 语句修改 account_profiles 表中的对应字段\n' +
        '3. 如果画像不存在，使用 INSERT 创建\n\n' +
        '数据库路径：.sprouty-ai/db/creator.db（相对于工作目录）\n\n' +
        '可修改字段：niche, sub_niche, persona, target_audience, tone, keywords, bio, content_pillars, posting_frequency, best_posting_time, style_references, taboo_topics, pillar_weights\n\n' +
        '重要：\n' +
        '- 必须使用 sqlite3 命令操作数据库\n' +
        '- 文本中的单引号需要转义为两个单引号\n' +
        '- 修改完成后告知用户已更新',
    },
    example: '把目标受众改为 25-35 岁职场女性',
    overridePlaceholder: '描述你想修改的画像内容',
    model: 'haiku',
    systemPromptPreset: 'default',
    inlineExecution: true,
  }),

  'creator-media-create-content': (location) => ({
    context: {
      label: '新建内容',
      filePath: `${location}/.sprouty-ai/db/creator.db`,
      context:
        '用户想要为当前活跃项目创建新的内容条目。\n\n' +
        '操作步骤：\n' +
        '1. 先读取当前活跃项目：\n' +
        '   sqlite3 .sprouty-ai/db/creator.db "SELECT id, name, platform FROM projects WHERE is_active = 1;"\n' +
        '2. 生成 UUID：CONTENT_ID=$(uuidgen | tr \'[:upper:]\' \'[:lower:]\')\n' +
        '3. 使用 sqlite3 INSERT 创建内容：\n' +
        '   sqlite3 .sprouty-ai/db/creator.db "INSERT INTO contents (id, project_id, title, content_type, status, body) VALUES (...)"\n\n' +
        '数据库路径：.sprouty-ai/db/creator.db（相对于工作目录）\n\n' +
        '字段说明：\n' +
        '- content_type 取值：article / short-video / video / image-text / thread / story\n' +
        '- status 取值：idea / creating / review / scheduled / published\n' +
        '- body：内容正文（Markdown 格式）\n' +
        '- 可选字段：hook, outline, tags, scheduled_at, metadata\n\n' +
        '重要：\n' +
        '- 必须使用 sqlite3 命令操作数据库\n' +
        '- 根据用户描述智能推断 content_type\n' +
        '- 创建完成后告知用户',
    },
    example: '写一篇关于防晒霜选购指南的小红书笔记',
    overridePlaceholder: '描述你想创建的内容',
    model: 'haiku',
    systemPromptPreset: 'default',
    inlineExecution: true,
  }),

  'creator-media-edit-scheduled-task': (location) => ({
    context: {
      label: '定时任务',
      filePath: `${location}/.sprouty-ai/db/creator.db`,
      context:
        '用户想要通过 AI 辅助管理定时任务。\n\n' +
        '操作步骤：\n' +
        '1. 先读取当前定时任务列表：\n' +
        '   sqlite3 -header -column .sprouty-ai/db/creator.db "SELECT id, name, task_type, schedule_mode, cron_expression, interval_seconds, enabled, status FROM scheduled_tasks ORDER BY created_at DESC;"\n' +
        '2. 根据用户需求，使用 sqlite3 INSERT/UPDATE/DELETE 操作 scheduled_tasks 表\n\n' +
        '数据库路径：.sprouty-ai/db/creator.db（相对于工作目录）\n\n' +
        '字段说明：\n' +
        '- task_type：review / publish / collect / custom\n' +
        '- schedule_mode：cron / interval / once\n' +
        '- cron_expression：Cron 表达式（当 schedule_mode=cron 时）\n' +
        '- interval_seconds：间隔秒数（当 schedule_mode=interval 时）\n' +
        '- scheduled_at：执行时间（当 schedule_mode=once 时）\n' +
        '- enabled：1=启用 0=禁用\n' +
        '- status：active / paused / error / completed\n' +
        '- payload：JSON 格式的任务配置\n\n' +
        '重要：\n' +
        '- 必须使用 sqlite3 命令操作数据库\n' +
        '- 新建任务时 id 使用 UUID\n' +
        '- 操作完成后告知用户',
    },
    example: '创建一个每天早上8点执行的采集任务',
    overridePlaceholder: '描述你想创建或修改的定时任务',
    model: 'haiku',
    systemPromptPreset: 'default',
    inlineExecution: true,
  }),
}

/**
 * Get full edit config by key. Returns both context (for agent) and example (for UI).
 *
 * @param key - The edit context key
 * @param location - Base path (e.g., workspace root path)
 *
 * @example
 * const { context, example } = getEditConfig('workspace-permissions', workspace.rootPath)
 */
export function getEditConfig(key: EditContextKey, location: string): EditConfig {
  const factory = EDIT_CONFIGS[key]
  if (!factory) {
    throw new Error(`Unknown edit context key: ${key}. Add it to EDIT_CONFIGS in EditPopover.tsx`)
  }
  return factory(location)
}

/**
 * Optional secondary action button displayed on the left side of the popover footer.
 * Styled as plain text with underline on hover - typically used for "Edit File" actions.
 */
export interface SecondaryAction {
  /** Button label (e.g., "Edit File") */
  label: string
  /** Click handler - typically opens a file for manual editing */
  onClick: () => void
}

export interface EditPopoverProps {
  /** Trigger element that opens the popover */
  trigger: React.ReactNode
  /** Example text shown in placeholder (e.g., "Allow 'make build' command") */
  example?: string
  /** Context passed to the new chat session */
  context: EditContext
  /** Permission mode for the new session (default: 'allow-all' for fast execution) */
  permissionMode?: 'safe' | 'ask' | 'allow-all'
  /**
   * Working directory for the new session:
   * - 'none' (default): No working directory (session folder only) - best for config edits
   * - 'user_default': Use workspace's configured default
   * - Absolute path string: Use this specific path
   */
  workingDirectory?: string | 'user_default' | 'none'
  /** Model override for mini agent (e.g., 'haiku', 'sonnet') */
  model?: string
  /** System prompt preset for mini agent (e.g., 'mini' for focused edits) */
  systemPromptPreset?: 'default' | 'mini'
  /** Width of the popover (default: 320) */
  width?: number
  /** Additional className for the trigger */
  triggerClassName?: string
  /** Side of the popover relative to trigger */
  side?: 'top' | 'right' | 'bottom' | 'left'
  /** Alignment of the popover */
  align?: 'start' | 'center' | 'end'
  /** Optional secondary action button on the left (e.g., "Edit File") */
  secondaryAction?: SecondaryAction
  /** Optional custom placeholder - overrides the default "Describe what you'd like to change" */
  overridePlaceholder?: string
  /**
   * Controlled open state - when provided, the popover becomes controlled.
   * Use this when opening the popover programmatically (e.g., from context menus).
   */
  open?: boolean
  /** Callback when open state changes (for controlled mode) */
  onOpenChange?: (open: boolean) => void
  /**
   * When true, prevents the popover from closing when clicking outside.
   * Useful for context menu triggered popovers where focus management is tricky.
   */
  modal?: boolean
  /**
   * Default value to pre-fill the input with.
   * Useful when the user types something (e.g., "#Test") and clicks "Add new label" -
   * the input can be pre-filled with "Add new label Test".
   */
  defaultValue?: string
  /**
   * When true, executes the mini agent inline within the popover instead of
   * opening a new window. Best for quick config edits with mini agents.
   */
  inlineExecution?: boolean
}

/**
 * Result from buildEditPrompt containing both the full prompt and badge metadata
 * for hiding the XML context in the UI while keeping it in the actual message.
 */
interface EditPromptResult {
  /** Full prompt including XML metadata and user instructions */
  prompt: string
  /** Badge marking the hidden metadata section */
  badges: ContentBadge[]
}

/**
 * Build the prompt that will be sent to the agent.
 * Uses XML-like tags for clear structure.
 *
 * Returns both the prompt and a context badge that marks the metadata section
 * so it can be hidden in the UI while still being sent to the agent.
 *
 * @param context - The edit context with label, filePath, and optional context
 * @param userInstructions - User's instructions (can be empty string for pre-filled context only)
 *
 * @example
 * // With user instructions (for EditPopover submit)
 * const { prompt, badges } = buildEditPrompt(context, "Add a Blocked status")
 *
 * // Without user instructions (for context menu - opens window with context pre-filled)
 * const { prompt, badges } = buildEditPrompt(context, "")
 */
export function buildEditPrompt(context: EditContext, userInstructions: string): EditPromptResult {
  // Build the metadata section (will be hidden by badge)
  // Simple structure: label (for display/context), file (where to edit), optional context
  const metadataSection = `<edit_request>
<label>${context.label}</label>
<file>${context.filePath}</file>
${context.context ? `<context>${context.context}</context>\n` : ''}</edit_request>

`

  // Badge display: just the label (no "Edit:" prefix for cleaner appearance)
  const collapsedLabel = context.label

  // Full prompt = metadata + user instructions
  const prompt = metadataSection + userInstructions

  // Create badge marking the metadata section (start=0, end=metadata length)
  const badge: ContentBadge = {
    type: 'context',
    label: collapsedLabel,
    rawText: metadataSection,
    start: 0,
    end: metadataSection.length,
    collapsedLabel,
  }

  return { prompt, badges: [badge] }
}

export function EditPopover({
  trigger,
  example,
  context,
  permissionMode = 'allow-all',
  workingDirectory = 'none', // Default to session folder for config edits
  model,
  systemPromptPreset,
  width = 400, // Default 400px for compact chat embedding
  triggerClassName,
  side = 'bottom',
  align = 'end',
  secondaryAction,
  overridePlaceholder,
  open: controlledOpen,
  onOpenChange: controlledOnOpenChange,
  modal = false,
  defaultValue = '',
  inlineExecution = false,
}: EditPopoverProps) {
  const t = useT()
  const { onOpenFile, onOpenUrl } = usePlatform()
  const workspace = useActiveWorkspace()

  // Build placeholder: for inline execution use rotating array, otherwise build descriptive string with i18n
  // overridePlaceholder allows contexts like add-source/add-skill to say "add" instead of "change"
  const COMPACT_PLACEHOLDERS = [
    t('告诉我您想改什么'),
    t('描述您的更新'),
    t('您想让我修改什么？'),
  ]
  const placeholder = inlineExecution
    ? COMPACT_PLACEHOLDERS
    : (() => {
        const basePlaceholder = overridePlaceholder ? t(overridePlaceholder) : t('描述您想要进行的更改...')
        const exampleText = example ? t(example) : undefined
        return exampleText
          ? `${basePlaceholder.replace(/\.{3}$/, '')}，${t('例如')}："${exampleText}"`
          : basePlaceholder
      })()

  // Support both controlled and uncontrolled modes:
  // - Uncontrolled (default): internal state manages open/close
  // - Controlled: parent manages state via open/onOpenChange props
  const [internalOpen, setInternalOpen] = useState(false)
  const isControlled = controlledOpen !== undefined
  const open = isControlled ? controlledOpen : internalOpen
  const setOpen = (value: boolean) => {
    if (isControlled) {
      controlledOnOpenChange?.(value)
    } else {
      setInternalOpen(value)
    }
  }

  // Use App context for session management (same code path as main chat)
  const { onCreateSession, onSendMessage } = useAppShellContext()

  // Session ID for inline execution (created on first message)
  const [inlineSessionId, setInlineSessionId] = useState<string | null>(null)

  // Get session data from Jotai atom (same as main chat - includes optimistic updates)
  // Pass empty string when no session yet - atom returns null for unknown IDs
  const inlineSession = useSession(inlineSessionId || '')

  // Model state for ChatDisplay (starts with prop value, can be changed by user)
  const [currentModel, setCurrentModel] = useState(model || 'haiku')

  // Create a stub session for ChatDisplay when no real session exists yet
  // This allows showing the input before the first message is sent
  const stubSession = useMemo((): Session => ({
    id: 'pending',
    workspaceId: workspace?.id || '',
    workspaceName: workspace?.name || '',
    messages: [],
    isProcessing: false,
    lastMessageAt: Date.now(),
  }), [workspace?.id, workspace?.name])

  // Use real session if available, otherwise stub
  const displaySession = inlineSession || stubSession

  // Track processing state for close prevention and backdrop
  const isProcessing = displaySession.isProcessing

  // Use existing escape interrupt context for double-ESC flow
  // This shows the "Press Esc again to interrupt" overlay in the input field
  const { handleEscapePress } = useEscapeInterrupt()

  // Reset inline session when popover closes
  const resetInlineSession = useCallback(() => {
    setInlineSessionId(null)
  }, [])

  // Stop/cancel generation for the inline session
  const handleStopGeneration = useCallback(() => {
    if (inlineSessionId && isProcessing) {
      window.electronAPI.cancelProcessing(inlineSessionId, false)
    }
  }, [inlineSessionId, isProcessing])

  // Handle ESC key during generation:
  // Uses EscapeInterruptContext for double-ESC flow (shows overlay, then interrupts)
  const handleEscapeKeyDown = useCallback((e: KeyboardEvent) => {
    if (!isProcessing) {
      // Not processing - allow normal close behavior
      return
    }

    // Prevent default close behavior during processing
    e.preventDefault()

    // Use context's double-ESC handler
    // Returns true if this is the second press (should interrupt)
    const shouldInterrupt = handleEscapePress()
    if (shouldInterrupt) {
      handleStopGeneration()
    }
  }, [isProcessing, handleEscapePress, handleStopGeneration])

  // Handle click outside during generation:
  // Show the ESC overlay via context, prevent closing
  const handleInteractOutside = useCallback((e: Event) => {
    if (isProcessing) {
      // Prevent close during processing
      e.preventDefault()
      // Show the ESC overlay so user knows how to cancel
      handleEscapePress()
    }
  }, [isProcessing, handleEscapePress])

  // Drag state for movable popover
  const [dragOffset, setDragOffset] = useState({ x: 0, y: 0 })
  const [isDragging, setIsDragging] = useState(false)
  const dragStartRef = useRef({ x: 0, y: 0, offsetX: 0, offsetY: 0 })
  const dragOffsetRef = useRef({ x: 0, y: 0 })
  const popoverRef = useRef<HTMLDivElement>(null)

  // Resize state for dynamic sizing
  const [containerSize, setContainerSize] = useState({ width: width || 400, height: 480 })
  const [isResizing, setIsResizing] = useState(false)
  const resizeStartRef = useRef({ x: 0, y: 0, width: 0, height: 0 })

  // Dynamic max dimensions based on viewport space
  const [maxDimensions, setMaxDimensions] = useState({
    maxWidth: width || 400,
    maxHeight: 480
  })

  // Calculate available space dynamically when popover opens
  useEffect(() => {
    if (!open) return

    const calculateMaxDimensions = () => {
      // Find the trigger element - use more specific selector
      // Radix adds data-state to the trigger when popover is open
      const triggerElement = document.querySelector('[data-slot="popover-trigger"][data-state="open"]')
      if (!triggerElement) {
        // Fallback: if we can't find the trigger, use generous defaults
        setMaxDimensions({ maxWidth: 600, maxHeight: 600 })
        return
      }

      const rect = triggerElement.getBoundingClientRect()
      const viewportWidth = window.innerWidth
      const viewportHeight = window.innerHeight

      // Calculate available space in each direction (留 10px padding)
      const availableRight = viewportWidth - rect.right - 10
      const availableLeft = rect.left - 10
      const availableBottom = viewportHeight - rect.bottom - 10
      const availableTop = rect.top - 10

      // Calculate max width based on align prop
      let maxWidth: number
      if (align === 'end') {
        // 右对齐：从按钮右边缘向左延伸
        maxWidth = Math.min(rect.right - 10, 600)
      } else if (align === 'start') {
        // 左对齐：从按钮左边缘向右延伸
        maxWidth = Math.min(viewportWidth - rect.left - 10, 600)
      } else {
        // 居中：取两侧较大值
        maxWidth = Math.min(Math.max(availableRight, availableLeft) * 2, 600)
      }

      // 确保最小宽度
      maxWidth = Math.max(maxWidth, 320)

      // Calculate max height based on side prop
      const maxHeight = side === 'bottom'
        ? Math.min(availableBottom, 600)
        : Math.min(availableTop, 600)

      setMaxDimensions({
        maxWidth,
        maxHeight: Math.max(maxHeight, 400)
      })
    }

    // 延迟计算，确保 DOM 已渲染
    const timer = setTimeout(calculateMaxDimensions, 0)
    return () => clearTimeout(timer)
  }, [open, side, align])

  // Reset drag position and size when popover opens
  useEffect(() => {
    if (open) {
      dragOffsetRef.current = { x: 0, y: 0 }
      setDragOffset({ x: 0, y: 0 })
      // 使用动态计算的最大尺寸
      setContainerSize({
        width: Math.min(width || 400, maxDimensions.maxWidth),
        height: Math.min(480, maxDimensions.maxHeight)
      })
    }
  }, [open, width, maxDimensions])

  // Handle drag events
  const handleDragStart = useCallback((e: React.MouseEvent) => {
    e.preventDefault()
    setIsDragging(true)
    dragStartRef.current = {
      x: e.clientX,
      y: e.clientY,
      offsetX: dragOffset.x,
      offsetY: dragOffset.y,
    }
  }, [dragOffset])

  useEffect(() => {
    if (!isDragging) return

    const handleMouseMove = (e: MouseEvent) => {
      const rect = popoverRef.current?.getBoundingClientRect()
      if (!rect) return

      const MARGIN = 20
      const MARGIN_TOP = 50 // Keep below header (chevrons, menu button)
      const curr = dragOffsetRef.current
      const baseX = rect.left - curr.x
      const baseY = rect.top - curr.y

      const newX = dragStartRef.current.offsetX + e.clientX - dragStartRef.current.x
      const newY = dragStartRef.current.offsetY + e.clientY - dragStartRef.current.y

      const clampedX = Math.max(MARGIN - baseX, Math.min(window.innerWidth - MARGIN - rect.width - baseX, newX))
      const clampedY = Math.max(MARGIN_TOP - baseY, Math.min(window.innerHeight - MARGIN - rect.height - baseY, newY))

      dragOffsetRef.current = { x: clampedX, y: clampedY }
      setDragOffset({ x: clampedX, y: clampedY })
    }

    const handleMouseUp = () => {
      setIsDragging(false)
    }

    document.addEventListener('mousemove', handleMouseMove)
    document.addEventListener('mouseup', handleMouseUp)
    return () => {
      document.removeEventListener('mousemove', handleMouseMove)
      document.removeEventListener('mouseup', handleMouseUp)
    }
  }, [isDragging])

  // Resize handlers
  const handleResizeStart = useCallback((e: React.MouseEvent) => {
    e.preventDefault()
    e.stopPropagation()
    setIsResizing(true)
    resizeStartRef.current = {
      x: e.clientX,
      y: e.clientY,
      width: containerSize.width,
      height: containerSize.height,
    }
  }, [containerSize])

  useEffect(() => {
    if (!isResizing) return

    const handleMouseMove = (e: MouseEvent) => {
      const deltaX = e.clientX - resizeStartRef.current.x
      const deltaY = e.clientY - resizeStartRef.current.y
      setContainerSize({
        width: Math.max(300, resizeStartRef.current.width + deltaX),
        height: Math.max(250, resizeStartRef.current.height + deltaY),
      })
    }

    const handleMouseUp = () => {
      setIsResizing(false)
    }

    document.addEventListener('mousemove', handleMouseMove)
    document.addEventListener('mouseup', handleMouseUp)
    return () => {
      document.removeEventListener('mousemove', handleMouseMove)
      document.removeEventListener('mouseup', handleMouseUp)
    }
  }, [isResizing])

  // Reset state when popover opens
  useEffect(() => {
    if (open) {
      setCurrentModel(model || 'haiku')
      resetInlineSession()
    }
  }, [open, model, resetInlineSession])

  // Handle sending message from ChatDisplay (inline mode)
  // Creates hidden session on first message, then uses App context for sending
  const handleInlineSendMessage = useCallback(async (message: string) => {
    const { prompt, badges } = buildEditPrompt(context, message)

    // Create session on first message
    let sessionId = inlineSessionId
    if (!sessionId && workspace?.id) {
      const createOptions: CreateSessionOptions = {
        model: model || 'haiku',
        systemPromptPreset: systemPromptPreset || 'mini',
        permissionMode,
        workingDirectory,
        hidden: true, // Hidden sessions use same App code path but don't appear in list
      }
      const newSession = await onCreateSession(workspace.id, createOptions)
      sessionId = newSession.id
      setInlineSessionId(sessionId)
    }

    // Send message via App context (includes optimistic user message update)
    // Pass badges to hide the <edit_request> XML metadata in the user message bubble
    if (sessionId) {
      onSendMessage(sessionId, prompt, undefined, undefined, badges)
    }
  }, [context, inlineSessionId, workspace?.id, model, systemPromptPreset, permissionMode, workingDirectory, onCreateSession, onSendMessage])

  // Legacy mode: navigates to chat in the same window
  const handleLegacySendMessage = useCallback((message: string) => {
    const { prompt, badges } = buildEditPrompt(context, message)
    const encodedInput = encodeURIComponent(prompt)
    const encodedBadges = encodeURIComponent(JSON.stringify(badges))

    const workdirParam = workingDirectory ? `&workdir=${encodeURIComponent(workingDirectory)}` : ''
    const modelParam = model ? `&model=${encodeURIComponent(model)}` : ''
    const systemPromptParam = systemPromptPreset ? `&systemPrompt=${encodeURIComponent(systemPromptPreset)}` : ''
    // Navigate in same window by omitting window=focused parameter
    const url = `sproutyai://action/new-chat?input=${encodedInput}&send=true&mode=${permissionMode}&badges=${encodedBadges}${workdirParam}${modelParam}${systemPromptParam}`

    window.electronAPI.openUrl(url)
    setOpen(false)
  }, [context, workingDirectory, model, systemPromptPreset, permissionMode, setOpen])

  return (
    <>
      {/* Full-screen backdrop - rendered BEHIND the popover during processing */}
      <AnimatePresence>
        {open && isProcessing && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.5, ease: 'easeInOut' }}
            className="fixed inset-0 bg-black/5 z-40"
          />
        )}
      </AnimatePresence>

      <Popover open={open} onOpenChange={setOpen} modal={modal}>
        <PopoverTrigger asChild className={triggerClassName}>
          {trigger}
        </PopoverTrigger>
        <PopoverContent
            side={side}
            align={align}
            className="p-0"
            style={{
              width: Math.min(containerSize.width, maxDimensions.maxWidth),
              height: Math.min(containerSize.height, maxDimensions.maxHeight),
              background: 'transparent',
              border: 'none',
              boxShadow: 'none',
            }}
            onInteractOutside={handleInteractOutside}
            onEscapeKeyDown={handleEscapeKeyDown}
          >
            {/* Container - size inherited from PopoverContent for Radix collision detection */}
            <div
              ref={popoverRef}
              className="relative bg-foreground-2 overflow-hidden w-full h-full"
              style={{
                width: '100%',
                height: '100%',
                transform: `translate(${dragOffset.x}px, ${dragOffset.y}px)`,
                borderRadius: 16,
                boxShadow: '0 4px 24px rgba(0, 0, 0, 0.12), 0 0 0 1px rgba(0, 0, 0, 0.05)',
              }}
            >
              {/* Drag handle - floating overlay */}
              <div
                onMouseDown={handleDragStart}
                className={cn(
                  "absolute top-0 left-1/2 -translate-x-1/2 z-50 px-4 py-2 cursor-grab rounded pointer-events-auto titlebar-no-drag",
                  isDragging && "cursor-grabbing"
                )}
              >
                <GripHorizontal className="w-4 h-4 text-muted-foreground/30" />
              </div>

              {/* Content area - always uses compact ChatDisplay */}
              <div className="flex-1 flex flex-col bg-foreground-2" style={{ height: '100%' }}>
                <ChatDisplay
                  session={displaySession}
                  onSendMessage={inlineExecution ? handleInlineSendMessage : handleLegacySendMessage}
                  onOpenFile={onOpenFile || (() => {})}
                  onOpenUrl={onOpenUrl || (() => {})}
                  currentModel={currentModel}
                  onModelChange={setCurrentModel}
                  compactMode={true}
                  placeholder={placeholder}
                  emptyStateLabel={context.label}
                />
              </div>

              {/* Bottom-right resize handle - invisible hit area */}
              <div
                onMouseDown={handleResizeStart}
                className="absolute -bottom-2 -right-2 w-6 h-6 cursor-nwse-resize pointer-events-auto z-50"
              />
            </div>
          </PopoverContent>
      </Popover>
    </>
  )
}

/**
 * Standard Edit button styled for use with EditPopover.
 * Use this as the trigger prop for consistent styling across the app.
 *
 * Uses forwardRef to properly work with Radix's asChild pattern,
 * which requires the child to accept ref and spread props.
 *
 * @example
 * <EditPopover
 *   trigger={<EditButton />}
 *   context={getEditContext('workspace-permissions', { workspacePath })}
 * />
 */
export const EditButton = React.forwardRef<
  HTMLButtonElement,
  React.ComponentPropsWithoutRef<typeof Button>
>(function EditButton({ className, ...props }, ref) {
  const t = useT()
  return (
    <Button
      ref={ref}
      variant="ghost"
      size="sm"
      // Merge our base styles with any className from asChild props
      className={cn("h-8 px-3 rounded-[6px] bg-background shadow-minimal text-foreground/70 hover:text-foreground", className)}
      {...props}
    >
      {t('编辑')}
    </Button>
  )
})
