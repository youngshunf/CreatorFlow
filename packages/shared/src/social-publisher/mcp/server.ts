/**
 * social-publisher MCP Server 入口
 *
 * stdio 传输，注册 8 个工具
 * 只负责浏览器平台交互，不操作数据库
 */

import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  ListToolsRequestSchema,
  CallToolRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import { BrowserManager } from '../browser/browser-manager.ts';
import { FingerprintManager } from '../browser/fingerprint.ts';
import { getAdapter } from '../adapters/registry.ts';
import type { PublishContent, BrowserFingerprint, SocialPlatform } from '../types.ts';

/** 从环境变量获取工作区路径 */
const workspacePath = process.env.WORKSPACE_PATH || process.cwd();

/** 浏览器管理器实例 */
const browserManager = new BrowserManager(workspacePath);

/** 指纹管理器实例 */
const fingerprintManager = new FingerprintManager();

/** 创建 MCP Server */
const server = new Server(
  { name: 'social-publisher', version: '1.0.0' },
  { capabilities: { tools: {} } },
);

// ============================================================
// 工具定义
// ============================================================

server.setRequestHandler(ListToolsRequestSchema, async () => ({
  tools: [
    {
      name: 'social_login',
      description: '打开反检测浏览器窗口，让用户手动登录社交平台。登录成功后自动保存 Profile。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: "平台标识: 'xiaohongshu' | 'douyin' | 'wechat' | 'bilibili' | 'zhihu' | 'weibo' | 'x' | 'toutiao' | 'sina' | 'sohu'" },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置（可选，无则自动生成）' },
        },
        required: ['platform', 'profile_path'],
      },
    },
    {
      name: 'social_check_login',
      description: '使用无头浏览器检查指定账号的登录状态是否有效。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置' },
        },
        required: ['platform', 'profile_path'],
      },
    },
    {
      name: 'social_generate_fingerprint',
      description: '为指定平台生成一致性浏览器指纹配置。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
        },
        required: ['platform'],
      },
    },
    {
      name: 'social_publish',
      description: '通过反检测浏览器将内容发布到指定平台。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置' },
          content: {
            type: 'object',
            description: '发布内容',
            properties: {
              type: { type: 'string', description: "'image_text' | 'video' | 'article' | 'text'" },
              title: { type: 'string', description: '标题' },
              text: { type: 'string', description: '正文' },
              media_files: { type: 'array', items: { type: 'string' }, description: '本地文件绝对路径' },
              tags: { type: 'array', items: { type: 'string' }, description: '标签' },
              cover_image: { type: 'string', description: '封面图路径' },
            },
            required: ['type', 'text', 'media_files', 'tags'],
          },
        },
        required: ['platform', 'profile_path', 'content'],
      },
    },
    {
      name: 'social_publish_via_api',
      description: '通过平台 API 发布内容（适用于有 API 的平台，如公众号、X）。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          api_config: {
            type: 'object',
            description: 'API 配置',
            properties: {
              method: { type: 'string', description: "'wechat_api' | 'x_api'" },
              credentials: { type: 'object', description: 'API 凭证' },
            },
            required: ['method', 'credentials'],
          },
          content: {
            type: 'object',
            description: '发布内容',
            properties: {
              type: { type: 'string' },
              title: { type: 'string' },
              text: { type: 'string' },
              media_files: { type: 'array', items: { type: 'string' } },
              tags: { type: 'array', items: { type: 'string' } },
            },
            required: ['type', 'text', 'media_files', 'tags'],
          },
        },
        required: ['platform', 'api_config', 'content'],
      },
    },
    {
      name: 'social_scrape_metrics',
      description: '通过浏览器抓取已发布内容的数据指标。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置' },
          publish_url: { type: 'string', description: '发布内容的 URL' },
        },
        required: ['platform', 'profile_path', 'publish_url'],
      },
    },
    {
      name: 'social_scrape_comments',
      description: '通过浏览器抓取已发布内容的评论列表。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置' },
          publish_url: { type: 'string', description: '发布内容的 URL' },
          limit: { type: 'number', description: '最大评论数，默认 20' },
        },
        required: ['platform', 'profile_path', 'publish_url'],
      },
    },
    {
      name: 'social_reply_comment',
      description: '通过浏览器回复指定评论。',
      inputSchema: {
        type: 'object' as const,
        properties: {
          platform: { type: 'string', description: '平台标识' },
          profile_path: { type: 'string', description: '浏览器 Profile 目录路径' },
          fingerprint: { type: 'object', description: '浏览器指纹配置' },
          publish_url: { type: 'string', description: '发布内容的 URL' },
          comment_id: { type: 'string', description: '评论 ID' },
          reply_text: { type: 'string', description: '回复内容' },
        },
        required: ['platform', 'profile_path', 'publish_url', 'comment_id', 'reply_text'],
      },
    },
  ],
}));

// ============================================================
// 工具处理
// ============================================================

/** 构造 MCP 工具返回值 */
function toolResult(data: unknown) {
  return { content: [{ type: 'text' as const, text: JSON.stringify(data, null, 2) }] };
}

/** 构造错误返回值 */
function toolError(message: string) {
  return { content: [{ type: 'text' as const, text: JSON.stringify({ error: message }) }], isError: true };
}

server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;

  try {
    switch (name) {
      // ---- 登录管理 ----
      case 'social_login': {
        const { platform, fingerprint } = args as {
          platform: string; fingerprint?: BrowserFingerprint;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }
        const result = await browserManager.launchForLogin(adapter, fingerprint);
        return toolResult({ success: true, profile_saved: true });
      }

      case 'social_check_login': {
        const { platform, fingerprint } = args as {
          platform: string; fingerprint?: BrowserFingerprint;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }

        const context = await browserManager.launchHeadless(platform as SocialPlatform, fingerprint);
        try {
          const page = await context.newPage();
          await page.goto(adapter.getLoginUrl(), { waitUntil: 'domcontentloaded' });
          const loggedIn = await adapter.detectLoginSuccess(page);
          return toolResult({ logged_in: loggedIn });
        } finally {
          await browserManager.close();
        }
      }

      case 'social_generate_fingerprint': {
        const { platform } = args as { platform: string };
        const fingerprint = fingerprintManager.generateFingerprint(platform as SocialPlatform);
        return toolResult({ fingerprint });
      }

      // ---- 发布 ----
      case 'social_publish': {
        const { platform, fingerprint, content } = args as {
          platform: string; fingerprint?: BrowserFingerprint; content: PublishContent;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }

        const context = await browserManager.launchHeadless(platform as SocialPlatform, fingerprint);
        try {
          const result = await (adapter as any).publish(context, content);
          return toolResult(result);
        } finally {
          await browserManager.close();
        }
      }

      case 'social_publish_via_api': {
        const { platform, api_config, content } = args as {
          platform: string; api_config: { method: string; credentials: unknown }; content: PublishContent;
        };
        // API 发布由外部 MCP（baoyu-post-to-wechat 等）处理
        // 这里作为统一入口的预留
        return toolError(`API 发布请使用对应的平台 MCP（如 baoyu-post-to-wechat）。平台: ${platform}, 方法: ${api_config.method}`);
      }

      // ---- 数据抓取 ----
      case 'social_scrape_metrics': {
        const { platform, fingerprint, publish_url } = args as {
          platform: string; fingerprint?: BrowserFingerprint; publish_url: string;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }

        const context = await browserManager.launchHeadless(platform as SocialPlatform, fingerprint);
        try {
          const result = await (adapter as any).scrape(context, publish_url);
          return toolResult(result);
        } finally {
          await browserManager.close();
        }
      }

      case 'social_scrape_comments': {
        const { platform, fingerprint, publish_url, limit = 20 } = args as {
          platform: string; fingerprint?: BrowserFingerprint;
          publish_url: string; limit?: number;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }

        const context = await browserManager.launchHeadless(platform as SocialPlatform, fingerprint);
        try {
          const comments = await (adapter as any).fetchComments(context, publish_url, limit);
          return toolResult({ success: true, comments });
        } finally {
          await browserManager.close();
        }
      }

      // ---- 评论操作 ----
      case 'social_reply_comment': {
        const { platform, fingerprint, publish_url, comment_id, reply_text } = args as {
          platform: string; fingerprint?: BrowserFingerprint;
          publish_url: string; comment_id: string; reply_text: string;
        };
        const adapter = getAdapter(platform);
        if (!adapter) {
          return toolError(`不支持的平台: ${platform}`);
        }

        const context = await browserManager.launchHeadless(platform as SocialPlatform, fingerprint);
        try {
          const page = await context.newPage();
          await page.goto(publish_url, { waitUntil: 'domcontentloaded' });
          const result = await adapter.replyComment(page, comment_id, reply_text);
          return toolResult(result);
        } finally {
          await browserManager.close();
        }
      }

      default:
        return toolError(`未知工具: ${name}`);
    }
  } catch (error) {
    return toolError(`工具执行失败: ${error instanceof Error ? error.message : String(error)}`);
  }
});

// ============================================================
// 启动 Server
// ============================================================

async function main() {
  const transport = new StdioServerTransport();
  await server.connect(transport);
}

main().catch((error) => {
  console.error('social-publisher MCP Server 启动失败:', error);
  process.exit(1);
});
