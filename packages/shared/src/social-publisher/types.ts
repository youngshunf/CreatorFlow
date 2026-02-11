/**
 * social-publisher 类型定义
 *
 * 社交平台自动发布系统的核心类型，包括平台定义、浏览器指纹、
 * 发布内容、发布结果、数据指标、评论等。
 */

import type { Page, BrowserContext } from 'playwright';

// ============================================================
// 平台枚举
// ============================================================

/** 支持的社交平台 */
export type SocialPlatform =
  | 'xiaohongshu'
  | 'douyin'
  | 'wechat'
  | 'bilibili'
  | 'zhihu'
  | 'weibo'
  | 'x'
  | 'toutiao'
  | 'sina'
  | 'sohu';

// ============================================================
// 浏览器指纹
// ============================================================

/** 浏览器指纹配置 */
export interface BrowserFingerprint {
  /** 用户代理字符串 */
  userAgent: string;
  /** 视口尺寸 */
  viewport: { width: number; height: number };
  /** 语言区域 */
  locale: string;
  /** 时区标识 */
  timezone: string;
  /** WebGL 厂商 */
  webglVendor: string;
  /** WebGL 渲染器 */
  webglRenderer: string;
  /** 操作系统平台标识 */
  platform: string;
  /** CPU 核心数 */
  hardwareConcurrency: number;
  /** 设备内存（GB） */
  deviceMemory: number;
}

// ============================================================
// 发布内容
// ============================================================

/** 发布内容类型 */
export type PublishContentType = 'image_text' | 'video' | 'article' | 'text';

/** 发布内容 */
export interface PublishContent {
  /** 内容类型 */
  type: PublishContentType;
  /** 标题（部分平台必填） */
  title?: string;
  /** 正文内容 */
  text: string;
  /** 媒体文件路径列表 */
  media_files: string[];
  /** 标签列表 */
  tags: string[];
  /** 封面图路径 */
  cover_image?: string;
}

// ============================================================
// 发布结果
// ============================================================

/** 发布结果 */
export interface PublishResult {
  /** 是否成功 */
  success: boolean;
  /** 发布后的平台链接 */
  platform_url?: string;
  /** 平台侧帖子 ID */
  platform_post_id?: string;
  /** 错误信息 */
  error?: string;
}

// ============================================================
// 数据指标
// ============================================================

/** 数据指标详情 */
export interface MetricsData {
  /** 浏览量 */
  views: number;
  /** 点赞数 */
  likes: number;
  /** 评论数 */
  comments: number;
  /** 分享数 */
  shares: number;
  /** 收藏数 */
  favorites: number;
  /** 原始 JSON（平台特有字段） */
  raw_json?: string;
}

/** 数据指标采集结果 */
export interface MetricsResult {
  /** 是否成功 */
  success: boolean;
  /** 指标数据 */
  metrics?: MetricsData;
  /** 错误信息 */
  error?: string;
}

// ============================================================
// 评论
// ============================================================

/** 评论数据 */
export interface CommentResult {
  /** 评论 ID */
  id: string;
  /** 作者昵称 */
  author: string;
  /** 作者头像 URL */
  author_avatar?: string;
  /** 评论内容 */
  content: string;
  /** 点赞数 */
  likes: number;
  /** 发布时间 */
  created_at: string;
}

// ============================================================
// 平台适配器接口
// ============================================================

/** 平台适配器 — 每个平台需实现此接口 */
export interface PlatformAdapter {
  /** 平台标识 */
  platform: SocialPlatform;

  /** 获取登录页 URL */
  getLoginUrl(): string;

  /** 获取发布页 URL */
  getPublishUrl(): string;

  /** 检测登录是否成功（通过页面状态判断） */
  detectLoginSuccess(page: Page): Promise<boolean>;

  /** 执行发布流程 */
  executePublish(page: Page, content: PublishContent): Promise<PublishResult>;

  /** 采集数据指标 */
  scrapeMetrics(page: Page): Promise<MetricsResult>;

  /** 采集评论列表 */
  scrapeComments(page: Page, limit: number): Promise<CommentResult[]>;

  /** 回复评论 */
  replyComment(
    page: Page,
    commentId: string,
    replyText: string,
  ): Promise<{ success: boolean; error?: string }>;
}
