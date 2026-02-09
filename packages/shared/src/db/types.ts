/**
 * 自媒体创作 APP v2.0 — 数据库类型定义
 *
 * 所有表的 TypeScript 接口 + 枚举常量
 */

// ============================================================
// 枚举常量
// ============================================================

/** 主平台 */
export type Platform =
  | 'xiaohongshu'
  | 'douyin'
  | 'wechat'
  | 'bilibili'
  | 'zhihu'
  | 'weibo'
  | 'x';

/** 内容状态流转 */
export type ContentStatus =
  | 'idea'
  | 'researching'
  | 'scripting'
  | 'creating'
  | 'reviewing'
  | 'scheduled'
  | 'published'
  | 'archived';

/** 内容类型 */
export type ContentType =
  | 'video'
  | 'image-text'
  | 'article'
  | 'short-video'
  | 'live';

/** 流水线模式 */
export type PipelineMode = 'auto' | 'semi-auto' | 'manual';

/** 选题来源 */
export type TopicSource = 'hot_topic' | 'competitor' | 'manual' | 'evergreen';

/** 发布记录状态 */
export type PublishStatus = 'pending' | 'publishing' | 'success' | 'failed' | 'deleted';

/** 发布方式 */
export type PublishMethod = 'api' | 'browser-agent' | 'manual';

/** 爆款模式分类 */
export type ViralPatternCategory =
  | 'hook'
  | 'structure'
  | 'title'
  | 'cta'
  | 'visual'
  | 'rhythm';

/** 爆款模式来源 */
export type ViralPatternSource = 'competitor_analysis' | 'manual' | 'ai_discovered';

/** 热点缓存状态 */
export type TopicCacheStatus = 'new' | 'selected' | 'dismissed' | 'expired';

/** 热点缓存来源 */
export type TopicCacheSource = 'newsnow' | 'trendradar';

/** 平台账号登录状态 */
export type AuthStatus = 'not_logged_in' | 'logged_in' | 'expired' | 'error';

/** 平台账号登录方式 */
export type AuthMethod = 'cookie' | 'oauth' | 'api_key';

/** 发布频率 */
export type PostingFrequency = 'daily' | '3_per_week' | 'weekly' | 'biweekly' | 'monthly';

// ============================================================
// 表接口
// ============================================================

/** schema_version — 版本管理表 */
export interface SchemaVersion {
  version: number;
  applied_at: string;
  description: string | null;
}

/** projects — 账号项目表 */
export interface Project {
  id: string;
  name: string;
  description: string | null;
  platform: Platform;
  platforms: string | null;       // JSON: string[]
  avatar_path: string | null;
  is_active: number;              // SQLite boolean: 0 | 1
  created_at: string;
  updated_at: string;
}

/** account_profiles — 账号画像表 */
export interface AccountProfile {
  id: string;
  project_id: string;
  niche: string;
  sub_niche: string | null;
  persona: string | null;
  target_audience: string | null;
  tone: string | null;
  keywords: string | null;        // JSON: string[]
  bio: string | null;
  content_pillars: string | null; // JSON: string[]
  posting_frequency: string | null;
  best_posting_time: string | null; // JSON: Record<string, string>
  style_references: string | null;  // JSON: string[]
  taboo_topics: string | null;     // JSON: string[]
  created_at: string;
  updated_at: string;
}

/** platform_accounts — 平台账号表 */
export interface PlatformAccount {
  id: string;
  project_id: string;
  platform: Platform;
  platform_uid: string | null;
  nickname: string | null;
  avatar_url: string | null;
  bio: string | null;
  home_url: string | null;
  followers: number;
  following: number;
  total_likes: number;
  total_favorites: number;
  total_comments: number;
  total_posts: number;
  metrics_json: string | null;     // JSON: 平台特有指标
  metrics_updated_at: string | null;
  auth_status: AuthStatus;
  auth_method: string | null;
  auth_data: string | null;        // 加密存储
  auth_expires_at: string | null;
  last_login_at: string | null;
  last_login_check: string | null;
  is_primary: number;              // SQLite boolean: 0 | 1
  notes: string | null;
  created_at: string;
  updated_at: string;
}

/** competitors — 竞品账号表 */
export interface Competitor {
  id: string;
  project_id: string;
  name: string;
  platform: Platform;
  url: string | null;
  follower_count: number | null;
  avg_likes: number | null;
  content_style: string | null;
  strengths: string | null;
  notes: string | null;
  tags: string | null;             // JSON: string[]
  last_analyzed: string | null;
  created_at: string;
  updated_at: string;
}

/** contents — 内容创作表 */
export interface Content {
  id: string;
  project_id: string;
  title: string | null;
  topic: string | null;
  topic_source: TopicSource | null;
  source_topic_id: string | null;
  script_path: string | null;
  status: ContentStatus;
  content_type: ContentType | null;
  target_platforms: string | null;  // JSON: string[]
  pipeline_mode: PipelineMode;
  pipeline_state: string | null;    // JSON: 流水线状态快照
  viral_pattern_id: string | null;
  tags: string | null;              // JSON: string[]
  scheduled_at: string | null;
  files: string | null;             // JSON: string[]
  metadata: string | null;          // JSON: 扩展元数据
  created_at: string;
  updated_at: string;
}

/** publish_records — 发布记录表 */
export interface PublishRecord {
  id: string;
  content_id: string;
  platform_account_id: string | null;
  platform: Platform;
  publish_url: string | null;
  status: PublishStatus;
  method: PublishMethod | null;
  error_message: string | null;
  published_at: string | null;
  views: number;
  likes: number;
  comments: number;
  shares: number;
  favorites: number;
  metrics_json: string | null;      // JSON: 完整指标快照
  metrics_updated_at: string | null;
  created_at: string;
  updated_at: string;
}

/** viral_patterns — 爆款模式库 */
export interface ViralPattern {
  id: string;
  project_id: string | null;       // NULL = 全局模式
  platform: Platform | null;       // NULL = 全平台
  category: ViralPatternCategory;
  name: string;
  description: string | null;
  template: string | null;
  examples: string | null;          // JSON: Array<{title, url, metrics}>
  source: ViralPatternSource | null;
  usage_count: number;
  success_rate: number | null;
  tags: string | null;              // JSON: string[]
  created_at: string;
  updated_at: string;
}

/** topic_cache — 热点缓存表 */
export interface TopicCache {
  id: string;
  source: TopicCacheSource;
  source_id: string;
  title: string;
  url: string | null;
  heat_score: number | null;
  relevance_score: number | null;
  category: string | null;
  keywords: string | null;          // JSON: string[]
  status: TopicCacheStatus;
  fetched_at: string;
  expires_at: string | null;
  created_at: string;
}

// ============================================================
// 创建/更新输入类型（省略自动生成字段）
// ============================================================

export type CreateProject = Omit<Project, 'created_at' | 'updated_at'>;
export type UpdateProject = Partial<Omit<Project, 'id' | 'created_at'>> & { updated_at?: string };

export type CreateAccountProfile = Omit<AccountProfile, 'created_at' | 'updated_at'>;
export type UpdateAccountProfile = Partial<Omit<AccountProfile, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

export type CreatePlatformAccount = Omit<PlatformAccount, 'created_at' | 'updated_at'>;
export type UpdatePlatformAccount = Partial<Omit<PlatformAccount, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

export type CreateCompetitor = Omit<Competitor, 'created_at' | 'updated_at'>;
export type UpdateCompetitor = Partial<Omit<Competitor, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

export type CreateContent = Omit<Content, 'created_at' | 'updated_at'>;
export type UpdateContent = Partial<Omit<Content, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

export type CreatePublishRecord = Omit<PublishRecord, 'created_at' | 'updated_at'>;
export type UpdatePublishRecord = Partial<Omit<PublishRecord, 'id' | 'content_id' | 'created_at'>> & { updated_at?: string };

export type CreateViralPattern = Omit<ViralPattern, 'created_at' | 'updated_at'>;
export type UpdateViralPattern = Partial<Omit<ViralPattern, 'id' | 'created_at'>> & { updated_at?: string };

export type CreateTopicCache = Omit<TopicCache, 'created_at'>;
