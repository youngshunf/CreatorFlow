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
  | 'x'
  | 'toutiao'
  | 'sina'
  | 'sohu';

/** 平台元数据 */
export interface PlatformMeta {
  id: Platform
  label: string
  shortLabel: string
  color: string
  desc: string
}

/** 所有平台配置（单一数据源） */
export const PLATFORM_LIST: PlatformMeta[] = [
  { id: 'xiaohongshu', label: '小红书', shortLabel: '小红书', color: 'text-red-500 border-red-200 dark:border-red-800', desc: '图文/短视频种草' },
  { id: 'douyin', label: '抖音', shortLabel: '抖音', color: 'text-pink-500 border-pink-200 dark:border-pink-800', desc: '短视频/直播' },
  { id: 'bilibili', label: 'B站', shortLabel: 'B站', color: 'text-blue-400 border-blue-200 dark:border-blue-800', desc: '中长视频/专栏' },
  { id: 'wechat', label: '微信公众号', shortLabel: '微信', color: 'text-green-500 border-green-200 dark:border-green-800', desc: '公众号文章' },
  { id: 'zhihu', label: '知乎', shortLabel: '知乎', color: 'text-blue-600 border-blue-200 dark:border-blue-800', desc: '问答/专栏文章' },
  { id: 'weibo', label: '微博', shortLabel: '微博', color: 'text-orange-500 border-orange-200 dark:border-orange-800', desc: '图文/短视频' },
  { id: 'x', label: 'X (Twitter)', shortLabel: 'X', color: 'text-foreground border-border', desc: '推文/长文' },
  { id: 'toutiao', label: '今日头条', shortLabel: '头条', color: 'text-red-600 border-red-200 dark:border-red-800', desc: '图文/短视频/微头条' },
  { id: 'sina', label: '新浪', shortLabel: '新浪', color: 'text-orange-600 border-orange-200 dark:border-orange-800', desc: '新闻/博客' },
  { id: 'sohu', label: '搜狐', shortLabel: '搜狐', color: 'text-yellow-600 border-yellow-200 dark:border-yellow-800', desc: '图文/视频号' },
]

/** 平台 id → 元数据 快速查找 */
export const PLATFORM_MAP: Record<Platform, PlatformMeta> =
  PLATFORM_LIST.reduce((m, p) => { m[p.id] = p; return m }, {} as Record<Platform, PlatformMeta>)

/** 所有平台 id 列表（用于 prompt 等场景） */
export const PLATFORM_IDS: string = PLATFORM_LIST.map(p => p.id).join(' / ')

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
export type AuthMethod = 'cookie' | 'oauth' | 'api_key' | 'browser_profile';

/** 发布频率 */
export type PostingFrequency = 'daily' | '3_per_week' | 'weekly' | 'biweekly' | 'monthly';

export const POSTING_FREQUENCY_LIST: { id: PostingFrequency; label: string }[] = [
  { id: 'daily', label: '每天' },
  { id: '3_per_week', label: '每周3次' },
  { id: 'weekly', label: '每周1次' },
  { id: 'biweekly', label: '每两周1次' },
  { id: 'monthly', label: '每月1次' },
];

/** 采集任务状态 */
export type ReviewTaskStatus = 'pending' | 'executing' | 'completed' | 'failed' | 'cancelled';

/** 采集时间点类型 */
export type ReviewType = '1h' | '6h' | '24h' | '72h' | '168h';

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
  pillar_weights: string | null;   // JSON: Record<string, number> 内容支柱动态权重
  pillar_weights_updated_at: string | null;
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
  profile_path: string | null;     // 浏览器 Profile 目录路径
  fingerprint_id: string | null;   // 关联的指纹配置 ID
  login_check_interval: number;    // 登录检查间隔秒数，默认 3600
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
  review_summary: string | null;    // JSON: 复盘摘要
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
  next_review_at: string | null;     // 下次采集时间
  review_count: number;              // 已完成采集次数
  review_schedule: string | null;    // JSON: 采集计划
  feedback_processed: number;        // SQLite boolean: 0 | 1，反馈是否已回写
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
  historical_score: number | null;   // 同类选题历史成功率
  locked_by_content_id: string | null; // 并行冲突控制，锁定选题
  fetched_at: string;
  expires_at: string | null;
  created_at: string;
}

/** review_tasks — 采集调度任务表 */
export interface ReviewTask {
  id: string;
  publish_record_id: string;
  scheduled_at: string;
  executed_at: string | null;
  status: ReviewTaskStatus;
  review_type: ReviewType;
  result_snapshot: string | null;  // JSON: 采集结果快照
  error_message: string | null;
  retry_count: number;
  created_at: string;
  updated_at: string;
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

export type CreateReviewTask = Omit<ReviewTask, 'created_at' | 'updated_at'>;
export type UpdateReviewTask = Partial<Omit<ReviewTask, 'id' | 'publish_record_id' | 'created_at'>> & { updated_at?: string };

// ============================================================
// 内容版本管理
// ============================================================

/** 版本阶段 */
export type VersionStage = 'script' | 'content' | 'adapted';

/** 变更来源 */
export type ChangeSource = 'auto' | 'user_edit' | 'rollback';

/** content_versions — 内容版本管理表 */
export interface ContentVersion {
  id: string;
  content_id: string;
  version_number: number;
  stage: VersionStage;
  title: string | null;
  content_snapshot: string;          // JSON: 内容快照
  files_snapshot: string | null;     // JSON: 文件列表快照
  change_source: ChangeSource | null;
  change_description: string | null;
  created_by: string;
  created_at: string;
}

export type CreateContentVersionInput = Omit<ContentVersion, 'created_at'>;

// ============================================================
// 发布队列
// ============================================================

/** 发布队列状态 */
export type PublishQueueStatus = 'queued' | 'processing' | 'completed' | 'failed' | 'cancelled';

/** 队列优先级常量 */
export const QueuePriority = {
  /** 热点内容 */
  HOT: 10,
  /** 计划发布 */
  SCHEDULED: 5,
  /** 常青内容 */
  EVERGREEN: 1,
} as const;

/** publish_queue — 发布队列表 */
export interface PublishQueueItem {
  id: string;
  content_id: string;
  platform_account_id: string;
  priority: number;
  status: PublishQueueStatus;
  scheduled_at: string | null;
  started_at: string | null;
  completed_at: string | null;
  error_message: string | null;
  retry_count: number;
  max_retries: number;
  created_at: string;
  updated_at: string;
}

export type CreatePublishQueueInput = Omit<PublishQueueItem, 'created_at' | 'updated_at'>;
export type UpdatePublishQueueInput = Partial<Omit<PublishQueueItem, 'id' | 'content_id' | 'created_at'>> & { updated_at?: string };

// ============================================================
// 草稿
// ============================================================

/** 草稿 — 未进入创作流水线的内容片段 */
export interface Draft {
  id: string;
  project_id: string;
  title: string | null;
  content: string;
  content_type: ContentType | null;
  media: string;                    // JSON: string[]
  tags: string | null;              // JSON: string[]
  target_platforms: string | null;  // JSON: Platform[]
  metadata: string | null;          // JSON: 扩展元数据
  created_at: string;
  updated_at: string;
}

export type CreateDraft = Omit<Draft, 'created_at' | 'updated_at'>;
export type UpdateDraft = Partial<Omit<Draft, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

// ============================================================
// 素材库
// ============================================================

/** 素材文件类型 */
export type MediaFileType = 'image' | 'video';

/** media_files — 素材库 */
export interface MediaFile {
  id: string;
  project_id: string;
  type: MediaFileType;
  path: string;
  filename: string;
  size: number;
  width: number | null;
  height: number | null;
  duration: number | null;          // 视频时长（秒）
  thumbnail: string | null;
  tags: string | null;              // JSON: string[]
  description: string | null;
  created_at: string;
}

export type CreateMediaFile = Omit<MediaFile, 'created_at'>;

// ============================================================
// 热榜快照
// ============================================================

/** 热榜数据来源 */
export type HotTopicFetchSource = 'newsnow' | 'self-hosted' | 'ai-scout';

/** hot_topics — 热榜快照表 */
export interface HotTopic {
  id: string;
  platform_id: string;
  platform_name: string;
  title: string;
  url: string | null;
  rank: number | null;
  heat_score: number | null;
  fetch_source: HotTopicFetchSource;
  fetched_at: string;
  batch_date: string;
  created_at: string;
}

export type CreateHotTopic = Omit<HotTopic, 'created_at'>;

// ============================================================
// 选题推荐
// ============================================================

/** 选题推荐状态: 0=待选 1=已采纳 2=已忽略 */
export type TopicRecommendStatus = 0 | 1 | 2;

/** recommended_topics — 选题推荐表（精简版，详情存 Markdown 文件） */
export interface RecommendedTopic {
  id: string;
  project_id: string;
  title: string;
  reason: string | null;             // 推荐理由
  potential_score: number;
  heat_index: number;
  status: TopicRecommendStatus;
  content_id: string | null;
  md_file_path: string | null;      // Markdown 文件相对路径
  source_uid: string | null;
  batch_date: string | null;
  created_at: string;
  updated_at: string | null;
}

export type CreateRecommendedTopic = Omit<RecommendedTopic, 'created_at' | 'updated_at'>;
export type UpdateRecommendedTopic = Partial<Omit<RecommendedTopic, 'id' | 'project_id' | 'created_at'>> & { updated_at?: string };

// ============================================================
// 视频内容元数据（存储在 contents.metadata JSON 字段）
// ============================================================

/** 视频渲染状态 */
export type VideoRenderStatus = 'not_started' | 'rendering' | 'completed' | 'failed';

/** 视频内容元数据 — 存储在 contents.metadata JSON 字段 */
export interface ContentVideoMetadata {
  /** 关联的 VideoProject.id */
  video_project_id?: string;
  /** 冗余存储项目名称，列表展示用 */
  video_project_name?: string;
  /** 使用的模板 ID */
  video_template_id?: string;
  /** 渲染状态 */
  video_render_status?: VideoRenderStatus;
  /** 最终渲染输出路径 */
  video_output_path?: string;
  /** 视频时长（秒） */
  video_duration?: number;
  /** 视频分辨率 */
  video_resolution?: { width: number; height: number };
}
