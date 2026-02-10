/**
 * 自媒体创作 APP v2.0 — 数据库模块统一导出
 */

// 连接管理
export {
  type DatabaseDriver,
  type PreparedStatement,
  type DriverFactory,
  CreatorMediaDB,
  createCreatorMediaDB,
  getCachedConnection,
  closeCachedConnection,
  closeAllConnections,
} from './connection.ts';

// Schema
export { SCHEMA_SQL, INITIAL_VERSION_SQL, CURRENT_SCHEMA_VERSION } from './schema.ts';

// 种子数据
export { getSeedSQL } from './seed.ts';

// 迁移
export { getCurrentVersion, initializeSchema, migrate } from './migrations.ts';

// 上下文生成
export { generateProjectContext } from './context-generator.ts';

// 类型（重新导出方便使用）
export type {
  Platform,
  ContentStatus,
  ContentType,
  PipelineMode,
  TopicSource,
  PublishStatus,
  PublishMethod,
  ViralPatternCategory,
  ViralPatternSource,
  TopicCacheStatus,
  TopicCacheSource,
  AuthStatus,
  AuthMethod,
  PostingFrequency,
  Project,
  AccountProfile,
  PlatformAccount,
  Competitor,
  Content,
  PublishRecord,
  ViralPattern,
  TopicCache,
  SchemaVersion,
  CreateProject,
  UpdateProject,
  CreateAccountProfile,
  UpdateAccountProfile,
  CreatePlatformAccount,
  UpdatePlatformAccount,
  CreateCompetitor,
  UpdateCompetitor,
  CreateContent,
  UpdateContent,
  CreatePublishRecord,
  UpdatePublishRecord,
  CreateViralPattern,
  UpdateViralPattern,
  CreateTopicCache,
  ReviewTaskStatus,
  ReviewType,
  ReviewTask,
  CreateReviewTask,
  UpdateReviewTask,
  VersionStage,
  ChangeSource,
  ContentVersion,
  CreateContentVersionInput,
  PublishQueueStatus,
  PublishQueueItem,
  CreatePublishQueueInput,
  UpdatePublishQueueInput,
} from './types.ts';

export { QueuePriority } from './types.ts';
