/**
 * 自媒体创作 APP v2.0 — SQL 建表语句
 *
 * 8 张业务表 + 1 张版本表 + 索引
 */

/** 当前 Schema 版本 */
export const CURRENT_SCHEMA_VERSION = 9;

/** 完整建表 SQL */
export const SCHEMA_SQL = `
-- 版本管理表
CREATE TABLE IF NOT EXISTS schema_version (
  version     INTEGER PRIMARY KEY,
  applied_at  DATETIME DEFAULT CURRENT_TIMESTAMP,
  description TEXT
);

-- 账号项目表
CREATE TABLE IF NOT EXISTS projects (
  id            TEXT PRIMARY KEY,
  name          TEXT NOT NULL,
  description   TEXT,
  platform      TEXT NOT NULL,
  platforms     TEXT,
  avatar_path   TEXT,
  is_active     BOOLEAN DEFAULT 0,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_active_project ON projects(is_active) WHERE is_active = 1;

-- 账号画像表
CREATE TABLE IF NOT EXISTS account_profiles (
  id                TEXT PRIMARY KEY,
  project_id        TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  niche             TEXT NOT NULL,
  sub_niche         TEXT,
  persona           TEXT,
  target_audience   TEXT,
  tone              TEXT,
  keywords          TEXT,
  bio               TEXT,
  content_pillars   TEXT,
  posting_frequency TEXT,
  best_posting_time TEXT,
  style_references  TEXT,
  taboo_topics      TEXT,
  pillar_weights    TEXT,
  pillar_weights_updated_at DATETIME,
  created_at        DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at        DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_profile_project ON account_profiles(project_id);

-- 平台账号表
CREATE TABLE IF NOT EXISTS platform_accounts (
  id              TEXT PRIMARY KEY,
  project_id      TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  platform        TEXT NOT NULL,
  platform_uid    TEXT,
  nickname        TEXT,
  avatar_url      TEXT,
  bio             TEXT,
  home_url        TEXT,
  followers       INTEGER DEFAULT 0,
  following       INTEGER DEFAULT 0,
  total_likes     INTEGER DEFAULT 0,
  total_favorites INTEGER DEFAULT 0,
  total_comments  INTEGER DEFAULT 0,
  total_posts     INTEGER DEFAULT 0,
  metrics_json    TEXT,
  metrics_updated_at DATETIME,
  auth_status     TEXT DEFAULT 'not_logged_in',
  auth_method     TEXT,
  auth_data       TEXT,
  auth_expires_at DATETIME,
  last_login_at   DATETIME,
  last_login_check DATETIME,
  is_primary      BOOLEAN DEFAULT 0,
  notes           TEXT,
  profile_path    TEXT,
  fingerprint_id  TEXT,
  login_check_interval INTEGER DEFAULT 3600,
  created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_account_project ON platform_accounts(project_id);
CREATE INDEX IF NOT EXISTS idx_account_platform ON platform_accounts(platform);
CREATE UNIQUE INDEX IF NOT EXISTS idx_account_primary ON platform_accounts(project_id, platform, is_primary)
  WHERE is_primary = 1;

-- 竞品账号表
CREATE TABLE IF NOT EXISTS competitors (
  id            TEXT PRIMARY KEY,
  project_id    TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  name          TEXT NOT NULL,
  platform      TEXT NOT NULL,
  url           TEXT,
  follower_count INTEGER,
  avg_likes     INTEGER,
  content_style TEXT,
  strengths     TEXT,
  notes         TEXT,
  tags          TEXT,
  last_analyzed DATETIME,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_competitor_project ON competitors(project_id);

-- 内容创作表
CREATE TABLE IF NOT EXISTS contents (
  id              TEXT PRIMARY KEY,
  project_id      TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  title           TEXT,
  topic           TEXT,
  topic_source    TEXT,
  source_topic_id TEXT,
  script_path     TEXT,
  status          TEXT DEFAULT 'idea',
  content_type    TEXT,
  target_platforms TEXT,
  pipeline_mode   TEXT DEFAULT 'semi-auto',
  pipeline_state  TEXT,
  viral_pattern_id TEXT,
  tags            TEXT,
  scheduled_at    DATETIME,
  files           TEXT,
  metadata        TEXT,
  review_summary  TEXT,
  created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_content_project ON contents(project_id);
CREATE INDEX IF NOT EXISTS idx_content_status ON contents(status);
CREATE INDEX IF NOT EXISTS idx_content_scheduled ON contents(scheduled_at) WHERE scheduled_at IS NOT NULL;

-- 发布记录表
CREATE TABLE IF NOT EXISTS publish_records (
  id            TEXT PRIMARY KEY,
  content_id    TEXT NOT NULL REFERENCES contents(id) ON DELETE CASCADE,
  platform_account_id TEXT REFERENCES platform_accounts(id),
  platform      TEXT NOT NULL,
  publish_url   TEXT,
  status        TEXT DEFAULT 'pending',
  method        TEXT,
  error_message TEXT,
  published_at  DATETIME,
  views         INTEGER DEFAULT 0,
  likes         INTEGER DEFAULT 0,
  comments      INTEGER DEFAULT 0,
  shares        INTEGER DEFAULT 0,
  favorites     INTEGER DEFAULT 0,
  metrics_json  TEXT,
  metrics_updated_at DATETIME,
  next_review_at  DATETIME,
  review_count    INTEGER DEFAULT 0,
  review_schedule TEXT,
  feedback_processed BOOLEAN DEFAULT 0,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_publish_content ON publish_records(content_id);
CREATE INDEX IF NOT EXISTS idx_publish_platform ON publish_records(platform);
CREATE INDEX IF NOT EXISTS idx_publish_status ON publish_records(status);

-- 爆款模式库
CREATE TABLE IF NOT EXISTS viral_patterns (
  id            TEXT PRIMARY KEY,
  project_id    TEXT REFERENCES projects(id) ON DELETE SET NULL,
  platform      TEXT,
  category      TEXT NOT NULL,
  name          TEXT NOT NULL,
  description   TEXT,
  template      TEXT,
  examples      TEXT,
  source        TEXT,
  usage_count   INTEGER DEFAULT 0,
  success_rate  REAL,
  tags          TEXT,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_pattern_project ON viral_patterns(project_id);
CREATE INDEX IF NOT EXISTS idx_pattern_category ON viral_patterns(category);
CREATE INDEX IF NOT EXISTS idx_pattern_platform ON viral_patterns(platform);

-- 热点缓存表
CREATE TABLE IF NOT EXISTS topic_cache (
  id            TEXT PRIMARY KEY,
  source        TEXT NOT NULL,
  source_id     TEXT NOT NULL,
  title         TEXT NOT NULL,
  url           TEXT,
  heat_score    REAL,
  relevance_score REAL,
  category      TEXT,
  keywords      TEXT,
  status        TEXT DEFAULT 'new',
  historical_score REAL,
  locked_by_content_id TEXT,
  fetched_at    DATETIME NOT NULL,
  expires_at    DATETIME,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_topic_source ON topic_cache(source, source_id);
CREATE INDEX IF NOT EXISTS idx_topic_status ON topic_cache(status);
CREATE INDEX IF NOT EXISTS idx_topic_relevance ON topic_cache(relevance_score DESC);

-- 采集调度任务表
CREATE TABLE IF NOT EXISTS review_tasks (
  id                TEXT PRIMARY KEY,
  publish_record_id TEXT NOT NULL REFERENCES publish_records(id) ON DELETE CASCADE,
  scheduled_at      DATETIME NOT NULL,
  executed_at       DATETIME,
  status            TEXT DEFAULT 'pending',
  review_type       TEXT NOT NULL,
  result_snapshot   TEXT,
  error_message     TEXT,
  retry_count       INTEGER DEFAULT 0,
  created_at        DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at        DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_review_scheduled ON review_tasks(scheduled_at) WHERE status = 'pending';
CREATE INDEX IF NOT EXISTS idx_review_publish ON review_tasks(publish_record_id);

-- 内容版本管理表
CREATE TABLE IF NOT EXISTS content_versions (
  id                 TEXT PRIMARY KEY,
  content_id         TEXT NOT NULL REFERENCES contents(id) ON DELETE CASCADE,
  version_number     INTEGER NOT NULL,
  stage              TEXT NOT NULL,
  title              TEXT,
  content_snapshot   TEXT NOT NULL,
  files_snapshot     TEXT,
  change_source      TEXT,
  change_description TEXT,
  created_by         TEXT DEFAULT 'system',
  created_at         DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_version_content ON content_versions(content_id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_version_number ON content_versions(content_id, version_number);

-- 发布队列表
CREATE TABLE IF NOT EXISTS publish_queue (
  id                  TEXT PRIMARY KEY,
  content_id          TEXT NOT NULL REFERENCES contents(id) ON DELETE CASCADE,
  platform_account_id TEXT NOT NULL REFERENCES platform_accounts(id),
  priority            INTEGER DEFAULT 0,
  status              TEXT DEFAULT 'queued',
  scheduled_at        DATETIME,
  started_at          DATETIME,
  completed_at        DATETIME,
  error_message       TEXT,
  retry_count         INTEGER DEFAULT 0,
  max_retries         INTEGER DEFAULT 3,
  created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at          DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_queue_status ON publish_queue(status, priority DESC, scheduled_at);
CREATE INDEX IF NOT EXISTS idx_queue_platform ON publish_queue(platform_account_id, status);

-- 草稿表
CREATE TABLE IF NOT EXISTS drafts (
  id              TEXT PRIMARY KEY,
  project_id      TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  title           TEXT,
  content         TEXT NOT NULL,
  content_type    TEXT,
  media           TEXT NOT NULL DEFAULT '[]',
  tags            TEXT,
  target_platforms TEXT,
  metadata        TEXT,
  created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_draft_project ON drafts(project_id);

-- 素材库表
CREATE TABLE IF NOT EXISTS media_files (
  id              TEXT PRIMARY KEY,
  project_id      TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  type            TEXT NOT NULL,
  path            TEXT NOT NULL,
  filename        TEXT NOT NULL,
  size            INTEGER NOT NULL,
  width           INTEGER,
  height          INTEGER,
  duration        INTEGER,
  thumbnail       TEXT,
  tags            TEXT,
  description     TEXT,
  created_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_media_project ON media_files(project_id);
CREATE INDEX IF NOT EXISTS idx_media_type ON media_files(type);

-- 热榜快照表
CREATE TABLE IF NOT EXISTS hot_topics (
  id            TEXT PRIMARY KEY,
  platform_id   TEXT NOT NULL,
  platform_name TEXT NOT NULL,
  title         TEXT NOT NULL,
  url           TEXT,
  rank          INTEGER,
  heat_score    REAL,
  fetch_source  TEXT NOT NULL,
  fetched_at    TEXT NOT NULL,
  batch_date    TEXT NOT NULL,
  created_at    TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_hot_topics_platform ON hot_topics(platform_id);
CREATE INDEX IF NOT EXISTS idx_hot_topics_batch ON hot_topics(batch_date);
CREATE UNIQUE INDEX IF NOT EXISTS idx_hot_topics_dedup ON hot_topics(platform_id, url, batch_date);

-- 选题推荐表
CREATE TABLE IF NOT EXISTS recommended_topics (
  id                TEXT PRIMARY KEY,
  project_id        TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  title             TEXT NOT NULL,
  industry_id       INTEGER,
  potential_score   REAL NOT NULL,
  heat_index        REAL NOT NULL,
  reason            TEXT,
  keywords          TEXT,
  platform_heat     TEXT,
  heat_sources      TEXT,
  trend             TEXT,
  industry_tags     TEXT,
  target_audience   TEXT,
  creative_angles   TEXT,
  content_outline   TEXT,
  format_suggestions TEXT,
  material_clues    TEXT,
  risk_notes        TEXT,
  source_info       TEXT,
  md_file_path      TEXT,
  batch_date        TEXT,
  source_uid        TEXT,
  status            INTEGER NOT NULL DEFAULT 0,
  content_id        TEXT,
  created_at        TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at        TEXT
);

CREATE INDEX IF NOT EXISTS idx_rec_topics_project ON recommended_topics(project_id);
CREATE INDEX IF NOT EXISTS idx_rec_topics_status ON recommended_topics(status);
CREATE INDEX IF NOT EXISTS idx_rec_topics_batch ON recommended_topics(batch_date);
CREATE UNIQUE INDEX IF NOT EXISTS idx_rec_topics_dedup ON recommended_topics(source_uid);
`;

/** 初始版本记录 SQL */
export const INITIAL_VERSION_SQL = `
INSERT OR IGNORE INTO schema_version (version, description) VALUES (7, '新增热榜快照表与选题推荐表');
INSERT OR IGNORE INTO schema_version (version, description) VALUES (8, '移除定时任务表（迁移至 hooks.json）');
INSERT OR IGNORE INTO schema_version (version, description) VALUES (9, '选题推荐表新增 md_file_path 字段');
`;
