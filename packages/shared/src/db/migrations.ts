/**
 * 自媒体创作 APP v2.0 — 数据库迁移管理
 *
 * 版本检测 + Schema 初始化 + 增量迁移
 */

import type { CreatorMediaDB } from './connection.ts';
import { SCHEMA_SQL, INITIAL_VERSION_SQL, CURRENT_SCHEMA_VERSION } from './schema.ts';
import { getSeedSQL } from './seed.ts';

/** 迁移定义 */
interface Migration {
  version: number;
  description: string;
  up: (db: CreatorMediaDB) => void;
}

/** 迁移注册表（后续版本在此追加） */
const MIGRATIONS: Migration[] = [
  {
    version: 2,
    description: '添加反馈闭环字段：支柱权重、历史评分、选题锁定、复盘摘要、采集调度',
    up: (db) => {
      // account_profiles: 内容支柱动态权重
      db.exec('ALTER TABLE account_profiles ADD COLUMN pillar_weights TEXT');
      db.exec('ALTER TABLE account_profiles ADD COLUMN pillar_weights_updated_at DATETIME');

      // topic_cache: 历史成功率 + 并行冲突控制
      db.exec('ALTER TABLE topic_cache ADD COLUMN historical_score REAL');
      db.exec('ALTER TABLE topic_cache ADD COLUMN locked_by_content_id TEXT');

      // contents: 复盘摘要
      db.exec('ALTER TABLE contents ADD COLUMN review_summary TEXT');

      // publish_records: 采集调度 + 反馈回写
      db.exec('ALTER TABLE publish_records ADD COLUMN next_review_at DATETIME');
      db.exec('ALTER TABLE publish_records ADD COLUMN review_count INTEGER DEFAULT 0');
      db.exec('ALTER TABLE publish_records ADD COLUMN review_schedule TEXT');
      db.exec('ALTER TABLE publish_records ADD COLUMN feedback_processed BOOLEAN DEFAULT 0');
    },
  },
  {
    version: 3,
    description: '创建采集调度任务表 review_tasks',
    up: (db) => {
      db.exec(`
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
        )
      `);
      db.exec('CREATE INDEX IF NOT EXISTS idx_review_scheduled ON review_tasks(scheduled_at) WHERE status = \'pending\'');
      db.exec('CREATE INDEX IF NOT EXISTS idx_review_publish ON review_tasks(publish_record_id)');
    },
  },
  {
    version: 4,
    description: '创建内容版本管理表 content_versions',
    up: (db) => {
      db.exec(`
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
        )
      `);
      db.exec('CREATE INDEX IF NOT EXISTS idx_version_content ON content_versions(content_id)');
      db.exec('CREATE UNIQUE INDEX IF NOT EXISTS idx_version_number ON content_versions(content_id, version_number)');
    },
  },
  {
    version: 5,
    description: '平台账号添加浏览器 Profile 字段 + 创建发布队列表 publish_queue',
    up: (db) => {
      // platform_accounts: 浏览器 Profile 相关字段
      db.exec('ALTER TABLE platform_accounts ADD COLUMN profile_path TEXT');
      db.exec('ALTER TABLE platform_accounts ADD COLUMN fingerprint_id TEXT');
      db.exec('ALTER TABLE platform_accounts ADD COLUMN login_check_interval INTEGER DEFAULT 3600');

      // 发布队列表
      db.exec(`
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
        )
      `);
      db.exec('CREATE INDEX IF NOT EXISTS idx_queue_status ON publish_queue(status, priority DESC, scheduled_at)');
      db.exec('CREATE INDEX IF NOT EXISTS idx_queue_platform ON publish_queue(platform_account_id, status)');
    },
  },
];

/**
 * 获取当前数据库 Schema 版本
 * 如果 schema_version 表不存在，返回 0
 */
export function getCurrentVersion(db: CreatorMediaDB): number {
  try {
    const row = db.prepare<{ version: number }>(
      'SELECT MAX(version) as version FROM schema_version'
    ).get();
    return row?.version ?? 0;
  } catch {
    // 表不存在
    return 0;
  }
}

/**
 * 初始化数据库 Schema（首次创建）
 * 创建所有表 + 插入版本记录 + 种子数据
 */
export function initializeSchema(db: CreatorMediaDB): void {
  db.exec(SCHEMA_SQL);
  db.exec(INITIAL_VERSION_SQL);
  db.exec(getSeedSQL());
}

/**
 * 执行增量迁移
 * 从当前版本迁移到最新版本
 */
export function migrate(db: CreatorMediaDB): { from: number; to: number; applied: number } {
  const currentVersion = getCurrentVersion(db);

  if (currentVersion === 0) {
    // 全新数据库，执行完整初始化
    initializeSchema(db);
    return { from: 0, to: CURRENT_SCHEMA_VERSION, applied: 1 };
  }

  // 筛选需要执行的迁移
  const pendingMigrations = MIGRATIONS.filter(m => m.version > currentVersion)
    .sort((a, b) => a.version - b.version);

  if (pendingMigrations.length === 0) {
    return { from: currentVersion, to: currentVersion, applied: 0 };
  }

  // 在事务中执行所有迁移
  db.transaction(() => {
    for (const migration of pendingMigrations) {
      migration.up(db);
      db.prepare(
        'INSERT INTO schema_version (version, description) VALUES (?, ?)'
      ).run(migration.version, migration.description);
    }
  });

  const newVersion = pendingMigrations[pendingMigrations.length - 1].version;
  return { from: currentVersion, to: newVersion, applied: pendingMigrations.length };
}
