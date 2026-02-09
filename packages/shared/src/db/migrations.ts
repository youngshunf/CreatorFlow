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
  // 示例：未来 v2 迁移
  // {
  //   version: 2,
  //   description: '添加 xxx 字段',
  //   up: (db) => { db.exec('ALTER TABLE ...'); }
  // },
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
