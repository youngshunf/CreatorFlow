/**
 * Electron 主进程 — better-sqlite3 驱动实现
 *
 * 实现 DatabaseDriver 接口，包装 better-sqlite3
 * 注意：better-sqlite3 使用懒加载，避免启动时加载原生模块导致崩溃
 */

import { existsSync, mkdirSync } from 'fs';
import { dirname } from 'path';
import { getWorkspaceDbPath } from '@sprouty-ai/shared/workspaces';
import {
  type DatabaseDriver,
  type PreparedStatement,
  createCreatorMediaDB,
  getCachedConnection,
  closeCachedConnection,
  type CreatorMediaDB,
  migrate,
} from '@sprouty-ai/shared/db';

/** 懒加载 better-sqlite3（缓存模块引用） */
let _Database: any = null;
function loadBetterSqlite3() {
  if (!_Database) {
    try {
      // eslint-disable-next-line @typescript-eslint/no-require-imports
      _Database = require('better-sqlite3');
    } catch (err) {
      throw new Error(
        `better-sqlite3 原生模块加载失败，请运行 npx electron-rebuild -f -w better-sqlite3 重新编译。原始错误: ${err}`
      );
    }
  }
  return _Database;
}

/** better-sqlite3 驱动实现 */
class BetterSqlite3Driver implements DatabaseDriver {
  private db: any;

  constructor(dbPath: string) {
    const Database = loadBetterSqlite3();
    this.db = new Database(dbPath);
  }

  exec(sql: string): void {
    this.db.exec(sql);
  }

  prepare<T = unknown>(sql: string): PreparedStatement<T> {
    const stmt = this.db.prepare(sql);
    return {
      run: (...params: unknown[]) => stmt.run(...params) as { changes: number; lastInsertRowid: number | bigint },
      get: (...params: unknown[]) => stmt.get(...params) as T | undefined,
      all: (...params: unknown[]) => stmt.all(...params) as T[],
    };
  }

  close(): void {
    this.db.close();
  }

  pragma(name: string, value?: unknown): unknown {
    if (value !== undefined) {
      return this.db.pragma(`${name} = ${value}`);
    }
    return this.db.pragma(name);
  }

  transaction<R>(fn: () => R): R {
    const txn = this.db.transaction(fn);
    return txn();
  }
}

/** better-sqlite3 驱动工厂 */
function betterSqlite3Factory(dbPath: string): DatabaseDriver {
  return new BetterSqlite3Driver(dbPath);
}

/**
 * 初始化 creator-media 数据库
 * 创建目录 + 初始化 schema + 种子数据
 */
export function initCreatorMediaDB(workspaceRoot: string): CreatorMediaDB {
  const dbPath = getWorkspaceDbPath(workspaceRoot);

  // 确保 db 目录存在
  const dbDir = dirname(dbPath);
  if (!existsSync(dbDir)) {
    mkdirSync(dbDir, { recursive: true });
  }

  // 创建连接（带 WAL + foreign_keys）
  const db = createCreatorMediaDB(dbPath, betterSqlite3Factory);

  // 执行迁移（首次会创建所有表 + 种子数据）
  migrate(db);

  return db;
}

/**
 * 获取 creator-media 数据库连接（单例缓存）
 * 如果连接不存在，自动初始化
 */
export function getCreatorMediaDB(workspaceRoot: string): CreatorMediaDB {
  const dbPath = getWorkspaceDbPath(workspaceRoot);

  // 尝试获取缓存连接
  const cached = getCachedConnection(dbPath);
  if (cached) return cached;

  // 不存在则初始化
  return initCreatorMediaDB(workspaceRoot);
}

/**
 * 关闭 creator-media 数据库连接
 */
export function closeCreatorMediaDB(workspaceRoot: string): void {
  const dbPath = getWorkspaceDbPath(workspaceRoot);
  closeCachedConnection(dbPath);
}
