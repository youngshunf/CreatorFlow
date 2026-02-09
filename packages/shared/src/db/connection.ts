/**
 * 自媒体创作 APP v2.0 — 数据库连接管理
 *
 * DatabaseDriver 接口抽象 + CreatorMediaDB 封装 + 连接缓存
 */

// ============================================================
// 驱动接口（由 Electron 主进程注入具体实现）
// ============================================================

/** 预编译语句接口 */
export interface PreparedStatement<T = unknown> {
  run(...params: unknown[]): { changes: number; lastInsertRowid: number | bigint };
  get(...params: unknown[]): T | undefined;
  all(...params: unknown[]): T[];
}

/** 数据库驱动接口 */
export interface DatabaseDriver {
  /** 执行原始 SQL（建表、多语句等） */
  exec(sql: string): void;
  /** 预编译 SQL 语句 */
  prepare<T = unknown>(sql: string): PreparedStatement<T>;
  /** 关闭数据库连接 */
  close(): void;
  /** 执行 PRAGMA 命令 */
  pragma(name: string, value?: unknown): unknown;
  /** 事务执行 */
  transaction<R>(fn: () => R): R;
}

// ============================================================
// CreatorMediaDB 封装
// ============================================================

/** 数据库封装类，持有 driver 实例 */
export class CreatorMediaDB {
  constructor(public readonly driver: DatabaseDriver) {}

  /** 执行原始 SQL */
  exec(sql: string): void {
    this.driver.exec(sql);
  }

  /** 预编译 SQL 语句 */
  prepare<T = unknown>(sql: string): PreparedStatement<T> {
    return this.driver.prepare<T>(sql);
  }

  /** 执行 PRAGMA */
  pragma(name: string, value?: unknown): unknown {
    return this.driver.pragma(name, value);
  }

  /** 事务执行 */
  transaction<R>(fn: () => R): R {
    return this.driver.transaction(fn);
  }

  /** 关闭连接 */
  close(): void {
    this.driver.close();
  }
}

// ============================================================
// 工厂函数 + 连接缓存
// ============================================================

/** 驱动工厂函数类型 */
export type DriverFactory = (dbPath: string) => DatabaseDriver;

/** 连接缓存（按 dbPath 单例） */
const connectionCache = new Map<string, CreatorMediaDB>();

/**
 * 创建或获取 CreatorMediaDB 实例
 * @param dbPath - 数据库文件路径
 * @param driverFactory - 驱动工厂函数（由调用方提供具体实现）
 */
export function createCreatorMediaDB(dbPath: string, driverFactory: DriverFactory): CreatorMediaDB {
  const cached = connectionCache.get(dbPath);
  if (cached) return cached;

  const driver = driverFactory(dbPath);

  // 启用 WAL 模式和外键约束
  driver.pragma('journal_mode', 'WAL');
  driver.pragma('foreign_keys', 'ON');

  const db = new CreatorMediaDB(driver);
  connectionCache.set(dbPath, db);
  return db;
}

/**
 * 获取已缓存的连接
 * @param dbPath - 数据库文件路径
 */
export function getCachedConnection(dbPath: string): CreatorMediaDB | undefined {
  return connectionCache.get(dbPath);
}

/**
 * 关闭并移除缓存的连接
 * @param dbPath - 数据库文件路径
 */
export function closeCachedConnection(dbPath: string): void {
  const cached = connectionCache.get(dbPath);
  if (cached) {
    cached.close();
    connectionCache.delete(dbPath);
  }
}

/**
 * 关闭所有缓存的连接
 */
export function closeAllConnections(): void {
  for (const [path, db] of connectionCache) {
    db.close();
    connectionCache.delete(path);
  }
}
