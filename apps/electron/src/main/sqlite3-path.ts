/**
 * sqlite3 跨平台路径管理
 *
 * 内置各平台 sqlite3 二进制，启动时加入 PATH
 * 解决 Windows 用户没有 sqlite3 命令的问题
 */

import path from 'path';

/**
 * 获取 sqlite3 二进制所在目录
 */
export function getSqlite3BinDir(): string {
  const platform = process.platform;
  const arch = process.arch;
  const ext = platform === 'win32' ? '.exe' : '';
  const binName = `sqlite3-${platform}-${arch}${ext}`;

  // 判断是否打包环境
  const isPackaged = typeof process.resourcesPath === 'string'
    && !process.resourcesPath.includes('node_modules');

  if (isPackaged) {
    return path.join(process.resourcesPath!, 'bin');
  }
  // 开发模式: 从项目根目录的 resources/bin/ 读取
  return path.join(__dirname, '../../../../resources/bin');
}

/**
 * 将 sqlite3 二进制目录加入 PATH
 * 在应用启动时调用
 */
export function addSqlite3ToPath(): void {
  const binDir = getSqlite3BinDir();
  const sep = path.delimiter;
  process.env.PATH = `${binDir}${sep}${process.env.PATH}`;
}
