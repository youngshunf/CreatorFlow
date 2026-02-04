/**
 * Bun 路径解析器
 * 提供获取内置 Bun 可执行文件路径的统一接口
 *
 * 从 sessions.ts 提取的 Bun 路径解析逻辑，支持：
 * - 打包模式：使用 vendor/bun/ 目录下的 Bun 可执行文件
 * - 开发模式：使用系统 PATH 中的 bun 命令
 * - 跨平台支持：Windows、macOS、Linux
 *
 * Requirements: 1.1, 1.2, 1.3, 1.4, 1.5, 1.6
 */

import { app } from 'electron'
import { join } from 'path'
import { existsSync } from 'fs'
import { execFile } from 'child_process'
import { promisify } from 'util'

const execFileAsync = promisify(execFile)

/**
 * Bun 路径解析器配置选项
 */
export interface BunPathResolverOptions {
  /** 是否为打包模式（默认自动检测 app.isPackaged） */
  isPackaged?: boolean
  /** 应用基础路径（默认自动检测） */
  basePath?: string
  /** 资源路径（Windows 专用，默认 process.resourcesPath） */
  resourcesPath?: string
}

/**
 * Bun 路径解析器接口
 */
export interface BunPathResolver {
  /**
   * 获取 Bun 可执行文件路径
   * @throws Error 如果 Bun 不存在
   */
  getBunPath(): string

  /**
   * 检查 Bun 是否可用
   */
  isBunAvailable(): boolean

  /**
   * 获取 Bun 版本
   */
  getBunVersion(): Promise<string>
}

/**
 * Bun 路径解析错误
 */
export class BunPathError extends Error {
  constructor(
    message: string,
    public readonly code: 'BUN_NOT_FOUND' | 'BUN_EXECUTION_FAILED',
    public readonly path?: string
  ) {
    super(message)
    this.name = 'BunPathError'
  }
}

/**
 * 获取平台特定的 Bun 二进制文件名
 */
function getBunBinaryName(): string {
  return process.platform === 'win32' ? 'bun.exe' : 'bun'
}

/**
 * 获取打包模式下的 Bun 基础路径
 *
 * Windows: 使用 process.resourcesPath (extraResources) 以避免 EBUSY 错误
 * macOS/Linux: 使用应用基础路径
 */
function getPackagedBunBasePath(basePath: string, resourcesPath: string): string {
  return process.platform === 'win32' ? resourcesPath : basePath
}

/**
 * 创建 Bun 路径解析器
 *
 * @param options - 配置选项
 * @returns BunPathResolver 实例
 *
 * @example
 * ```typescript
 * // 使用默认配置（自动检测）
 * const resolver = createBunPathResolver()
 * const bunPath = resolver.getBunPath()
 *
 * // 自定义配置（用于测试）
 * const resolver = createBunPathResolver({
 *   isPackaged: true,
 *   basePath: '/path/to/app',
 *   resourcesPath: '/path/to/resources'
 * })
 * ```
 */
export function createBunPathResolver(options?: BunPathResolverOptions): BunPathResolver {
  // 解析配置选项，使用默认值
  const isPackaged = options?.isPackaged ?? app.isPackaged

  // 获取应用基础路径
  // 在打包模式下，app.getAppPath() 返回 app.asar 路径
  // 我们需要其父目录作为基础路径
  const defaultBasePath = isPackaged
    ? join(app.getAppPath(), '..')
    : app.getAppPath()
  const basePath = options?.basePath ?? defaultBasePath

  // 获取资源路径（Windows 专用）
  const resourcesPath = options?.resourcesPath ?? process.resourcesPath

  // 缓存计算结果
  let cachedBunPath: string | null = null
  let cachedAvailability: boolean | null = null

  /**
   * 计算 Bun 可执行文件路径
   */
  function computeBunPath(): string {
    if (isPackaged) {
      // 打包模式：使用 vendor/bun/ 目录下的 Bun
      const bunBinary = getBunBinaryName()
      const bunBasePath = getPackagedBunBasePath(basePath, resourcesPath)
      return join(bunBasePath, 'vendor', 'bun', bunBinary)
    } else {
      // 开发模式：使用系统 PATH 中的 bun 命令
      return 'bun'
    }
  }

  /**
   * 检查 Bun 是否存在
   */
  function checkBunExists(bunPath: string): boolean {
    if (isPackaged) {
      // 打包模式：检查文件是否存在
      return existsSync(bunPath)
    } else {
      // 开发模式：假设系统 bun 可用
      // 实际可用性通过 getBunVersion() 验证
      return true
    }
  }

  return {
    getBunPath(): string {
      // 使用缓存
      if (cachedBunPath !== null) {
        return cachedBunPath
      }

      const bunPath = computeBunPath()

      // 验证 Bun 是否存在
      if (!checkBunExists(bunPath)) {
        throw new BunPathError(
          `Bundled Bun runtime not found at ${bunPath}. The app package may be corrupted.`,
          'BUN_NOT_FOUND',
          bunPath
        )
      }

      cachedBunPath = bunPath
      return bunPath
    },

    isBunAvailable(): boolean {
      // 使用缓存
      if (cachedAvailability !== null) {
        return cachedAvailability
      }

      try {
        const bunPath = computeBunPath()
        cachedAvailability = checkBunExists(bunPath)
        return cachedAvailability
      } catch {
        cachedAvailability = false
        return false
      }
    },

    async getBunVersion(): Promise<string> {
      const bunPath = this.getBunPath()

      try {
        const { stdout } = await execFileAsync(bunPath, ['--version'], {
          timeout: 5000, // 5 秒超时
        })
        return stdout.trim()
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error)
        throw new BunPathError(
          `Failed to get Bun version: ${message}`,
          'BUN_EXECUTION_FAILED',
          bunPath
        )
      }
    },
  }
}

/**
 * 默认的 Bun 路径解析器实例（单例）
 * 使用自动检测的配置
 */
let defaultResolver: BunPathResolver | null = null

/**
 * 获取默认的 Bun 路径解析器实例
 *
 * @returns 默认的 BunPathResolver 实例
 *
 * @example
 * ```typescript
 * import { getDefaultBunPathResolver } from './bun-path'
 *
 * const resolver = getDefaultBunPathResolver()
 * const bunPath = resolver.getBunPath()
 * ```
 */
export function getDefaultBunPathResolver(): BunPathResolver {
  if (defaultResolver === null) {
    defaultResolver = createBunPathResolver()
  }
  return defaultResolver
}

/**
 * 重置默认解析器（用于测试）
 */
export function resetDefaultBunPathResolver(): void {
  defaultResolver = null
}
