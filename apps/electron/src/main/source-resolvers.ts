/**
 * Electron 专用的命令/路径解析器
 *
 * 集中管理 Electron 主进程中 MCP source 的命令和工作目录解析逻辑，
 * 供 sessions.ts 和 ipc.ts 共用，避免重复的临时 workaround。
 */

import { app } from 'electron'
import { join } from 'path'
import { getDefaultBunPathResolver } from './bun-path'
import { resolveCommand, resolveCwd } from '@sprouty-ai/shared/sources'
import type { CommandResolver, CwdResolver } from '@sprouty-ai/shared/sources'

/**
 * 获取 Electron 环境的命令解析器。
 * 对 'bun' 命令使用打包的 bundled bun，其余走默认 which 查找。
 */
export function getElectronCommandResolver(): CommandResolver {
  return (command: string) => {
    if (command === 'bun') {
      const bunPath = getDefaultBunPathResolver().getBunPath()
      // 开发模式下 getBunPath() 返回裸 'bun'，需要 resolveCommand 解析为绝对路径
      return bunPath.startsWith('/') ? bunPath : resolveCommand(bunPath)
    }
    return resolveCommand(command)
  }
}

/**
 * 获取 Electron 环境的工作目录解析器。
 * 正确处理 'app:' 前缀：
 * - 打包模式：process.resourcesPath 指向应用资源目录
 * - 开发模式：process.resourcesPath 指向 Electron npm 包内部（不可用），
 *   改用 apps/electron/dist/ 下的资源目录
 */
export function getElectronCwdResolver(): CwdResolver {
  return (cwd: string | undefined, workspaceRootPath: string) => {
    if (cwd?.startsWith('app:')) {
      const relativePath = cwd.slice(4)
      const appRoot = app.isPackaged
        ? (process.resourcesPath ?? join(app.getAppPath(), '..'))
        : join(process.cwd(), 'apps', 'electron', 'dist')
      return join(appRoot, relativePath)
    }
    return resolveCwd(cwd, workspaceRootPath)
  }
}
