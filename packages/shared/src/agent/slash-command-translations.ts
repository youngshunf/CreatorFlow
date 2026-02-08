/**
 * SDK Slash Command 翻译映射
 *
 * 两层翻译：
 * 1. SDK 内置命令 - 硬编码在 slash-command-data.ts 中
 * 2. Plugin 命令 - 从 plugin 目录的 .claude-plugin/translations.json 加载（需要 fs）
 */

import { existsSync, readFileSync, readdirSync, statSync } from 'fs'
import { join, basename } from 'path'
import { SDK_COMMAND_TRANSLATIONS, type SlashCommandTranslation } from './slash-command-data.ts'

// Re-export pure data (no fs dependency) for convenience
export { SDK_COMMAND_TRANSLATIONS, getCommandDisplay, type SlashCommandTranslation } from './slash-command-data.ts'

/**
 * 从 plugin 目录加载翻译文件
 *
 * 翻译文件格式: {plugin-dir}/.claude-plugin/translations.json
 * ```json
 * {
 *   "locale": "zh-CN",
 *   "commands": {
 *     "command-name": { "label": "中文名", "description": "描述" }
 *   }
 * }
 * ```
 */
export function loadPluginTranslations(
  pluginPaths: Array<{ name: string; path: string }>
): Record<string, SlashCommandTranslation> {
  const translations: Record<string, SlashCommandTranslation> = {}
  for (const plugin of pluginPaths) {
    const translationFile = join(plugin.path, '.claude-plugin', 'translations.json')
    if (!existsSync(translationFile)) continue
    try {
      const data = JSON.parse(readFileSync(translationFile, 'utf-8'))
      if (data.commands && typeof data.commands === 'object') {
        for (const [cmdName, trans] of Object.entries(data.commands)) {
          const t = trans as { label?: string; description?: string }
          if (t.label) {
            // Store with plugin-prefixed key: pluginName:commandName
            translations[`${plugin.name}:${cmdName}`] = {
              label: t.label,
              description: t.description ?? '',
            }
          }
        }
      }
    } catch {
      // 静默忽略解析错误
    }
  }
  return translations
}

/**
 * 扫描单个插件目录，发现已安装插件及其命令
 *
 * 目录结构:
 * {pluginBaseDir}/
 *   plugin.json          (plugin manifest)
 *   {pluginName}/        (installed plugin)
 *     commands/
 *       {command}.md
 *     .claude-plugin/
 *       plugin.json
 *       translations.json
 */
function scanPluginDirectory(
  pluginBaseDir: string,
  commands: Array<{ name: string; description: string; argumentHint: string }>,
  translations: Record<string, SlashCommandTranslation>,
): void {
  if (!existsSync(pluginBaseDir)) return

  try {
    const entries = readdirSync(pluginBaseDir)
    for (const entry of entries) {
      const pluginDir = join(pluginBaseDir, entry)
      // Skip files (like plugin.json, .DS_Store)
      try {
        if (!statSync(pluginDir).isDirectory()) continue
      } catch {
        continue
      }

      const pluginName = entry
      const commandsDir = join(pluginDir, 'commands')

      // Load plugin translations
      const translationFile = join(pluginDir, '.claude-plugin', 'translations.json')
      if (existsSync(translationFile)) {
        try {
          const data = JSON.parse(readFileSync(translationFile, 'utf-8'))
          if (data.commands && typeof data.commands === 'object') {
            for (const [cmdName, trans] of Object.entries(data.commands)) {
              const t = trans as { label?: string; description?: string }
              if (t.label) {
                translations[`${pluginName}:${cmdName}`] = {
                  label: t.label,
                  description: t.description ?? '',
                }
              }
            }
          }
        } catch {
          // 静默忽略解析错误
        }
      }

      // Scan commands directory
      if (!existsSync(commandsDir)) continue
      try {
        const cmdFiles = readdirSync(commandsDir)
        for (const file of cmdFiles) {
          if (!file.endsWith('.md')) continue
          const cmdName = basename(file, '.md')
          const fullName = `${pluginName}:${cmdName}`

          // Skip if already added (e.g. from global plugin)
          if (commands.some(c => c.name === fullName)) continue

          // Read first line of .md file for description (fallback)
          let description = ''
          try {
            const content = readFileSync(join(commandsDir, file), 'utf-8')
            const firstLine = content.split('\n').find(l => l.trim() && !l.startsWith('#') && !l.startsWith('---'))
            description = firstLine?.trim() ?? ''
          } catch {
            // ignore
          }

          // Use translation description if available
          const trans = translations[fullName]
          commands.push({
            name: fullName,
            description: trans?.description ?? description,
            argumentHint: '',
          })
        }
      } catch {
        // 静默忽略目录读取错误
      }
    }
  } catch {
    // 静默忽略目录读取错误
  }
}

/**
 * 扫描工作区和全局插件目录，发现所有已安装插件及其命令
 *
 * 返回: SDK 内置命令 + 全局插件命令 + 工作区插件命令的完整列表，以及合并后的翻译
 */
export function scanWorkspaceCommands(
  workspaceRootPath: string,
  globalPluginDataPath?: string,
): {
  commands: Array<{ name: string; description: string; argumentHint: string }>
  translations: Record<string, SlashCommandTranslation>
} {
  // Start with SDK built-in commands
  const commands: Array<{ name: string; description: string; argumentHint: string }> = Object.entries(SDK_COMMAND_TRANSLATIONS).map(([name, trans]) => ({
    name,
    description: trans.description,
    argumentHint: '',
  }))

  const translations: Record<string, SlashCommandTranslation> = {}

  // Scan global plugins (if provided)
  if (globalPluginDataPath) {
    scanPluginDirectory(join(globalPluginDataPath, '.claude-plugin'), commands, translations)
  }

  // Scan workspace plugins
  scanPluginDirectory(join(workspaceRootPath, '.sprouty-ai', '.claude-plugin'), commands, translations)

  return { commands, translations }
}
