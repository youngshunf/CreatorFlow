/**
 * SDK Slash Command 静态数据
 *
 * 纯数据文件，不依赖 Node.js API（fs/path），可安全在渲染进程中使用。
 * 需要 fs 的函数（loadPluginTranslations）保留在 slash-command-translations.ts 中。
 */

export interface SlashCommandTranslation {
  label: string
  description: string
}

/**
 * SDK 内置命令的中文翻译
 */
export const SDK_COMMAND_TRANSLATIONS: Record<string, SlashCommandTranslation> = {
  'compact': { label: '压缩对话', description: '压缩对话历史以释放上下文空间' },
  'commit': { label: '提交代码', description: '将当前更改提交到 Git 仓库' },
  'review-pr': { label: '审查 PR', description: '审查 GitHub Pull Request' },
  'help': { label: '帮助', description: '获取使用帮助' },
  'init': { label: '初始化项目', description: '初始化项目配置' },
}

/**
 * 获取命令的语义化显示信息
 *
 * 优先级：Plugin 翻译 > SDK 翻译 > SDK 原始描述 > 命令名
 */
export function getCommandDisplay(
  cmd: { name: string; description: string; argumentHint: string },
  pluginTranslations?: Record<string, SlashCommandTranslation>
): {
  label: string
  description: string
  hasArgs: boolean
} {
  const pluginTrans = pluginTranslations?.[cmd.name]
  const sdkTrans = SDK_COMMAND_TRANSLATIONS[cmd.name]
  const translation = pluginTrans ?? sdkTrans
  return {
    label: translation?.label ?? cmd.name,
    description: translation?.description ?? cmd.description ?? '',
    hasArgs: !!cmd.argumentHint,
  }
}
