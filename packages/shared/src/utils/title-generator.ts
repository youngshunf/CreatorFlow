/**
 * Session title generator utility.
 * Uses Claude Agent SDK query() for all auth types (API Key, Claude OAuth).
 */

import { query } from '@anthropic-ai/claude-agent-sdk';
import { getDefaultOptions } from '../agent/options.ts';
import { SUMMARIZATION_MODEL } from '../config/models.ts';
import { resolveModelId } from '../config/storage.ts';
import { loadPreferences } from '../config/preferences.ts';

/**
 * Check if the language setting indicates Chinese
 */
function isChinese(language: string): boolean {
  const lower = language.toLowerCase();
  return lower.startsWith('zh') || lower === '中文' || lower.includes('chinese');
}

/**
 * Get the prompt for title generation based on locale
 */
function getTitlePrompt(userSnippet: string, locale: string): string {
  if (isChinese(locale)) {
    return [
      '用户想做什么？请只回复一个简短的任务描述（2-5个词）。',
      '以动词开头。只使用纯文本，不要使用 markdown 格式。',
      '示例："修复认证问题"、"添加深色模式"、"重构API层"、"解释代码结构"',
      '',
      '用户: ' + userSnippet,
      '',
      '任务:',
    ].join('\n');
  }
  
  return [
    'What is the user trying to do? Reply with ONLY a short task description (2-5 words).',
    'Start with a verb. Use plain text only - no markdown.',
    'Examples: "Fix authentication bug", "Add dark mode", "Refactor API layer", "Explain codebase structure"',
    '',
    'User: ' + userSnippet,
    '',
    'Task:',
  ].join('\n');
}

/**
 * Get the prompt for title regeneration based on locale
 */
function getRegenerateTitlePrompt(userContext: string, assistantSnippet: string, locale: string): string {
  if (isChinese(locale)) {
    return [
      '根据这些最近的消息，当前对话的焦点是什么？',
      '请只回复一个简短的任务描述（2-5个词）。',
      '以动词开头。只使用纯文本，不要使用 markdown 格式。',
      '示例："修复认证问题"、"添加深色模式"、"重构API层"、"解释代码结构"',
      '',
      '最近的用户消息:',
      userContext,
      '',
      '最新的助手回复:',
      assistantSnippet,
      '',
      '当前焦点:',
    ].join('\n');
  }
  
  return [
    'Based on these recent messages, what is the current focus of this conversation?',
    'Reply with ONLY a short task description (2-5 words).',
    'Start with a verb. Use plain text only - no markdown.',
    'Examples: "Fix authentication bug", "Add dark mode", "Refactor API layer", "Explain codebase structure"',
    '',
    'Recent user messages:',
    userContext,
    '',
    'Latest assistant response:',
    assistantSnippet,
    '',
    'Current focus:',
  ].join('\n');
}

/**
 * Get the current language from preferences or default to '中文'
 */
function getCurrentLanguage(): string {
  try {
    const prefs = loadPreferences();
    return prefs.language || '中文';
  } catch {
    return '中文';
  }
}

/**
 * Generate a task-focused title (2-5 words) from the user's first message.
 * Extracts what the user is trying to accomplish, framing conversations as tasks.
 * Uses SDK query() which handles all auth types via getDefaultOptions().
 *
 * @param userMessage - The user's first message
 * @param language - Optional language override (defaults to user preference or '中文')
 * @returns Generated task title, or null if generation fails
 */
export async function generateSessionTitle(
  userMessage: string,
  language?: string
): Promise<string | null> {
  try {
    const effectiveLanguage = language || getCurrentLanguage();
    const userSnippet = userMessage.slice(0, 500);
    const prompt = getTitlePrompt(userSnippet, effectiveLanguage);

    const defaultOptions = getDefaultOptions();
    const options = {
      ...defaultOptions,
      model: resolveModelId(SUMMARIZATION_MODEL),
      maxTurns: 1,
    };

    let title = '';

    for await (const message of query({ prompt, options })) {
      if (message.type === 'assistant') {
        for (const block of message.message.content) {
          if (block.type === 'text') {
            title += block.text;
          }
        }
      }
    }

    const trimmed = title.trim();

    // Validate: reasonable length, not empty
    if (trimmed && trimmed.length > 0 && trimmed.length < 100) {
      return trimmed;
    }

    return null;
  } catch (error) {
    console.error('[title-generator] Failed to generate title:', error);
    return null;
  }
}

/**
 * Regenerate a session title based on recent messages.
 * Uses the most recent user messages to capture what the session has evolved into,
 * rather than just the initial topic.
 *
 * @param recentUserMessages - The last few user messages (most recent context)
 * @param lastAssistantResponse - The most recent assistant response
 * @param language - Optional language override (defaults to user preference or '中文')
 * @returns Generated title reflecting current session focus, or null if generation fails
 */
export async function regenerateSessionTitle(
  recentUserMessages: string[],
  lastAssistantResponse: string,
  language?: string
): Promise<string | null> {
  try {
    const effectiveLanguage = language || getCurrentLanguage();
    // Combine recent user messages, taking up to 300 chars from each
    const userContext = recentUserMessages
      .map((msg) => msg.slice(0, 300))
      .join('\n\n');
    const assistantSnippet = lastAssistantResponse.slice(0, 500);
    const prompt = getRegenerateTitlePrompt(userContext, assistantSnippet, effectiveLanguage);

    const defaultOptions = getDefaultOptions();
    const options = {
      ...defaultOptions,
      model: resolveModelId(SUMMARIZATION_MODEL),
      maxTurns: 1,
    };

    let title = '';

    for await (const message of query({ prompt, options })) {
      if (message.type === 'assistant') {
        for (const block of message.message.content) {
          if (block.type === 'text') {
            title += block.text;
          }
        }
      }
    }

    const trimmed = title.trim();

    if (trimmed && trimmed.length > 0 && trimmed.length < 100) {
      return trimmed;
    }

    return null;
  } catch (error) {
    console.error('[title-generator] Failed to regenerate title:', error);
    return null;
  }
}
