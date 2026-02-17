/**
 * PromptHandler - Processes prompt hooks for App events
 *
 * Subscribes to App events and collects prompt hooks to be executed.
 * Prompts are queued and delivered via callback for the caller to execute.
 * Execution records are written to the database for UI display.
 */

import { randomUUID } from 'crypto';
import { join } from 'path';
import { createLogger } from '../../utils/debug.ts';
import type { EventBus, BaseEventPayload } from '../event-bus.ts';
import type { HookHandler, PromptHandlerOptions, HooksConfigProvider } from './types.ts';
import type { HookEvent, PromptHookDefinition, PendingPrompt, AppEvent } from '../types.ts';
import { matcherMatches, buildEnvFromPayload, expandEnvVars, parsePromptReferences } from '../utils.ts';
import { getCachedConnection } from '../../db/connection.ts';
import type { ScheduledTaskExecution } from '../../db/types.ts';

const log = createLogger('prompt-handler');

// App events that support prompt hooks
const APP_EVENTS: AppEvent[] = [
  'LabelAdd', 'LabelRemove', 'LabelConfigChange',
  'PermissionModeChange', 'FlagChange', 'SessionStatusChange',
  'SchedulerTick'
];

// ============================================================================
// PromptHandler Implementation
// ============================================================================

export class PromptHandler implements HookHandler {
  private readonly options: PromptHandlerOptions;
  private readonly configProvider: HooksConfigProvider;
  private bus: EventBus | null = null;
  private boundHandler: ((event: HookEvent, payload: BaseEventPayload) => Promise<void>) | null = null;

  constructor(options: PromptHandlerOptions, configProvider: HooksConfigProvider) {
    this.options = options;
    this.configProvider = configProvider;
  }

  /**
   * Subscribe to App events on the bus.
   */
  subscribe(bus: EventBus): void {
    this.bus = bus;
    this.boundHandler = this.handleEvent.bind(this);
    bus.onAny(this.boundHandler);
    log.debug(`[PromptHandler] Subscribed to event bus`);
  }

  /**
   * Handle an event by processing matching prompt hooks.
   */
  private async handleEvent(event: HookEvent, payload: BaseEventPayload): Promise<void> {
    // Only process App events for prompt hooks
    if (!APP_EVENTS.includes(event as AppEvent)) {
      return;
    }

    const matchers = this.configProvider.getMatchersForEvent(event);
    if (matchers.length === 0) return;

    // Find matching prompt hooks
    const promptHooks: Array<{ prompt: PromptHookDefinition; labels?: string[]; permissionMode?: 'safe' | 'ask' | 'allow-all' }> = [];

    for (const matcher of matchers) {
      if (!matcherMatches(matcher, event, payload as unknown as Record<string, unknown>)) continue;

      for (const hook of matcher.hooks) {
        if (hook.type === 'prompt') {
          promptHooks.push({ prompt: hook, labels: matcher.labels, permissionMode: matcher.permissionMode });
        }
      }
    }

    if (promptHooks.length === 0) return;

    log.debug(`[PromptHandler] Processing ${promptHooks.length} prompts for ${event}`);

    // Build environment variables
    const env = buildEnvFromPayload(event, payload);

    // Process prompts
    const pendingPrompts: PendingPrompt[] = [];

    for (const { prompt, labels, permissionMode } of promptHooks) {
      // Expand environment variables in the prompt
      const expandedPrompt = expandEnvVars(prompt.prompt, env);

      // Parse references
      const references = parsePromptReferences(expandedPrompt);

      // Expand labels
      const expandedLabels = labels?.map(label => expandEnvVars(label, env));

      pendingPrompts.push({
        sessionId: this.options.sessionId,
        prompt: expandedPrompt,
        mentions: references.mentions,
        labels: expandedLabels,
        permissionMode,
      });
    }

    // Deliver prompts via callback
    if (pendingPrompts.length > 0 && this.options.onPromptsReady) {
      const startTime = Date.now();
      log.debug(`[PromptHandler] Delivering ${pendingPrompts.length} prompts`);

      // 写入数据库执行记录
      const dbPath = join(this.options.workspaceRootPath, '.sprouty-ai', 'db', 'creator.db');
      const db = getCachedConnection(dbPath);

      if (db) {
        for (const p of pendingPrompts) {
          const executionId = randomUUID();
          const now = new Date().toISOString();

          const record: Omit<ScheduledTaskExecution, 'created_at'> = {
            id: executionId,
            task_id: 'prompt-hook', // 可以从 matcher 中提取更具体的 task_id
            task_name: event,
            trigger_event: event,
            trigger_time: now,
            started_at: now,
            completed_at: null,
            status: 'pending',
            result_summary: p.prompt.substring(0, 200),
            result_detail: JSON.stringify({ prompt: p.prompt, labels: p.labels }),
            error_message: null,
            duration_ms: null,
          };

          try {
            db.run(
              `INSERT INTO scheduled_task_executions
               (id, task_id, task_name, trigger_event, trigger_time, started_at, completed_at, status, result_summary, result_detail, error_message, duration_ms)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`,
              [
                record.id,
                record.task_id,
                record.task_name,
                record.trigger_event,
                record.trigger_time,
                record.started_at,
                record.completed_at,
                record.status,
                record.result_summary,
                record.result_detail,
                record.error_message,
                record.duration_ms,
              ]
            );
          } catch (err) {
            log.error('[PromptHandler] Failed to write execution record:', err);
          }
        }
      }

      this.options.onPromptsReady(pendingPrompts);

      // 保留 events.jsonl 日志用于调试
      for (const p of pendingPrompts) {
        this.options.eventLogger?.logResult({
          event,
          hookType: 'prompt',
          prompt: p.prompt,
          success: true,
          durationMs: Date.now() - startTime,
          workspaceId: this.options.workspaceId,
          sessionId: this.options.sessionId,
        });
      }
    }
  }

  /**
   * Clean up resources.
   */
  dispose(): void {
    if (this.bus && this.boundHandler) {
      this.bus.offAny(this.boundHandler);
      this.boundHandler = null;
    }
    this.bus = null;
    log.debug(`[PromptHandler] Disposed`);
  }
}
