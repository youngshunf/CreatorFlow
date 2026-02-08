import { HeadlessRunner, type HeadlessConfig, type HeadlessEvent } from '@sprouty-ai/shared/headless'
import { WorkspaceService } from './workspace'
import { SubscriptionService } from './subscription'

export interface AgentExecuteOptions {
  sessionId?: string
  model?: string
  permissionPolicy?: 'deny-all' | 'allow-safe' | 'allow-all'
}

export interface AgentEvent {
  type: string
  [key: string]: any
}

export class AgentService {
  /**
   * Execute an agent query for a user
   * Reuses the HeadlessRunner from @sprouty-ai/shared
   */
  static async *execute(
    userId: string,
    prompt: string,
    options?: AgentExecuteOptions
  ): AsyncGenerator<AgentEvent> {
    // 1. Check user quota
    const hasQuota = await SubscriptionService.checkQuota(userId)
    if (!hasQuota) {
      yield {
        type: 'error',
        message: 'Usage quota exceeded. Please upgrade your plan.',
        code: 'QUOTA_EXCEEDED',
      }
      return
    }
    
    // 2. Ensure user workspace exists
    const workspacePath = await WorkspaceService.ensureUserWorkspace(userId)
    
    // 3. Create HeadlessRunner config
    const config: HeadlessConfig = {
      workspace: {
        id: userId,
        rootPath: workspacePath,
        name: 'Cloud Workspace',
        createdAt: Date.now(),
      },
      prompt,
      model: options?.model,
      sessionId: options?.sessionId,
      sessionResume: !!options?.sessionId,
      // Cloud default: allow-safe (read operations + safe bash commands)
      permissionPolicy: options?.permissionPolicy || 'allow-safe',
    }
    
    // 4. Create and run HeadlessRunner
    const runner = new HeadlessRunner(config)
    
    let usage: { inputTokens?: number; outputTokens?: number } | undefined
    
    try {
      for await (const event of runner.runStreaming()) {
        // Forward events to caller
        yield event as AgentEvent
        
        // Capture usage from complete event
        if (event.type === 'complete' && event.result?.usage) {
          usage = {
            inputTokens: event.result.usage.inputTokens,
            outputTokens: event.result.usage.outputTokens,
          }
        }
      }
      
      // Record usage after successful completion
      if (usage) {
        await SubscriptionService.recordUsage(userId, usage)
      }
    } catch (error) {
      yield {
        type: 'error',
        message: error instanceof Error ? error.message : 'Unknown error',
        code: 'EXECUTION_ERROR',
      }
    }
  }
  
  /**
   * Execute a simple query and return the result (non-streaming)
   */
  static async query(
    userId: string,
    prompt: string,
    options?: AgentExecuteOptions
  ): Promise<{ success: boolean; response?: string; error?: string }> {
    let response = ''
    let error: string | undefined
    
    for await (const event of this.execute(userId, prompt, options)) {
      if (event.type === 'text_complete') {
        response = event.text
      }
      if (event.type === 'error') {
        error = event.message
      }
    }
    
    if (error) {
      return { success: false, error }
    }
    
    return { success: true, response }
  }
}
