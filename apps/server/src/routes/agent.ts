import { Hono } from 'hono'
import { zValidator } from '@hono/zod-validator'
import { z } from 'zod'
import { authMiddleware, verifyToken } from '../middleware/auth'
import { quotaMiddleware } from '../middleware/quota'
import { AgentService } from '../services/agent'
import { WorkspaceService } from '../services/workspace'
import type { Env } from '../app'

const agent = new Hono<Env>()

// WebSocket streaming chat endpoint
agent.get('/chat', async (c) => {
  // Upgrade to WebSocket
  const upgradeHeader = c.req.header('Upgrade')
  if (upgradeHeader !== 'websocket') {
    return c.json({ error: 'Expected WebSocket upgrade' }, 426)
  }
  
  // Get the raw request for Bun's WebSocket upgrade
  const server = (c.env as any)?.server
  if (!server) {
    return c.json({ error: 'WebSocket not supported in this environment' }, 500)
  }
  
  const success = server.upgrade(c.req.raw, {
    data: { /* initial data */ },
  })
  
  if (!success) {
    return c.json({ error: 'WebSocket upgrade failed' }, 500)
  }
  
  // Return empty response - connection is upgraded
  return new Response(null, { status: 101 })
})

// REST API for non-streaming queries
agent.post(
  '/query',
  authMiddleware,
  quotaMiddleware,
  zValidator('json', z.object({
    prompt: z.string().min(1).max(100000),
    sessionId: z.string().optional(),
    model: z.string().optional(),
  })),
  async (c) => {
    const userId = c.get('userId')
    const { prompt, sessionId, model } = c.req.valid('json')
    
    const events: any[] = []
    let response = ''
    let usage = null
    
    for await (const event of AgentService.execute(userId, prompt, { sessionId, model })) {
      events.push(event)
      
      if (event.type === 'text_complete') {
        response = event.text
      }
      if (event.type === 'complete' && event.result?.usage) {
        usage = event.result.usage
      }
    }
    
    const completeEvent = events.find(e => e.type === 'complete')
    
    return c.json({
      success: completeEvent?.result?.success ?? false,
      response,
      sessionId: completeEvent?.result?.sessionId,
      usage,
    })
  }
)

// Server-Sent Events streaming (alternative to WebSocket)
agent.post(
  '/stream',
  authMiddleware,
  quotaMiddleware,
  zValidator('json', z.object({
    prompt: z.string().min(1).max(100000),
    sessionId: z.string().optional(),
    model: z.string().optional(),
  })),
  async (c) => {
    const userId = c.get('userId')
    const { prompt, sessionId, model } = c.req.valid('json')
    
    // Set up SSE headers
    c.header('Content-Type', 'text/event-stream')
    c.header('Cache-Control', 'no-cache')
    c.header('Connection', 'keep-alive')
    
    const stream = new ReadableStream({
      async start(controller) {
        const encoder = new TextEncoder()
        
        try {
          for await (const event of AgentService.execute(userId, prompt, { sessionId, model })) {
            const data = `data: ${JSON.stringify(event)}\n\n`
            controller.enqueue(encoder.encode(data))
          }
          
          controller.enqueue(encoder.encode('data: [DONE]\n\n'))
          controller.close()
        } catch (error) {
          const errorEvent = `data: ${JSON.stringify({ 
            type: 'error', 
            message: error instanceof Error ? error.message : 'Unknown error' 
          })}\n\n`
          controller.enqueue(encoder.encode(errorEvent))
          controller.close()
        }
      },
    })
    
    return new Response(stream, {
      headers: {
        'Content-Type': 'text/event-stream',
        'Cache-Control': 'no-cache',
        'Connection': 'keep-alive',
      },
    })
  }
)

// List sessions
agent.get('/sessions', authMiddleware, async (c) => {
  const userId = c.get('userId')
  
  const workspacePath = await WorkspaceService.ensureUserWorkspace(userId)
  const { listSessions } = await import('@creator-flow/shared/sessions')
  const sessions = listSessions(workspacePath)
  
  return c.json({
    sessions: sessions.map(s => ({
      id: s.id,
      title: (s as any).title || `Session ${s.id.slice(0, 8)}`,
      createdAt: s.createdAt,
      lastUsedAt: (s as any).lastUsedAt || s.createdAt,
      messageCount: (s as any).messageCount || 0,
    })),
  })
})

// Get session details
agent.get('/sessions/:sessionId', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const sessionId = c.req.param('sessionId')
  
  const workspacePath = await WorkspaceService.ensureUserWorkspace(userId)
  const { loadSession } = await import('@creator-flow/shared/sessions')
  const session = loadSession(workspacePath, sessionId)
  
  if (!session) {
    return c.json({ error: 'Session not found' }, 404)
  }
  
  return c.json({ session })
})

// Delete session
agent.delete('/sessions/:sessionId', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const sessionId = c.req.param('sessionId')
  
  const workspacePath = await WorkspaceService.ensureUserWorkspace(userId)
  const { deleteSession } = await import('@creator-flow/shared/sessions')
  
  await deleteSession(workspacePath, sessionId)
  
  return c.json({ success: true })
})

// Get workspace info
agent.get('/workspace', authMiddleware, async (c) => {
  const userId = c.get('userId')
  
  const workspacePath = await WorkspaceService.ensureUserWorkspace(userId)
  const { loadWorkspaceConfig } = await import('@creator-flow/shared/workspaces')
  
  try {
    const config = loadWorkspaceConfig(workspacePath)
    return c.json({
      workspace: {
        path: workspacePath,
        config,
      },
    })
  } catch {
    return c.json({
      workspace: {
        path: workspacePath,
        config: null,
      },
    })
  }
})

export { agent }

// WebSocket handler for Bun (called from index.ts)
export async function handleAgentWebSocket(ws: any, message: string) {
  try {
    const data = JSON.parse(message)
    
    // Verify token
    const payload = await verifyToken(data.token)
    if (!payload?.userId) {
      ws.send(JSON.stringify({ type: 'error', message: 'Unauthorized' }))
      ws.close()
      return
    }
    
    const userId = payload.userId as string
    
    // Execute agent
    for await (const event of AgentService.execute(
      userId,
      data.prompt,
      { sessionId: data.sessionId, model: data.model }
    )) {
      ws.send(JSON.stringify(event))
    }
  } catch (error) {
    ws.send(JSON.stringify({ 
      type: 'error', 
      message: error instanceof Error ? error.message : 'Unknown error' 
    }))
  }
}
