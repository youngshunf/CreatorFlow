import { Hono } from 'hono'
import { cors } from 'hono/cors'
import { logger } from 'hono/logger'
import { secureHeaders } from 'hono/secure-headers'
import { auth } from './routes/auth'
import { billing } from './routes/billing'
import { agent } from './routes/agent'
import { admin } from './routes/admin'

// Create Hono app with typed context
export type Env = {
  Variables: {
    userId: string
    userEmail: string
  }
}

const app = new Hono<Env>()

// Global middleware
app.use('*', logger())
app.use('*', secureHeaders())
app.use('*', cors({
  origin: ['http://localhost:5173', 'http://localhost:3001', 'sproutyai://'],
  credentials: true,
}))

// Health check
app.get('/health', (c) => c.json({ status: 'ok', timestamp: new Date().toISOString() }))

// API routes
app.route('/api/auth', auth)
app.route('/api/billing', billing)
app.route('/api/agent', agent)
app.route('/api/admin', admin)

// 404 handler
app.notFound((c) => c.json({ error: 'Not found' }, 404))

// Error handler
app.onError((err, c) => {
  console.error('Server error:', err)
  return c.json({ 
    error: err.message || 'Internal server error',
    ...(process.env.NODE_ENV === 'development' ? { stack: err.stack } : {})
  }, 500)
})

export { app }
