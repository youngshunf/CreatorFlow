import { Hono } from 'hono'
import { zValidator } from '@hono/zod-validator'
import { z } from 'zod'
import { sign } from 'hono/jwt'
import { UserService } from '../services/user'
import { authMiddleware } from '../middleware/auth'
import type { Env } from '../app'

const auth = new Hono<Env>()

// Register
auth.post(
  '/register',
  zValidator('json', z.object({
    email: z.string().email(),
    password: z.string().min(8),
    name: z.string().optional(),
  })),
  async (c) => {
    const { email, password, name } = c.req.valid('json')
    
    try {
      const user = await UserService.register({ email, password, name })
      const token = await sign(
        { userId: user.id, email: user.email },
        process.env.JWT_SECRET!
      )
      
      return c.json({
        user: {
          id: user.id,
          email: user.email,
          name: user.name,
        },
        token,
      })
    } catch (error) {
      if (error instanceof Error && error.message === 'Email already exists') {
        return c.json({ error: 'Email already exists' }, 409)
      }
      throw error
    }
  }
)

// Login
auth.post(
  '/login',
  zValidator('json', z.object({
    email: z.string().email(),
    password: z.string(),
  })),
  async (c) => {
    const { email, password } = c.req.valid('json')
    
    try {
      const user = await UserService.login(email, password)
      const token = await sign(
        { userId: user.id, email: user.email },
        process.env.JWT_SECRET!
      )
      
      return c.json({
        user: {
          id: user.id,
          email: user.email,
          name: user.name,
        },
        token,
      })
    } catch (error) {
      if (error instanceof Error && error.message === 'Invalid credentials') {
        return c.json({ error: 'Invalid email or password' }, 401)
      }
      throw error
    }
  }
)

// Get current user
auth.get('/me', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const user = await UserService.getById(userId)
  
  if (!user) {
    return c.json({ error: 'User not found' }, 404)
  }
  
  return c.json({
    user: {
      id: user.id,
      email: user.email,
      name: user.name,
      avatar: user.avatar,
      emailVerified: user.emailVerified,
      createdAt: user.createdAt,
    },
  })
})

// Update profile
auth.patch(
  '/me',
  authMiddleware,
  zValidator('json', z.object({
    name: z.string().optional(),
    avatar: z.string().url().optional(),
  })),
  async (c) => {
    const userId = c.get('userId')
    const updates = c.req.valid('json')
    
    const user = await UserService.update(userId, updates)
    
    return c.json({
      user: {
        id: user.id,
        email: user.email,
        name: user.name,
        avatar: user.avatar,
      },
    })
  }
)

// OAuth: Google
auth.get('/oauth/google', async (c) => {
  const redirectUri = `${process.env.APP_URL}/api/auth/oauth/google/callback`
  const scope = encodeURIComponent('openid email profile')
  const url = `https://accounts.google.com/o/oauth2/v2/auth?` +
    `client_id=${process.env.GOOGLE_CLIENT_ID}` +
    `&redirect_uri=${encodeURIComponent(redirectUri)}` +
    `&response_type=code` +
    `&scope=${scope}`
  
  return c.redirect(url)
})

auth.get('/oauth/google/callback', async (c) => {
  const code = c.req.query('code')
  
  if (!code) {
    return c.json({ error: 'Missing authorization code' }, 400)
  }
  
  try {
    const user = await UserService.handleGoogleOAuth(code)
    const token = await sign(
      { userId: user.id, email: user.email },
      process.env.JWT_SECRET!
    )
    
    // Redirect to desktop app via deep link
    const deepLinkScheme = process.env.DESKTOP_DEEPLINK_SCHEME || 'creatorflow'
    return c.redirect(`${deepLinkScheme}://auth?token=${token}`)
  } catch (error) {
    console.error('Google OAuth error:', error)
    return c.json({ error: 'OAuth authentication failed' }, 500)
  }
})

// Refresh token
auth.post('/refresh', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const userEmail = c.get('userEmail')
  
  const token = await sign(
    { userId, email: userEmail },
    process.env.JWT_SECRET!
  )
  
  return c.json({ token })
})

// Logout (client-side token removal, but we can track it)
auth.post('/logout', authMiddleware, async (c) => {
  // In a more complex setup, you might want to:
  // - Add the token to a blacklist
  // - Clear refresh tokens from database
  return c.json({ success: true })
})

export { auth }
