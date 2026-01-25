import { eq, like, or, sql, count } from 'drizzle-orm'
import { nanoid } from 'nanoid'
import { db, users, apiTokens, type User, type NewUser } from '../db'

// Simple password hashing using Bun's built-in crypto
async function hashPassword(password: string): Promise<string> {
  return await Bun.password.hash(password, {
    algorithm: 'bcrypt',
    cost: 10,
  })
}

async function verifyPassword(password: string, hash: string): Promise<boolean> {
  return await Bun.password.verify(password, hash)
}

export interface RegisterInput {
  email: string
  password: string
  name?: string
}

export interface GoogleOAuthProfile {
  email: string
  name?: string
  picture?: string
  sub: string  // Google user ID
}

export class UserService {
  /**
   * Register a new user with email/password
   */
  static async register(input: RegisterInput): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    // Check if email already exists
    const existing = await db.query.users.findFirst({
      where: eq(users.email, input.email),
    })
    
    if (existing) {
      throw new Error('Email already exists')
    }
    
    const passwordHash = await hashPassword(input.password)
    
    const newUser: NewUser = {
      id: nanoid(),
      email: input.email,
      passwordHash,
      name: input.name,
      provider: 'email',
      emailVerified: false,
    }
    
    const [user] = await db.insert(users).values(newUser).returning()
    return user!
  }
  
  /**
   * Login with email/password
   */
  static async login(email: string, password: string): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    const user = await db.query.users.findFirst({
      where: eq(users.email, email),
    })
    
    if (!user || !user.passwordHash) {
      throw new Error('Invalid credentials')
    }
    
    const valid = await verifyPassword(password, user.passwordHash)
    if (!valid) {
      throw new Error('Invalid credentials')
    }
    
    return user
  }
  
  /**
   * Handle Google OAuth callback
   */
  static async handleGoogleOAuth(code: string): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    // Exchange code for tokens
    const tokenResponse = await fetch('https://oauth2.googleapis.com/token', {
      method: 'POST',
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
      body: new URLSearchParams({
        code,
        client_id: process.env.GOOGLE_CLIENT_ID!,
        client_secret: process.env.GOOGLE_CLIENT_SECRET!,
        redirect_uri: `${process.env.APP_URL}/api/auth/oauth/google/callback`,
        grant_type: 'authorization_code',
      }),
    })
    
    if (!tokenResponse.ok) {
      throw new Error('Failed to exchange OAuth code')
    }
    
    const tokens = await tokenResponse.json()
    
    // Get user info
    const userInfoResponse = await fetch('https://www.googleapis.com/oauth2/v2/userinfo', {
      headers: { Authorization: `Bearer ${tokens.access_token}` },
    })
    
    if (!userInfoResponse.ok) {
      throw new Error('Failed to get user info')
    }
    
    const profile: GoogleOAuthProfile = await userInfoResponse.json()
    
    // Check if user exists
    let user = await db.query.users.findFirst({
      where: or(
        eq(users.email, profile.email),
        eq(users.providerId, profile.sub)
      ),
    })
    
    if (user) {
      // Update existing user
      const [updated] = await db.update(users)
        .set({
          name: user.name || profile.name,
          avatar: user.avatar || profile.picture,
          providerId: profile.sub,
          emailVerified: true,
          updatedAt: new Date(),
        })
        .where(eq(users.id, user.id))
        .returning()
      return updated!
    }
    
    // Create new user
    const newUser: NewUser = {
      id: nanoid(),
      email: profile.email,
      name: profile.name,
      avatar: profile.picture,
      provider: 'google',
      providerId: profile.sub,
      emailVerified: true,
    }
    
    const [created] = await db.insert(users).values(newUser).returning()
    return created!
  }
  
  /**
   * Get user by ID
   */
  static async getById(id: string): Promise<User | null> {
    if (!db) return null
    
    const user = await db.query.users.findFirst({
      where: eq(users.id, id),
    })
    
    return user ?? null
  }
  
  /**
   * Get user by API key
   */
  static async getByApiKey(apiKey: string): Promise<User | null> {
    if (!db) return null
    
    // Hash the API key to compare
    const tokenHash = await Bun.password.hash(apiKey, { algorithm: 'bcrypt', cost: 4 })
    
    const token = await db.query.apiTokens.findFirst({
      where: eq(apiTokens.tokenHash, tokenHash),
    })
    
    if (!token) return null
    
    // Update last used
    await db.update(apiTokens)
      .set({ lastUsedAt: new Date() })
      .where(eq(apiTokens.id, token.id))
    
    return this.getById(token.userId)
  }
  
  /**
   * Update user profile
   */
  static async update(id: string, updates: Partial<Pick<User, 'name' | 'avatar' | 'emailVerified'>>): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    const [user] = await db.update(users)
      .set({
        ...updates,
        updatedAt: new Date(),
      })
      .where(eq(users.id, id))
      .returning()
    
    if (!user) {
      throw new Error('User not found')
    }
    
    return user
  }
  
  /**
   * List users (admin)
   */
  static async list(options: { page: number; limit: number; search?: string }) {
    if (!db) return { users: [], total: 0 }
    
    const { page, limit, search } = options
    const offset = (page - 1) * limit
    
    let whereClause = undefined
    if (search) {
      whereClause = or(
        like(users.email, `%${search}%`),
        like(users.name, `%${search}%`)
      )
    }
    
    const [results, totalResult] = await Promise.all([
      db.query.users.findMany({
        where: whereClause,
        limit,
        offset,
        orderBy: (users, { desc }) => [desc(users.createdAt)],
      }),
      db.select({ count: count() }).from(users).where(whereClause),
    ])
    
    return {
      users: results,
      total: totalResult[0]?.count ?? 0,
      page,
      limit,
    }
  }
  
  /**
   * Count total users
   */
  static async count(): Promise<number> {
    if (!db) return 0
    
    const result = await db.select({ count: count() }).from(users)
    return result[0]?.count ?? 0
  }
}
