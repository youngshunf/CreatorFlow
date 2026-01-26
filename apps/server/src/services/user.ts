import { eq, like, or, sql, count } from 'drizzle-orm'
import { nanoid } from 'nanoid'
import { db, users, apiTokens, type User, type NewUser, type ClientType, type UserSource } from '../db'
import { SmsService } from './sms'

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
  source?: UserSource
  allowedClients?: ClientType[]
}

export interface PhoneLoginInput {
  phone: string
  code: string
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
    
    // Default: desktop registration only allows desktop access
    const source = input.source ?? 'desktop'
    const allowedClients = input.allowedClients ?? (source === 'admin' ? ['desktop', 'admin'] : ['desktop'])
    
    const newUser: NewUser = {
      id: nanoid(),
      email: input.email,
      passwordHash,
      name: input.name,
      provider: 'email',
      emailVerified: false,
      source,
      allowedClients,
    }
    
    const [user] = await db.insert(users).values(newUser).returning()
    return user!
  }
  
  /**
   * Login with email/password
   * @param client - which client is attempting to login
   */
  static async login(email: string, password: string, client: ClientType = 'desktop'): Promise<User> {
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
    
    // Check client access
    const allowedClients = user.allowedClients as ClientType[] | null
    if (!allowedClients?.includes(client)) {
      throw new Error('Client access denied')
    }
    
    return user
  }
  
  /**
   * Login or register with phone number (auto-registration)
   * If user doesn't exist, create new account
   * @param client - which client is attempting to login
   */
  static async loginOrRegisterByPhone(input: PhoneLoginInput, client: ClientType = 'desktop'): Promise<{ user: User; isNewUser: boolean }> {
    if (!db) throw new Error('Database not configured')
    
    // Verify SMS code first
    const isValid = await SmsService.verifyCode(input.phone, input.code)
    if (!isValid) {
      throw new Error('Invalid verification code')
    }
    
    // Check if user exists
    let user = await db.query.users.findFirst({
      where: eq(users.phone, input.phone),
    })
    
    if (user) {
      // Existing user - check client access
      const allowedClients = user.allowedClients as ClientType[] | null
      if (!allowedClients?.includes(client)) {
        throw new Error('Client access denied')
      }
      return { user, isNewUser: false }
    }
    
    // New user - auto register (desktop only by default)
    const newUser: NewUser = {
      id: nanoid(),
      phone: input.phone,
      provider: 'phone',
      phoneVerified: true,
      source: 'desktop',
      allowedClients: ['desktop'],
    }
    
    // If registering from admin, deny (admin users must be pre-created or granted access)
    if (client === 'admin') {
      throw new Error('Client access denied')
    }
    
    const [created] = await db.insert(users).values(newUser).returning()
    return { user: created!, isNewUser: true }
  }
  
  /**
   * Get user by phone number
   */
  static async getByPhone(phone: string): Promise<User | null> {
    if (!db) return null
    
    const user = await db.query.users.findFirst({
      where: eq(users.phone, phone),
    })
    
    return user ?? null
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
  
  /**
   * Grant client access to a user
   */
  static async grantClientAccess(userId: string, client: ClientType): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    const user = await this.getById(userId)
    if (!user) throw new Error('User not found')
    
    const currentClients = (user.allowedClients as ClientType[] | null) ?? ['desktop']
    if (currentClients.includes(client)) {
      return user // Already has access
    }
    
    const [updated] = await db.update(users)
      .set({
        allowedClients: [...currentClients, client],
        updatedAt: new Date(),
      })
      .where(eq(users.id, userId))
      .returning()
    
    return updated!
  }
  
  /**
   * Revoke client access from a user
   */
  static async revokeClientAccess(userId: string, client: ClientType): Promise<User> {
    if (!db) throw new Error('Database not configured')
    
    const user = await this.getById(userId)
    if (!user) throw new Error('User not found')
    
    const currentClients = (user.allowedClients as ClientType[] | null) ?? ['desktop']
    const newClients = currentClients.filter(c => c !== client)
    
    // Must have at least desktop access
    if (newClients.length === 0) {
      newClients.push('desktop')
    }
    
    const [updated] = await db.update(users)
      .set({
        allowedClients: newClients,
        updatedAt: new Date(),
      })
      .where(eq(users.id, userId))
      .returning()
    
    return updated!
  }
}
