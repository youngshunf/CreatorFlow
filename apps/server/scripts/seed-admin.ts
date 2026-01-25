/**
 * 初始化管理员账号脚本
 * 运行: bun run scripts/seed-admin.ts
 */

import { db } from '../src/db'
import { users, subscriptions } from '../src/db/schema'
import { eq } from 'drizzle-orm'

const ADMIN_EMAIL = 'admin@creatorflow.top'
const ADMIN_PASSWORD = '123456'

async function seedAdmin() {
  if (!db) {
    console.error('Database not initialized')
    process.exit(1)
  }
  
  console.log('Seeding admin user...')
  
  // Check if admin already exists
  const [existing] = await db.select().from(users).where(eq(users.email, ADMIN_EMAIL)).limit(1)
  
  if (existing) {
    console.log('Admin user already exists:', existing.email)
    return
  }
  
  // Hash password with bcrypt
  const passwordHash = await Bun.password.hash(ADMIN_PASSWORD, {
    algorithm: 'bcrypt',
    cost: 10,
  })
  
  const userId = crypto.randomUUID()
  
  // Create admin user
  await db.insert(users).values({
    id: userId,
    email: ADMIN_EMAIL,
    passwordHash,
    name: 'Admin',
    provider: 'email',
    emailVerified: true,
  })
  
  // Create team subscription for admin
  await db.insert(subscriptions).values({
    id: crypto.randomUUID(),
    userId,
    plan: 'team',
    status: 'active',
    currentPeriodStart: new Date(),
    currentPeriodEnd: new Date(Date.now() + 365 * 24 * 60 * 60 * 1000), // 1 year
    grantedBy: 'system',
    grantReason: 'Admin account',
  })
  
  console.log('✓ Admin user created successfully')
  console.log('  Email:', ADMIN_EMAIL)
  console.log('  Password:', ADMIN_PASSWORD)
}

seedAdmin()
  .then(() => process.exit(0))
  .catch((err) => {
    console.error('Failed to seed admin:', err)
    process.exit(1)
  })
