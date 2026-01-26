/**
 * 初始化管理员账号脚本
 * 运行: bun run scripts/seed-admin.ts
 */

import { db } from '../src/db'
import { users, subscriptions } from '../src/db/schema'
import { eq, or } from 'drizzle-orm'

// 默认管理员配置
const ADMINS = [
  {
    phone: '18611348367',
    password: '123456',
    name: '超级管理员',
    provider: 'phone' as const,
  },
  {
    email: 'admin@creatorflow.top',
    password: '123456',
    name: 'Admin',
    provider: 'email' as const,
  },
]

async function seedAdmin() {
  if (!db) {
    console.error('Database not initialized')
    process.exit(1)
  }
  
  console.log('Seeding admin users...')
  
  for (const admin of ADMINS) {
    const identifier = admin.email || admin.phone
    const identifierField = admin.email ? 'email' : 'phone'
    
    // Check if admin already exists
    const whereClause = admin.email 
      ? eq(users.email, admin.email)
      : eq(users.phone, admin.phone!)
    
    const [existing] = await db.select().from(users).where(whereClause).limit(1)
    
    if (existing) {
      // Update existing user to have admin access if not already
      const allowedClients = existing.allowedClients as string[] | null
      if (!allowedClients?.includes('admin')) {
        await db.update(users)
          .set({ 
            allowedClients: ['desktop', 'admin'],
            source: 'admin',
          })
          .where(eq(users.id, existing.id))
        console.log(`✓ Updated existing user with admin access: ${identifier}`)
      } else {
        console.log(`Admin user already exists: ${identifier}`)
      }
      continue
    }
    
    // Hash password with bcrypt
    const passwordHash = await Bun.password.hash(admin.password, {
      algorithm: 'bcrypt',
      cost: 10,
    })
    
    const userId = crypto.randomUUID()
    
    // Create admin user with full access
    await db.insert(users).values({
      id: userId,
      email: admin.email,
      phone: admin.phone,
      passwordHash,
      name: admin.name,
      provider: admin.provider,
      emailVerified: admin.email ? true : false,
      phoneVerified: admin.phone ? true : false,
      source: 'admin',
      allowedClients: ['desktop', 'admin'],
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
    
    console.log(`✓ Admin user created: ${identifier}`)
    console.log(`  ${identifierField}: ${identifier}`)
    console.log(`  Password: ${admin.password}`)
  }
  
  console.log('\n✓ All admin users seeded successfully')
}

seedAdmin()
  .then(() => process.exit(0))
  .catch((err) => {
    console.error('Failed to seed admin:', err)
    process.exit(1)
  })
