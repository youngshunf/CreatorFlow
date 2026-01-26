import { pgTable, text, timestamp, boolean, integer, jsonb, index } from 'drizzle-orm/pg-core'

// Client types
export type ClientType = 'desktop' | 'admin'
export type UserSource = 'desktop' | 'admin' | 'web'

// Users table
export const users = pgTable('users', {
  id: text('id').primaryKey(),
  email: text('email').unique(),  // null if phone registration
  phone: text('phone').unique(),  // null if email/OAuth registration
  passwordHash: text('password_hash'),  // null if OAuth only
  name: text('name'),
  avatar: text('avatar'),
  provider: text('provider'),  // 'email' | 'phone' | 'google' | 'apple'
  providerId: text('provider_id'),
  emailVerified: boolean('email_verified').default(false),
  phoneVerified: boolean('phone_verified').default(false),
  // Multi-client access control
  source: text('source').$type<UserSource>().default('desktop'),  // where the user registered
  allowedClients: jsonb('allowed_clients').$type<ClientType[]>().default(['desktop']),  // which clients user can access
  createdAt: timestamp('created_at').defaultNow().notNull(),
  updatedAt: timestamp('updated_at').defaultNow().notNull(),
}, (table) => ({
  emailIdx: index('users_email_idx').on(table.email),
  phoneIdx: index('users_phone_idx').on(table.phone),
  providerIdx: index('users_provider_idx').on(table.provider, table.providerId),
}))

// Subscriptions table
export const subscriptions = pgTable('subscriptions', {
  id: text('id').primaryKey(),
  userId: text('user_id').references(() => users.id).notNull(),
  plan: text('plan').notNull(),  // 'free' | 'pro' | 'team'
  status: text('status').notNull(),  // 'active' | 'canceled' | 'past_due' | 'canceling'
  stripeCustomerId: text('stripe_customer_id'),
  stripeSubscriptionId: text('stripe_subscription_id'),
  currentPeriodStart: timestamp('current_period_start'),
  currentPeriodEnd: timestamp('current_period_end'),
  cancelAtPeriodEnd: boolean('cancel_at_period_end').default(false),
  // Manual grants (for support/promos)
  grantedBy: text('granted_by'),
  grantReason: text('grant_reason'),
  createdAt: timestamp('created_at').defaultNow().notNull(),
  updatedAt: timestamp('updated_at').defaultNow().notNull(),
}, (table) => ({
  userIdIdx: index('subscriptions_user_id_idx').on(table.userId),
  statusIdx: index('subscriptions_status_idx').on(table.status),
  stripeSubIdx: index('subscriptions_stripe_sub_idx').on(table.stripeSubscriptionId),
}))

// API Tokens table (for programmatic access)
export const apiTokens = pgTable('api_tokens', {
  id: text('id').primaryKey(),
  userId: text('user_id').references(() => users.id).notNull(),
  name: text('name').notNull(),
  tokenHash: text('token_hash').notNull(),  // hashed token
  lastUsedAt: timestamp('last_used_at'),
  expiresAt: timestamp('expires_at'),
  createdAt: timestamp('created_at').defaultNow().notNull(),
}, (table) => ({
  userIdIdx: index('api_tokens_user_id_idx').on(table.userId),
  tokenHashIdx: index('api_tokens_token_hash_idx').on(table.tokenHash),
}))

// Usage tracking table
export const usageRecords = pgTable('usage_records', {
  id: text('id').primaryKey(),
  userId: text('user_id').references(() => users.id).notNull(),
  type: text('type').notNull(),  // 'message' | 'tokens'
  amount: integer('amount').notNull(),
  metadata: jsonb('metadata'),  // { model, sessionId, etc. }
  createdAt: timestamp('created_at').defaultNow().notNull(),
}, (table) => ({
  userIdIdx: index('usage_records_user_id_idx').on(table.userId),
  createdAtIdx: index('usage_records_created_at_idx').on(table.createdAt),
  typeIdx: index('usage_records_type_idx').on(table.type),
}))

// Payment orders table
export const orders = pgTable('orders', {
  id: text('id').primaryKey(),
  orderNo: text('order_no').notNull().unique(),  // 商户订单号
  userId: text('user_id').references(() => users.id).notNull(),
  plan: text('plan').notNull(),  // 'pro' | 'team'
  period: text('period').notNull(),  // 'monthly' | 'yearly'
  amount: integer('amount').notNull(),  // 金额（分）
  channel: text('channel').notNull(),  // 'wechat' | 'alipay'
  method: text('method').notNull(),  // 'native' | 'jsapi' | 'h5' | 'page' | 'wap'
  status: text('status').notNull(),  // 'pending' | 'paid' | 'closed' | 'expired'
  transactionId: text('transaction_id'),  // 第三方交易号
  paidAt: timestamp('paid_at'),
  expiredAt: timestamp('expired_at'),
  metadata: jsonb('metadata'),  // 额外信息
  createdAt: timestamp('created_at').defaultNow().notNull(),
  updatedAt: timestamp('updated_at').defaultNow().notNull(),
}, (table) => ({
  orderNoIdx: index('orders_order_no_idx').on(table.orderNo),
  userIdIdx: index('orders_user_id_idx').on(table.userId),
  statusIdx: index('orders_status_idx').on(table.status),
  channelIdx: index('orders_channel_idx').on(table.channel),
  createdAtIdx: index('orders_created_at_idx').on(table.createdAt),
}))

// Audit log table (for admin actions)
export const auditLogs = pgTable('audit_logs', {
  id: text('id').primaryKey(),
  actorId: text('actor_id').references(() => users.id),
  action: text('action').notNull(),  // 'subscription.grant', 'subscription.revoke', etc.
  targetId: text('target_id'),
  targetType: text('target_type'),  // 'user', 'subscription'
  details: jsonb('details'),
  createdAt: timestamp('created_at').defaultNow().notNull(),
}, (table) => ({
  actorIdIdx: index('audit_logs_actor_id_idx').on(table.actorId),
  actionIdx: index('audit_logs_action_idx').on(table.action),
  createdAtIdx: index('audit_logs_created_at_idx').on(table.createdAt),
}))

// Type exports
export type User = typeof users.$inferSelect
export type NewUser = typeof users.$inferInsert
export type Subscription = typeof subscriptions.$inferSelect
export type NewSubscription = typeof subscriptions.$inferInsert
export type ApiToken = typeof apiTokens.$inferSelect
export type NewApiToken = typeof apiTokens.$inferInsert
export type UsageRecord = typeof usageRecords.$inferSelect
export type NewUsageRecord = typeof usageRecords.$inferInsert
export type AuditLog = typeof auditLogs.$inferSelect
export type NewAuditLog = typeof auditLogs.$inferInsert
export type Order = typeof orders.$inferSelect
export type NewOrder = typeof orders.$inferInsert
