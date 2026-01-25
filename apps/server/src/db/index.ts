import { drizzle } from 'drizzle-orm/postgres-js'
import postgres from 'postgres'
import * as schema from './schema'

// Create postgres connection
const connectionString = process.env.DATABASE_URL

if (!connectionString) {
  console.warn('DATABASE_URL not set, using mock database')
}

// Create postgres client
const client = connectionString 
  ? postgres(connectionString, { 
      max: 10,
      idle_timeout: 20,
      connect_timeout: 10,
    })
  : null

// Create drizzle instance
export const db = client 
  ? drizzle(client, { schema })
  : null

// Export schema for convenience
export * from './schema'

// Health check function
export async function checkDatabaseConnection(): Promise<boolean> {
  if (!db || !client) {
    return false
  }
  
  try {
    await client`SELECT 1`
    return true
  } catch (error) {
    console.error('Database connection failed:', error)
    return false
  }
}

// Graceful shutdown
export async function closeDatabaseConnection(): Promise<void> {
  if (client) {
    await client.end()
  }
}
