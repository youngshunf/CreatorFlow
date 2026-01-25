import { useAuthStore } from './auth'

const API_BASE = '/api'

class ApiError extends Error {
  constructor(public status: number, message: string) {
    super(message)
    this.name = 'ApiError'
  }
}

async function request<T>(path: string, options: RequestInit = {}): Promise<T> {
  const token = useAuthStore.getState().token
  
  const response = await fetch(`${API_BASE}${path}`, {
    ...options,
    headers: {
      'Content-Type': 'application/json',
      ...(token ? { Authorization: `Bearer ${token}` } : {}),
      ...options.headers,
    },
  })
  
  if (!response.ok) {
    const error = await response.json().catch(() => ({ error: 'Request failed' }))
    
    if (response.status === 401) {
      useAuthStore.getState().logout()
    }
    
    throw new ApiError(response.status, error.error || 'Request failed')
  }
  
  return response.json()
}

// Auth
export async function adminLogin(email: string, password: string) {
  return request<{ user: any; token: string }>('/auth/login', {
    method: 'POST',
    body: JSON.stringify({ email, password }),
  })
}

// Dashboard
export async function getDashboardStats() {
  return request<{
    userCount: number
    activeSubscriptions: number
    revenue: { mrr: number }
    recentOrders: any[]
  }>('/admin/stats')
}

// Users
export async function getUsers(params: { page?: number; limit?: number; search?: string }) {
  const query = new URLSearchParams()
  if (params.page) query.set('page', String(params.page))
  if (params.limit) query.set('limit', String(params.limit))
  if (params.search) query.set('search', params.search)
  
  return request<{
    users: any[]
    total: number
    page: number
    limit: number
  }>(`/admin/users?${query}`)
}

export async function getUser(userId: string) {
  return request<{ user: any }>(`/admin/users/${userId}`)
}

// Subscriptions
export async function getSubscriptions(params: { page?: number; limit?: number; status?: string; plan?: string }) {
  const query = new URLSearchParams()
  if (params.page) query.set('page', String(params.page))
  if (params.limit) query.set('limit', String(params.limit))
  if (params.status) query.set('status', params.status)
  if (params.plan) query.set('plan', params.plan)
  
  return request<{
    subscriptions: any[]
    total: number
    page: number
    limit: number
  }>(`/admin/subscriptions?${query}`)
}

export async function grantSubscription(userId: string, plan: string, months: number, reason: string) {
  return request<{ subscription: any }>('/admin/subscriptions/grant', {
    method: 'POST',
    body: JSON.stringify({ userId, plan, months, reason }),
  })
}

// Orders
export async function getOrders(params: { page?: number; limit?: number; status?: string; channel?: string }) {
  const query = new URLSearchParams()
  if (params.page) query.set('page', String(params.page))
  if (params.limit) query.set('limit', String(params.limit))
  if (params.status) query.set('status', params.status)
  if (params.channel) query.set('channel', params.channel)
  
  return request<{
    orders: any[]
    total: number
    page: number
    limit: number
  }>(`/admin/orders?${query}`)
}
