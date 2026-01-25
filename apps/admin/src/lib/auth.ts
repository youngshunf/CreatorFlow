import { create } from 'zustand'
import { persist } from 'zustand/middleware'

interface AuthState {
  token: string | null
  isAuthenticated: boolean
  isAdmin: boolean
  login: (token: string) => void
  logout: () => void
}

export const useAuthStore = create<AuthState>()(
  persist(
    (set) => ({
      token: null,
      isAuthenticated: false,
      isAdmin: false,
      login: (token) => set({ token, isAuthenticated: true, isAdmin: true }),
      logout: () => set({ token: null, isAuthenticated: false, isAdmin: false }),
    }),
    {
      name: 'admin-auth',
    }
  )
)
