import { atom, useAtom } from 'jotai'

const STORAGE_KEY = 'admin-auth-token'

interface AuthState {
  token: string | null
  isAuthenticated: boolean
  isAdmin: boolean
}

function getInitialState(): AuthState {
  const token = localStorage.getItem(STORAGE_KEY)
  return {
    token,
    isAuthenticated: !!token,
    isAdmin: !!token,
  }
}

const authAtom = atom<AuthState>(getInitialState())

export function useAuthStore<T>(selector: (state: AuthState & { login: (token: string) => void; logout: () => void }) => T): T {
  const [state, setState] = useAtom(authAtom)
  
  const actions = {
    login: (token: string) => {
      // Write to localStorage immediately
      localStorage.setItem(STORAGE_KEY, token)
      setState({ token, isAuthenticated: true, isAdmin: true })
    },
    logout: () => {
      localStorage.removeItem(STORAGE_KEY)
      setState({ token: null, isAuthenticated: false, isAdmin: false })
    },
  }
  
  return selector({ ...state, ...actions })
}

// For API client (non-hook access)
export function getAuthState(): AuthState {
  const token = localStorage.getItem(STORAGE_KEY)
  return {
    token,
    isAuthenticated: !!token,
    isAdmin: !!token,
  }
}
