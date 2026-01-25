import { Link, useLocation } from 'react-router-dom'
import { useAuthStore } from '../lib/auth'
import clsx from 'clsx'

const navItems = [
  { path: '/', label: '仪表盘', icon: '📊' },
  { path: '/users', label: '用户管理', icon: '👥' },
  { path: '/subscriptions', label: '订阅管理', icon: '💳' },
  { path: '/orders', label: '订单管理', icon: '📋' },
]

export function Layout({ children }: { children: React.ReactNode }) {
  const location = useLocation()
  const logout = useAuthStore((s) => s.logout)
  
  return (
    <div className="min-h-screen flex">
      {/* Sidebar */}
      <aside className="w-64 bg-gray-900 text-white">
        <div className="p-4">
          <h1 className="text-xl font-bold">CreatorFlow</h1>
          <p className="text-sm text-gray-400">管理后台</p>
        </div>
        
        <nav className="mt-4">
          {navItems.map((item) => (
            <Link
              key={item.path}
              to={item.path}
              className={clsx(
                'flex items-center gap-3 px-4 py-3 transition-colors',
                location.pathname === item.path
                  ? 'bg-gray-800 text-white'
                  : 'text-gray-400 hover:bg-gray-800 hover:text-white'
              )}
            >
              <span>{item.icon}</span>
              <span>{item.label}</span>
            </Link>
          ))}
        </nav>
        
        <div className="absolute bottom-0 w-64 p-4">
          <button
            onClick={logout}
            className="w-full px-4 py-2 text-sm text-gray-400 hover:text-white hover:bg-gray-800 rounded transition-colors"
          >
            退出登录
          </button>
        </div>
      </aside>
      
      {/* Main content */}
      <main className="flex-1 p-8 overflow-auto">
        {children}
      </main>
    </div>
  )
}
