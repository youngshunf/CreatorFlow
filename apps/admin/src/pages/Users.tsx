import { useState } from 'react'
import { useQuery } from '@tanstack/react-query'
import { getUsers } from '../lib/api'
import { format } from 'date-fns'

export function Users() {
  const [page, setPage] = useState(1)
  const [search, setSearch] = useState('')
  const limit = 20
  
  const { data, isLoading, error } = useQuery({
    queryKey: ['users', page, search],
    queryFn: () => getUsers({ page, limit, search: search || undefined }),
  })
  
  const totalPages = Math.ceil((data?.total ?? 0) / limit)
  
  return (
    <div>
      <h1 className="text-2xl font-bold mb-8">用户管理</h1>
      
      {/* Search */}
      <div className="mb-6">
        <input
          type="text"
          placeholder="搜索邮箱或用户名..."
          value={search}
          onChange={(e) => {
            setSearch(e.target.value)
            setPage(1)
          }}
          className="w-full max-w-md px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
        />
      </div>
      
      {/* Table */}
      <div className="bg-white rounded-lg shadow overflow-hidden">
        <table className="w-full">
          <thead className="bg-gray-50">
            <tr className="text-left text-sm text-gray-500">
              <th className="px-6 py-3">用户</th>
              <th className="px-6 py-3">邮箱</th>
              <th className="px-6 py-3">登录方式</th>
              <th className="px-6 py-3">套餐</th>
              <th className="px-6 py-3">注册时间</th>
              <th className="px-6 py-3">操作</th>
            </tr>
          </thead>
          <tbody>
            {isLoading && (
              <tr>
                <td colSpan={6} className="px-6 py-8 text-center text-gray-500">
                  加载中...
                </td>
              </tr>
            )}
            {error && (
              <tr>
                <td colSpan={6} className="px-6 py-8 text-center text-red-600">
                  加载失败: {(error as Error).message}
                </td>
              </tr>
            )}
            {data?.users?.map((user: any) => (
              <tr key={user.id} className="border-t hover:bg-gray-50">
                <td className="px-6 py-4">
                  <div className="flex items-center gap-3">
                    {user.avatar ? (
                      <img src={user.avatar} alt="" className="w-8 h-8 rounded-full" />
                    ) : (
                      <div className="w-8 h-8 rounded-full bg-gray-200 flex items-center justify-center">
                        {user.name?.[0] || user.email[0].toUpperCase()}
                      </div>
                    )}
                    <span className="font-medium">{user.name || '未设置'}</span>
                  </div>
                </td>
                <td className="px-6 py-4 text-sm">{user.email}</td>
                <td className="px-6 py-4 text-sm">{user.provider || 'email'}</td>
                <td className="px-6 py-4">
                  <PlanBadge plan={user.subscription?.plan || 'free'} />
                </td>
                <td className="px-6 py-4 text-sm text-gray-500">
                  {format(new Date(user.createdAt), 'yyyy-MM-dd')}
                </td>
                <td className="px-6 py-4">
                  <button className="text-blue-600 hover:underline text-sm">
                    查看详情
                  </button>
                </td>
              </tr>
            ))}
            {data?.users?.length === 0 && (
              <tr>
                <td colSpan={6} className="px-6 py-8 text-center text-gray-500">
                  暂无用户
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
      
      {/* Pagination */}
      {totalPages > 1 && (
        <div className="mt-4 flex items-center justify-between">
          <p className="text-sm text-gray-500">
            共 {data?.total ?? 0} 条记录，第 {page} / {totalPages} 页
          </p>
          <div className="flex gap-2">
            <button
              onClick={() => setPage((p) => Math.max(1, p - 1))}
              disabled={page === 1}
              className="px-3 py-1 border rounded disabled:opacity-50"
            >
              上一页
            </button>
            <button
              onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
              disabled={page === totalPages}
              className="px-3 py-1 border rounded disabled:opacity-50"
            >
              下一页
            </button>
          </div>
        </div>
      )}
    </div>
  )
}

function PlanBadge({ plan }: { plan: string }) {
  const planConfig: Record<string, { label: string; className: string }> = {
    free: { label: 'Free', className: 'bg-gray-100 text-gray-800' },
    pro: { label: 'Pro', className: 'bg-blue-100 text-blue-800' },
    team: { label: 'Team', className: 'bg-purple-100 text-purple-800' },
  }
  
  const config = planConfig[plan] || planConfig.free
  
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config.className}`}>
      {config.label}
    </span>
  )
}
