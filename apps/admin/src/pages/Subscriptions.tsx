import { useState } from 'react'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { getSubscriptions, grantSubscription } from '../lib/api'
import { format } from 'date-fns'

export function Subscriptions() {
  const [page, setPage] = useState(1)
  const [statusFilter, setStatusFilter] = useState('')
  const [planFilter, setPlanFilter] = useState('')
  const [showGrantModal, setShowGrantModal] = useState(false)
  const limit = 20
  
  const queryClient = useQueryClient()
  
  const { data, isLoading, error } = useQuery({
    queryKey: ['subscriptions', page, statusFilter, planFilter],
    queryFn: () => getSubscriptions({
      page,
      limit,
      status: statusFilter || undefined,
      plan: planFilter || undefined,
    }),
  })
  
  const totalPages = Math.ceil((data?.total ?? 0) / limit)
  
  return (
    <div>
      <div className="flex items-center justify-between mb-8">
        <h1 className="text-2xl font-bold">订阅管理</h1>
        <button
          onClick={() => setShowGrantModal(true)}
          className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
        >
          手动授予订阅
        </button>
      </div>
      
      {/* Filters */}
      <div className="flex gap-4 mb-6">
        <select
          value={statusFilter}
          onChange={(e) => {
            setStatusFilter(e.target.value)
            setPage(1)
          }}
          className="px-3 py-2 border rounded-md"
        >
          <option value="">全部状态</option>
          <option value="active">活跃</option>
          <option value="canceled">已取消</option>
          <option value="past_due">欠费</option>
        </select>
        <select
          value={planFilter}
          onChange={(e) => {
            setPlanFilter(e.target.value)
            setPage(1)
          }}
          className="px-3 py-2 border rounded-md"
        >
          <option value="">全部套餐</option>
          <option value="free">Free</option>
          <option value="pro">Pro</option>
          <option value="team">Team</option>
        </select>
      </div>
      
      {/* Table */}
      <div className="bg-white rounded-lg shadow overflow-hidden">
        <table className="w-full">
          <thead className="bg-gray-50">
            <tr className="text-left text-sm text-gray-500">
              <th className="px-6 py-3">用户ID</th>
              <th className="px-6 py-3">套餐</th>
              <th className="px-6 py-3">状态</th>
              <th className="px-6 py-3">开始时间</th>
              <th className="px-6 py-3">到期时间</th>
              <th className="px-6 py-3">来源</th>
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
            {data?.subscriptions?.map((sub: any) => (
              <tr key={sub.id} className="border-t hover:bg-gray-50">
                <td className="px-6 py-4 text-sm font-mono">{sub.userId?.slice(0, 12)}...</td>
                <td className="px-6 py-4">
                  <PlanBadge plan={sub.plan} />
                </td>
                <td className="px-6 py-4">
                  <StatusBadge status={sub.status} />
                </td>
                <td className="px-6 py-4 text-sm">
                  {sub.currentPeriodStart ? format(new Date(sub.currentPeriodStart), 'yyyy-MM-dd') : '-'}
                </td>
                <td className="px-6 py-4 text-sm">
                  {sub.currentPeriodEnd ? format(new Date(sub.currentPeriodEnd), 'yyyy-MM-dd') : '-'}
                </td>
                <td className="px-6 py-4 text-sm text-gray-500">
                  {sub.grantedBy ? '手动授予' : '订阅购买'}
                </td>
              </tr>
            ))}
            {data?.subscriptions?.length === 0 && (
              <tr>
                <td colSpan={6} className="px-6 py-8 text-center text-gray-500">
                  暂无订阅
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
      
      {/* Grant Modal */}
      {showGrantModal && (
        <GrantModal
          onClose={() => setShowGrantModal(false)}
          onSuccess={() => {
            setShowGrantModal(false)
            queryClient.invalidateQueries({ queryKey: ['subscriptions'] })
          }}
        />
      )}
    </div>
  )
}

function GrantModal({ onClose, onSuccess }: { onClose: () => void; onSuccess: () => void }) {
  const [userId, setUserId] = useState('')
  const [plan, setPlan] = useState('pro')
  const [months, setMonths] = useState(1)
  const [reason, setReason] = useState('')
  
  const mutation = useMutation({
    mutationFn: () => grantSubscription(userId, plan, months, reason),
    onSuccess,
  })
  
  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center">
      <div className="bg-white rounded-lg p-6 w-full max-w-md">
        <h2 className="text-lg font-semibold mb-4">手动授予订阅</h2>
        
        <div className="space-y-4">
          <div>
            <label className="block text-sm font-medium mb-1">用户ID</label>
            <input
              type="text"
              value={userId}
              onChange={(e) => setUserId(e.target.value)}
              className="w-full px-3 py-2 border rounded-md"
              placeholder="输入用户ID"
            />
          </div>
          
          <div>
            <label className="block text-sm font-medium mb-1">套餐</label>
            <select
              value={plan}
              onChange={(e) => setPlan(e.target.value)}
              className="w-full px-3 py-2 border rounded-md"
            >
              <option value="pro">Pro</option>
              <option value="team">Team</option>
            </select>
          </div>
          
          <div>
            <label className="block text-sm font-medium mb-1">时长（月）</label>
            <input
              type="number"
              value={months}
              onChange={(e) => setMonths(Number(e.target.value))}
              min={1}
              max={24}
              className="w-full px-3 py-2 border rounded-md"
            />
          </div>
          
          <div>
            <label className="block text-sm font-medium mb-1">原因</label>
            <textarea
              value={reason}
              onChange={(e) => setReason(e.target.value)}
              className="w-full px-3 py-2 border rounded-md"
              rows={2}
              placeholder="授予原因..."
            />
          </div>
        </div>
        
        <div className="flex gap-3 mt-6">
          <button
            onClick={onClose}
            className="flex-1 px-4 py-2 border rounded-md hover:bg-gray-50"
          >
            取消
          </button>
          <button
            onClick={() => mutation.mutate()}
            disabled={!userId || mutation.isPending}
            className="flex-1 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50"
          >
            {mutation.isPending ? '处理中...' : '确认授予'}
          </button>
        </div>
        
        {mutation.isError && (
          <p className="mt-3 text-sm text-red-600">
            {(mutation.error as Error).message}
          </p>
        )}
      </div>
    </div>
  )
}

function PlanBadge({ plan }: { plan: string }) {
  const config: Record<string, { label: string; className: string }> = {
    free: { label: 'Free', className: 'bg-gray-100 text-gray-800' },
    pro: { label: 'Pro', className: 'bg-blue-100 text-blue-800' },
    team: { label: 'Team', className: 'bg-purple-100 text-purple-800' },
  }
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config[plan]?.className || config.free.className}`}>
      {config[plan]?.label || plan}
    </span>
  )
}

function StatusBadge({ status }: { status: string }) {
  const config: Record<string, { label: string; className: string }> = {
    active: { label: '活跃', className: 'bg-green-100 text-green-800' },
    canceled: { label: '已取消', className: 'bg-gray-100 text-gray-800' },
    past_due: { label: '欠费', className: 'bg-red-100 text-red-800' },
    canceling: { label: '取消中', className: 'bg-yellow-100 text-yellow-800' },
  }
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config[status]?.className || 'bg-gray-100 text-gray-800'}`}>
      {config[status]?.label || status}
    </span>
  )
}
