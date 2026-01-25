import { useState, useEffect } from 'react'
import { getOrders } from '../lib/api'
import { format } from 'date-fns'

export function Orders() {
  const [page, setPage] = useState(1)
  const [statusFilter, setStatusFilter] = useState('')
  const [channelFilter, setChannelFilter] = useState('')
  const [data, setData] = useState<any>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  const limit = 20
  
  useEffect(() => {
    setIsLoading(true)
    getOrders({
      page,
      limit,
      status: statusFilter || undefined,
      channel: channelFilter || undefined,
    })
      .then(setData)
      .catch(setError)
      .finally(() => setIsLoading(false))
  }, [page, statusFilter, channelFilter])
  
  const totalPages = Math.ceil((data?.total ?? 0) / limit)
  
  return (
    <div>
      <h1 className="text-2xl font-bold mb-8">订单管理</h1>
      
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
          <option value="pending">待支付</option>
          <option value="paid">已支付</option>
          <option value="closed">已关闭</option>
          <option value="expired">已过期</option>
        </select>
        <select
          value={channelFilter}
          onChange={(e) => {
            setChannelFilter(e.target.value)
            setPage(1)
          }}
          className="px-3 py-2 border rounded-md"
        >
          <option value="">全部渠道</option>
          <option value="wechat">微信支付</option>
          <option value="alipay">支付宝</option>
        </select>
      </div>
      
      {/* Table */}
      <div className="bg-white rounded-lg shadow overflow-hidden">
        <table className="w-full">
          <thead className="bg-gray-50">
            <tr className="text-left text-sm text-gray-500">
              <th className="px-6 py-3">订单号</th>
              <th className="px-6 py-3">用户ID</th>
              <th className="px-6 py-3">套餐</th>
              <th className="px-6 py-3">周期</th>
              <th className="px-6 py-3">金额</th>
              <th className="px-6 py-3">渠道</th>
              <th className="px-6 py-3">状态</th>
              <th className="px-6 py-3">创建时间</th>
            </tr>
          </thead>
          <tbody>
            {isLoading && (
              <tr>
                <td colSpan={8} className="px-6 py-8 text-center text-gray-500">
                  加载中...
                </td>
              </tr>
            )}
            {error && (
              <tr>
                <td colSpan={8} className="px-6 py-8 text-center text-red-600">
                  加载失败: {(error as Error).message}
                </td>
              </tr>
            )}
            {data?.orders?.map((order: any) => (
              <tr key={order.orderNo} className="border-t hover:bg-gray-50">
                <td className="px-6 py-4 text-sm font-mono">{order.orderNo}</td>
                <td className="px-6 py-4 text-sm font-mono">{order.userId?.slice(0, 12)}...</td>
                <td className="px-6 py-4">
                  <PlanBadge plan={order.plan} />
                </td>
                <td className="px-6 py-4 text-sm">
                  {order.period === 'monthly' ? '月付' : '年付'}
                </td>
                <td className="px-6 py-4 text-sm font-medium">
                  ¥{(order.amount / 100).toFixed(2)}
                </td>
                <td className="px-6 py-4">
                  <ChannelBadge channel={order.channel} />
                </td>
                <td className="px-6 py-4">
                  <StatusBadge status={order.status} />
                </td>
                <td className="px-6 py-4 text-sm text-gray-500">
                  {format(new Date(order.createdAt), 'MM-dd HH:mm')}
                </td>
              </tr>
            ))}
            {data?.orders?.length === 0 && (
              <tr>
                <td colSpan={8} className="px-6 py-8 text-center text-gray-500">
                  暂无订单
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
  const config: Record<string, { label: string; className: string }> = {
    pro: { label: 'Pro', className: 'bg-blue-100 text-blue-800' },
    team: { label: 'Team', className: 'bg-purple-100 text-purple-800' },
  }
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config[plan]?.className || 'bg-gray-100 text-gray-800'}`}>
      {config[plan]?.label || plan}
    </span>
  )
}

function ChannelBadge({ channel }: { channel: string }) {
  const config: Record<string, { label: string; className: string }> = {
    wechat: { label: '微信', className: 'bg-green-100 text-green-800' },
    alipay: { label: '支付宝', className: 'bg-blue-100 text-blue-800' },
  }
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config[channel]?.className || 'bg-gray-100 text-gray-800'}`}>
      {config[channel]?.label || channel}
    </span>
  )
}

function StatusBadge({ status }: { status: string }) {
  const config: Record<string, { label: string; className: string }> = {
    pending: { label: '待支付', className: 'bg-yellow-100 text-yellow-800' },
    paid: { label: '已支付', className: 'bg-green-100 text-green-800' },
    closed: { label: '已关闭', className: 'bg-gray-100 text-gray-800' },
    expired: { label: '已过期', className: 'bg-red-100 text-red-800' },
  }
  return (
    <span className={`px-2 py-1 rounded text-xs font-medium ${config[status]?.className || 'bg-gray-100 text-gray-800'}`}>
      {config[status]?.label || status}
    </span>
  )
}
