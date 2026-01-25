import { useState, useEffect } from 'react'
import { getDashboardStats } from '../lib/api'
import { format } from 'date-fns'

export function Dashboard() {
  const [data, setData] = useState<any>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  
  useEffect(() => {
    getDashboardStats()
      .then(setData)
      .catch(setError)
      .finally(() => setIsLoading(false))
  }, [])
  
  if (isLoading) {
    return <div className="animate-pulse">加载中...</div>
  }
  
  if (error) {
    return <div className="text-red-600">加载失败: {(error as Error).message}</div>
  }
  
  return (
    <div>
      <h1 className="text-2xl font-bold mb-8">仪表盘</h1>
      
      {/* Stats cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6 mb-8">
        <StatCard
          title="用户总数"
          value={data?.userCount ?? 0}
          icon="👥"
        />
        <StatCard
          title="活跃订阅"
          value={data?.activeSubscriptions ?? 0}
          icon="💳"
        />
        <StatCard
          title="月度收入"
          value={`¥${((data?.revenue?.mrr ?? 0) / 100).toFixed(2)}`}
          icon="💰"
        />
        <StatCard
          title="今日订单"
          value={data?.recentOrders?.length ?? 0}
          icon="📋"
        />
      </div>
      
      {/* Recent orders */}
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-lg font-semibold mb-4">最近订单</h2>
        <table className="w-full">
          <thead>
            <tr className="text-left text-sm text-gray-500 border-b">
              <th className="pb-2">订单号</th>
              <th className="pb-2">用户</th>
              <th className="pb-2">套餐</th>
              <th className="pb-2">金额</th>
              <th className="pb-2">状态</th>
              <th className="pb-2">时间</th>
            </tr>
          </thead>
          <tbody>
            {data?.recentOrders?.map((order: any) => (
              <tr key={order.orderNo} className="border-b last:border-0">
                <td className="py-3 text-sm font-mono">{order.orderNo}</td>
                <td className="py-3 text-sm">{order.userId?.slice(0, 8)}...</td>
                <td className="py-3 text-sm">{order.plan}</td>
                <td className="py-3 text-sm">¥{(order.amount / 100).toFixed(2)}</td>
                <td className="py-3">
                  <OrderStatus status={order.status} />
                </td>
                <td className="py-3 text-sm text-gray-500">
                  {format(new Date(order.createdAt), 'MM-dd HH:mm')}
                </td>
              </tr>
            ))}
            {(!data?.recentOrders || data.recentOrders.length === 0) && (
              <tr>
                <td colSpan={6} className="py-8 text-center text-gray-500">
                  暂无订单
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  )
}

function StatCard({ title, value, icon }: { title: string; value: string | number; icon: string }) {
  return (
    <div className="bg-white rounded-lg shadow p-6">
      <div className="flex items-center gap-4">
        <span className="text-3xl">{icon}</span>
        <div>
          <p className="text-sm text-gray-500">{title}</p>
          <p className="text-2xl font-bold">{value}</p>
        </div>
      </div>
    </div>
  )
}

function OrderStatus({ status }: { status: string }) {
  const statusConfig: Record<string, { label: string; className: string }> = {
    pending: { label: '待支付', className: 'bg-yellow-100 text-yellow-800' },
    paid: { label: '已支付', className: 'bg-green-100 text-green-800' },
    closed: { label: '已关闭', className: 'bg-gray-100 text-gray-800' },
    expired: { label: '已过期', className: 'bg-red-100 text-red-800' },
  }
  
  const config = statusConfig[status] || { label: status, className: 'bg-gray-100 text-gray-800' }
  
  return (
    <span className={`px-2 py-1 rounded text-xs ${config.className}`}>
      {config.label}
    </span>
  )
}
