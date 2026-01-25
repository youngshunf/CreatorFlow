import { Hono } from 'hono'
import { zValidator } from '@hono/zod-validator'
import { z } from 'zod'
import { authMiddleware } from '../middleware/auth'
import { SubscriptionService } from '../services/subscription'
import { PLANS } from '../config/plans'
import {
  createOrder,
  getOrderStatus,
  getUserOrders,
  handleWechatNotify,
  handleAlipayNotify,
} from '../services/payment'
import type { PaymentChannel, PaymentMethod } from '../config/payment'
import type { Env } from '../app'

const billing = new Hono<Env>()

// Get current subscription status
billing.get('/subscription', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const subscription = await SubscriptionService.getByUserId(userId)
  const usage = await SubscriptionService.getUsage(userId)
  
  const plan = subscription?.plan || 'free'
  const limits = PLANS[plan as keyof typeof PLANS]?.limits || PLANS.free!.limits
  
  return c.json({
    subscription: subscription ? {
      id: subscription.id,
      plan: subscription.plan,
      status: subscription.status,
      currentPeriodEnd: subscription.currentPeriodEnd,
    } : null,
    plan,
    limits,
    usage: {
      messagesToday: usage.messagesToday,
      tokensThisMonth: usage.tokensThisMonth,
    },
  })
})

// Get available plans with CNY prices
billing.get('/plans', async (c) => {
  return c.json({
    plans: Object.entries(PLANS).map(([key, plan]) => ({
      id: key,
      name: plan.name,
      priceUSD: plan.price,
      priceCNY: plan.priceCNY / 100,              // 转换为元
      yearlyPriceCNY: plan.yearlyPriceCNY / 100,  // 转换为元
      limits: plan.limits,
    })),
  })
})

// Create payment order
billing.post(
  '/order',
  authMiddleware,
  zValidator('json', z.object({
    plan: z.enum(['pro', 'team']),
    period: z.enum(['monthly', 'yearly']),
    channel: z.enum(['wechat', 'alipay']),
    method: z.enum(['native', 'jsapi', 'h5', 'page', 'wap', 'qr']),
  })),
  async (c) => {
    const userId = c.get('userId')
    const { plan, period, channel, method } = c.req.valid('json')
    
    const result = await createOrder({
      userId,
      plan,
      period,
      channel: channel as PaymentChannel,
      method: method as PaymentMethod,
    })
    
    if (!result.success) {
      return c.json({ error: result.error }, 400)
    }
    
    return c.json({
      orderId: result.orderId,
      orderNo: result.orderNo,
      // 微信支付
      codeUrl: result.codeUrl,
      prepayId: result.prepayId,
      mwebUrl: result.mwebUrl,
      // 支付宝
      formHtml: result.formHtml,
      redirectUrl: result.redirectUrl,
    })
  }
)

// Query order status
billing.get('/order/:orderNo', authMiddleware, async (c) => {
  const orderNo = c.req.param('orderNo')
  const userId = c.get('userId')
  
  const result = await getOrderStatus(orderNo)
  
  if (!result.success) {
    return c.json({ error: result.error }, 404)
  }
  
  // 验证订单属于当前用户
  if (result.order!.userId !== userId) {
    return c.json({ error: 'Order not found' }, 404)
  }
  
  return c.json({
    orderNo: result.order!.orderNo,
    status: result.order!.status,
    plan: result.order!.plan,
    period: result.order!.period,
    amount: result.order!.amount / 100,  // 转换为元
    paidAt: result.order!.paidAt,
    createdAt: result.order!.createdAt,
  })
})

// Get user order history
billing.get('/orders', authMiddleware, async (c) => {
  const userId = c.get('userId')
  const orders = await getUserOrders(userId)
  
  return c.json({
    orders: orders.map(order => ({
      orderNo: order.orderNo,
      status: order.status,
      plan: order.plan,
      period: order.period,
      amount: order.amount / 100,
      channel: order.channel,
      paidAt: order.paidAt,
      createdAt: order.createdAt,
    })),
  })
})

// WeChat pay callback
billing.post('/notify/wechat', async (c) => {
  const body = await c.req.text()
  const headers: Record<string, string> = {}
  
  // 提取微信签名相关 header
  const signatureHeaders = [
    'Wechatpay-Timestamp',
    'Wechatpay-Nonce',
    'Wechatpay-Signature',
    'Wechatpay-Serial',
  ]
  for (const key of signatureHeaders) {
    const value = c.req.header(key)
    if (value) headers[key] = value
  }
  
  const result = await handleWechatNotify(body, headers)
  
  if (!result.success) {
    return c.json({ code: 'FAIL', message: result.message }, 400)
  }
  
  return c.json({ code: 'SUCCESS', message: 'OK' })
})

// Alipay callback
billing.post('/notify/alipay', async (c) => {
  const formData = await c.req.parseBody()
  const params: Record<string, string> = {}
  
  for (const [key, value] of Object.entries(formData)) {
    if (typeof value === 'string') {
      params[key] = value
    }
  }
  
  const result = await handleAlipayNotify(params)
  
  if (!result.success) {
    return c.text('fail')
  }
  
  return c.text('success')
})

export { billing }
