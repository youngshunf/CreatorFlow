/**
 * 统一支付服务
 * 管理订单创建、查询、回调处理和订阅激活
 */

import { eq, and, desc } from 'drizzle-orm'
import { db } from '../db'
import { orders, subscriptions, type Order } from '../db/schema'

function getDb() {
  if (!db) throw new Error('Database not initialized')
  return db
}
import { getWechatPay } from '../lib/wechat-pay'
import { getAlipay } from '../lib/alipay'
import { getPlanPrice, getPlan } from '../config/plans'
import type { PaymentChannel, PaymentMethod, OrderStatus, PaymentResult, NotifyResult } from '../config/payment'

function generateOrderNo(): string {
  const date = new Date()
  const dateStr = date.toISOString().slice(0, 10).replace(/-/g, '')
  const random = Math.random().toString(36).substring(2, 10).toUpperCase()
  return `CF${dateStr}${random}`
}

function generateId(): string {
  return crypto.randomUUID()
}

export interface CreateOrderParams {
  userId: string
  plan: 'pro' | 'team'
  period: 'monthly' | 'yearly'
  channel: PaymentChannel
  method: PaymentMethod
  clientIp?: string  // H5支付需要
}

export interface CreateOrderResult {
  success: boolean
  orderId?: string
  orderNo?: string
  // 微信支付
  codeUrl?: string      // Native 扫码
  prepayId?: string     // JSAPI
  mwebUrl?: string      // H5 跳转
  // 支付宝
  formHtml?: string     // Page pay form
  redirectUrl?: string  // WAP redirect
  error?: string
}

/**
 * 创建支付订单
 */
export async function createOrder(params: CreateOrderParams): Promise<CreateOrderResult> {
  const { userId, plan, period, channel, method, clientIp } = params
  
  const planInfo = getPlan(plan)
  if (!planInfo || planInfo.priceCNY === 0) {
    return { success: false, error: 'Invalid plan' }
  }
  
  const amount = getPlanPrice(plan, period)
  const orderNo = generateOrderNo()
  const orderId = generateId()
  
  const subject = `CreatorFlow ${planInfo.name} - ${period === 'monthly' ? '月度订阅' : '年度订阅'}`
  const expireMinutes = 30
  const expiredAt = new Date(Date.now() + expireMinutes * 60 * 1000)
  
  // 创建订单记录
  await getDb().insert(orders).values({
    id: orderId,
    orderNo,
    userId,
    plan,
    period,
    amount,
    channel,
    method,
    status: 'pending',
    expiredAt,
  })
  
  // 调用支付渠道创建支付
  let paymentResult: PaymentResult
  
  if (channel === 'wechat') {
    const wechat = getWechatPay()
    switch (method) {
      case 'native':
        paymentResult = await wechat.createNativeOrder({ orderNo, amount, description: subject, expireMinutes })
        break
      case 'jsapi':
        // JSAPI 需要 openid，暂时返回错误
        return { success: false, error: 'JSAPI requires openid' }
      case 'h5':
        if (!clientIp) {
          return { success: false, error: 'H5 payment requires clientIp' }
        }
        paymentResult = await wechat.createH5Order({ orderNo, amount, description: subject, clientIp, expireMinutes })
        break
      default:
        return { success: false, error: 'Unsupported payment method' }
    }
  } else if (channel === 'alipay') {
    const alipay = getAlipay()
    switch (method) {
      case 'page':
        paymentResult = await alipay.createPagePay({ orderNo, amount, subject, expireMinutes })
        break
      case 'wap':
        paymentResult = await alipay.createWapPay({ orderNo, amount, subject, expireMinutes })
        break
      case 'qr':
        paymentResult = await alipay.createQrPay({ orderNo, amount, subject, expireMinutes })
        break
      default:
        return { success: false, error: 'Unsupported payment method' }
    }
  } else {
    return { success: false, error: 'Unsupported payment channel' }
  }
  
  if (!paymentResult.success) {
    // 支付创建失败，更新订单状态
    await getDb().update(orders)
      .set({ status: 'closed', updatedAt: new Date() })
      .where(eq(orders.id, orderId))
    
    return { success: false, error: paymentResult.error }
  }
  
  return {
    success: true,
    orderId,
    orderNo,
    codeUrl: paymentResult.codeUrl,
    prepayId: paymentResult.prepayId,
    mwebUrl: paymentResult.mwebUrl,
    formHtml: paymentResult.formHtml,
    redirectUrl: paymentResult.redirectUrl,
  }
}

/**
 * 查询订单状态
 */
export async function getOrderStatus(orderNo: string): Promise<{
  success: boolean
  order?: Order
  error?: string
}> {
  const [order] = await getDb().select()
    .from(orders)
    .where(eq(orders.orderNo, orderNo))
    .limit(1)
  
  if (!order) {
    return { success: false, error: 'Order not found' }
  }
  
  // 如果订单是 pending 状态，主动查询支付结果
  if (order.status === 'pending') {
    const queryResult = await queryPaymentStatus(order)
    if (queryResult.paid) {
      // 支付成功，处理回调
      await handlePaymentSuccess(order.orderNo, queryResult.transactionId!, queryResult.paidAmount!)
      const [updatedOrder] = await getDb().select()
        .from(orders)
        .where(eq(orders.orderNo, orderNo))
        .limit(1)
      return { success: true, order: updatedOrder }
    }
  }
  
  return { success: true, order }
}

/**
 * 主动查询支付状态
 */
async function queryPaymentStatus(order: Order): Promise<{
  paid: boolean
  transactionId?: string
  paidAmount?: number
}> {
  if (order.channel === 'wechat') {
    const wechat = getWechatPay()
    const result = await wechat.queryOrder(order.orderNo)
    if (result.success && result.status === 'SUCCESS') {
      return {
        paid: true,
        transactionId: result.transactionId,
        paidAmount: result.paidAmount,
      }
    }
  } else if (order.channel === 'alipay') {
    const alipay = getAlipay()
    const result = await alipay.queryOrder(order.orderNo)
    if (result.success && (result.status === 'TRADE_SUCCESS' || result.status === 'TRADE_FINISHED')) {
      return {
        paid: true,
        transactionId: result.transactionId,
        paidAmount: result.paidAmount,
      }
    }
  }
  
  return { paid: false }
}

/**
 * 处理微信支付回调
 */
export async function handleWechatNotify(body: string, headers: Record<string, string>): Promise<{
  success: boolean
  message?: string
}> {
  const wechat = getWechatPay()
  const result = wechat.verifyNotify(body, headers)
  
  if (!result.success) {
    console.error('Wechat notify verification failed:', result.error)
    return { success: false, message: result.error }
  }
  
  await handlePaymentSuccess(result.orderNo!, result.transactionId!, result.paidAmount!)
  return { success: true }
}

/**
 * 处理支付宝回调
 */
export async function handleAlipayNotify(params: Record<string, string>): Promise<{
  success: boolean
  message?: string
}> {
  const alipay = getAlipay()
  const result = alipay.verifyNotify(params)
  
  if (!result.success) {
    console.error('Alipay notify verification failed:', result.error)
    return { success: false, message: result.error }
  }
  
  await handlePaymentSuccess(result.orderNo!, result.transactionId!, result.paidAmount!)
  return { success: true }
}

/**
 * 支付成功后处理
 */
async function handlePaymentSuccess(orderNo: string, transactionId: string, paidAmount: number): Promise<void> {
  // 查找订单
  const [order] = await getDb().select()
    .from(orders)
    .where(eq(orders.orderNo, orderNo))
    .limit(1)
  
  if (!order) {
    console.error('Order not found:', orderNo)
    return
  }
  
  // 检查是否已处理
  if (order.status === 'paid') {
    console.log('Order already paid:', orderNo)
    return
  }
  
  // 验证金额
  if (paidAmount !== order.amount) {
    console.error('Amount mismatch:', { expected: order.amount, actual: paidAmount })
    // 仍然继续处理，但记录异常
  }
  
  const now = new Date()
  
  // 更新订单状态
  await getDb().update(orders)
    .set({
      status: 'paid',
      transactionId,
      paidAt: now,
      updatedAt: now,
    })
    .where(eq(orders.id, order.id))
  
  // 计算订阅周期
  const periodMonths = order.period === 'yearly' ? 12 : 1
  const periodEnd = new Date(now)
  periodEnd.setMonth(periodEnd.getMonth() + periodMonths)
  
  // 查找现有订阅
  const [existingSub] = await getDb().select()
    .from(subscriptions)
    .where(eq(subscriptions.userId, order.userId))
    .limit(1)
  
  if (existingSub) {
    // 更新现有订阅
    // 如果当前订阅有效，则延长；否则从现在开始
    const newPeriodStart = existingSub.status === 'active' && existingSub.currentPeriodEnd && existingSub.currentPeriodEnd > now
      ? existingSub.currentPeriodEnd
      : now
    const newPeriodEnd = new Date(newPeriodStart)
    newPeriodEnd.setMonth(newPeriodEnd.getMonth() + periodMonths)
    
    await getDb().update(subscriptions)
      .set({
        plan: order.plan,
        status: 'active',
        currentPeriodStart: newPeriodStart,
        currentPeriodEnd: newPeriodEnd,
        cancelAtPeriodEnd: false,
        updatedAt: now,
      })
      .where(eq(subscriptions.id, existingSub.id))
  } else {
    // 创建新订阅
    await getDb().insert(subscriptions).values({
      id: generateId(),
      userId: order.userId,
      plan: order.plan,
      status: 'active',
      currentPeriodStart: now,
      currentPeriodEnd: periodEnd,
    })
  }
  
  console.log(`Subscription activated for user ${order.userId}, plan: ${order.plan}, period: ${order.period}`)
}

/**
 * 获取用户订单历史
 */
export async function getUserOrders(userId: string, limit = 20): Promise<Order[]> {
  return getDb().select()
    .from(orders)
    .where(eq(orders.userId, userId))
    .orderBy(desc(orders.createdAt))
    .limit(limit)
}

/**
 * 关闭超时订单
 */
export async function closeExpiredOrders(): Promise<number> {
  const now = new Date()
  
  const expiredOrders = await getDb().select()
    .from(orders)
    .where(and(
      eq(orders.status, 'pending'),
    ))
  
  let closedCount = 0
  
  for (const order of expiredOrders) {
    if (order.expiredAt && order.expiredAt < now) {
      // 关闭第三方订单
      if (order.channel === 'wechat') {
        const wechat = getWechatPay()
        await wechat.closeOrder(order.orderNo)
      } else if (order.channel === 'alipay') {
        const alipay = getAlipay()
        await alipay.closeOrder(order.orderNo)
      }
      
      // 更新本地状态
      await getDb().update(orders)
        .set({ status: 'expired', updatedAt: now })
        .where(eq(orders.id, order.id))
      
      closedCount++
    }
  }
  
  return closedCount
}
