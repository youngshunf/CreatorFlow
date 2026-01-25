/**
 * 支付配置
 * 支持微信支付和支付宝
 */

export type PaymentChannel = 'wechat' | 'alipay'
// 微信: native(扫码), jsapi(公众号), h5(手机浏览器)
// 支付宝: page(PC网页), wap(手机网页), qr(扫码)
export type PaymentMethod = 'native' | 'jsapi' | 'h5' | 'app' | 'page' | 'wap' | 'qr'

export interface PaymentConfig {
  wechat: WechatPayConfig
  alipay: AlipayConfig
}

// 微信支付配置
export interface WechatPayConfig {
  appId: string           // 公众号/小程序 AppID
  mchId: string           // 商户号
  apiKey: string          // API密钥 (v2)
  apiV3Key: string        // APIv3密钥
  serialNo: string        // 证书序列号
  privateKey: string      // 商户私钥
  notifyUrl: string       // 支付回调地址
}

// 支付宝配置
export interface AlipayConfig {
  appId: string           // 应用ID
  privateKey: string      // 应用私钥
  alipayPublicKey: string // 支付宝公钥
  notifyUrl: string       // 支付回调地址
  returnUrl: string       // 支付成功跳转地址
  sandbox: boolean        // 是否沙箱环境
}

// 从环境变量加载配置
export function getPaymentConfig(): PaymentConfig {
  return {
    wechat: {
      appId: process.env.WECHAT_APP_ID || '',
      mchId: process.env.WECHAT_MCH_ID || '',
      apiKey: process.env.WECHAT_API_KEY || '',
      apiV3Key: process.env.WECHAT_API_V3_KEY || '',
      serialNo: process.env.WECHAT_SERIAL_NO || '',
      privateKey: process.env.WECHAT_PRIVATE_KEY || '',
      notifyUrl: `${process.env.APP_URL}/api/billing/notify/wechat`,
    },
    alipay: {
      appId: process.env.ALIPAY_APP_ID || '',
      privateKey: process.env.ALIPAY_PRIVATE_KEY || '',
      alipayPublicKey: process.env.ALIPAY_PUBLIC_KEY || '',
      notifyUrl: `${process.env.APP_URL}/api/billing/notify/alipay`,
      returnUrl: `${process.env.APP_URL}/billing/success`,
      sandbox: process.env.ALIPAY_SANDBOX === 'true',
    },
  }
}

// 订单状态
export type OrderStatus = 'pending' | 'paid' | 'failed' | 'refunded' | 'closed'

// 订单类型
export interface PaymentOrder {
  id: string
  userId: string
  orderNo: string          // 商户订单号
  transactionId?: string   // 第三方交易号
  channel: PaymentChannel
  method: PaymentMethod
  plan: string             // 订阅计划
  amount: number           // 金额（分）
  status: OrderStatus
  paidAt?: Date
  expireAt: Date           // 订单过期时间
  metadata?: Record<string, any>
  createdAt: Date
  updatedAt: Date
}

// 支付结果
export interface PaymentResult {
  success: boolean
  // Native/QR支付返回二维码链接
  codeUrl?: string
  // H5支付返回跳转链接
  mwebUrl?: string
  // JSAPI支付返回预支付ID
  prepayId?: string
  // 支付宝返回表单
  formHtml?: string
  // 支付宝返回跳转URL
  redirectUrl?: string
  // 错误信息
  error?: string
}

// 回调验证结果
export interface NotifyResult {
  success: boolean
  orderNo?: string
  transactionId?: string
  paidAmount?: number
  error?: string
}
