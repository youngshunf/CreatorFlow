/**
 * 微信支付服务
 * 支持 Native（扫码）、H5、JSAPI 支付
 * 文档: https://pay.weixin.qq.com/wiki/doc/apiv3/apis/index.shtml
 */

import crypto from 'crypto'
import { getPaymentConfig, type WechatPayConfig, type PaymentResult, type NotifyResult } from '../config/payment'

const WECHAT_PAY_API = 'https://api.mch.weixin.qq.com'

export class WechatPay {
  private config: WechatPayConfig
  
  constructor() {
    this.config = getPaymentConfig().wechat
  }
  
  /**
   * Native支付 - 生成支付二维码
   * 适用于PC网站扫码支付
   */
  async createNativeOrder(params: {
    orderNo: string
    amount: number       // 单位：分
    description: string
    expireMinutes?: number
  }): Promise<PaymentResult> {
    const { orderNo, amount, description, expireMinutes = 30 } = params
    
    const expireTime = new Date(Date.now() + expireMinutes * 60 * 1000)
      .toISOString()
      .replace(/\.\d{3}Z$/, '+08:00')
    
    const body = {
      appid: this.config.appId,
      mchid: this.config.mchId,
      description,
      out_trade_no: orderNo,
      time_expire: expireTime,
      notify_url: this.config.notifyUrl,
      amount: {
        total: amount,
        currency: 'CNY',
      },
    }
    
    try {
      const response = await this.request('POST', '/v3/pay/transactions/native', body)
      
      if (response.code_url) {
        return {
          success: true,
          codeUrl: response.code_url,
        }
      }
      
      return {
        success: false,
        error: response.message || 'Failed to create payment',
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      }
    }
  }
  
  /**
   * H5支付 - 手机浏览器支付
   * 适用于手机网页
   */
  async createH5Order(params: {
    orderNo: string
    amount: number
    description: string
    clientIp: string
    expireMinutes?: number
  }): Promise<PaymentResult> {
    const { orderNo, amount, description, clientIp, expireMinutes = 30 } = params
    
    const expireTime = new Date(Date.now() + expireMinutes * 60 * 1000)
      .toISOString()
      .replace(/\.\d{3}Z$/, '+08:00')
    
    const body = {
      appid: this.config.appId,
      mchid: this.config.mchId,
      description,
      out_trade_no: orderNo,
      time_expire: expireTime,
      notify_url: this.config.notifyUrl,
      amount: {
        total: amount,
        currency: 'CNY',
      },
      scene_info: {
        payer_client_ip: clientIp,
        h5_info: {
          type: 'Wap',
        },
      },
    }
    
    try {
      const response = await this.request('POST', '/v3/pay/transactions/h5', body)
      
      if (response.h5_url) {
        return {
          success: true,
          mwebUrl: response.h5_url,
        }
      }
      
      return {
        success: false,
        error: response.message || 'Failed to create H5 payment',
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      }
    }
  }
  
  /**
   * 查询订单状态
   */
  async queryOrder(orderNo: string): Promise<{
    success: boolean
    status?: 'SUCCESS' | 'NOTPAY' | 'CLOSED' | 'REFUND' | 'PAYERROR'
    transactionId?: string
    paidAmount?: number
    error?: string
  }> {
    try {
      const response = await this.request(
        'GET',
        `/v3/pay/transactions/out-trade-no/${orderNo}?mchid=${this.config.mchId}`
      )
      
      return {
        success: true,
        status: response.trade_state,
        transactionId: response.transaction_id,
        paidAmount: response.amount?.total,
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      }
    }
  }
  
  /**
   * 关闭订单
   */
  async closeOrder(orderNo: string): Promise<boolean> {
    try {
      await this.request('POST', `/v3/pay/transactions/out-trade-no/${orderNo}/close`, {
        mchid: this.config.mchId,
      })
      return true
    } catch {
      return false
    }
  }
  
  /**
   * 验证回调通知
   * @param body 请求体 JSON 字符串
   * @param headers 微信签名相关 header
   */
  verifyNotify(body: string, headers: Record<string, string>): NotifyResult {
    try {
      // 验证签名
      const timestamp = headers['wechatpay-timestamp']
      const nonce = headers['wechatpay-nonce']
      const signature = headers['wechatpay-signature']
      
      if (!timestamp || !nonce || !signature) {
        return { success: false, error: 'Missing signature headers' }
      }
      
      // TODO: 实现完整的签名验证（需要微信平台证书）
      // 这里简化处理，生产环境需要完整验证
      
      // 解密通知内容
      const data = JSON.parse(body)
      const resource = data.resource
      
      if (!resource) {
        return { success: false, error: 'Missing resource in notify' }
      }
      
      const decrypted = this.decryptResource(resource)
      const payment = JSON.parse(decrypted)
      
      if (payment.trade_state === 'SUCCESS') {
        return {
          success: true,
          orderNo: payment.out_trade_no,
          transactionId: payment.transaction_id,
          paidAmount: payment.amount?.total,
        }
      }
      
      return {
        success: false,
        orderNo: payment.out_trade_no,
        error: `Trade state: ${payment.trade_state}`,
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Verify failed',
      }
    }
  }
  
  /**
   * 解密回调通知资源
   */
  private decryptResource(resource: {
    algorithm: string
    ciphertext: string
    associated_data: string
    nonce: string
  }): string {
    const { ciphertext, associated_data, nonce } = resource
    
    const ciphertextBuffer = Buffer.from(ciphertext, 'base64')
    const authTag = ciphertextBuffer.slice(-16)
    const data = ciphertextBuffer.slice(0, -16)
    
    const decipher = crypto.createDecipheriv(
      'aes-256-gcm',
      Buffer.from(this.config.apiV3Key),
      Buffer.from(nonce)
    )
    
    decipher.setAuthTag(authTag)
    decipher.setAAD(Buffer.from(associated_data))
    
    const decrypted = Buffer.concat([
      decipher.update(data),
      decipher.final(),
    ])
    
    return decrypted.toString('utf8')
  }
  
  /**
   * 发送API请求
   */
  private async request(method: string, path: string, body?: object): Promise<any> {
    const url = `${WECHAT_PAY_API}${path}`
    const timestamp = Math.floor(Date.now() / 1000).toString()
    const nonce = crypto.randomBytes(16).toString('hex')
    const bodyStr = body ? JSON.stringify(body) : ''
    
    // 构建签名串
    const signStr = `${method}\n${path}\n${timestamp}\n${nonce}\n${bodyStr}\n`
    
    // 使用私钥签名
    const sign = crypto.createSign('RSA-SHA256')
    sign.update(signStr)
    const signature = sign.sign(this.config.privateKey, 'base64')
    
    // 构建 Authorization
    const authorization = `WECHATPAY2-SHA256-RSA2048 mchid="${this.config.mchId}",nonce_str="${nonce}",timestamp="${timestamp}",serial_no="${this.config.serialNo}",signature="${signature}"`
    
    const response = await fetch(url, {
      method,
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': authorization,
      },
      body: method !== 'GET' ? bodyStr : undefined,
    })
    
    if (!response.ok) {
      const errorData = await response.json().catch(() => ({} as Record<string, unknown>))
      throw new Error((errorData as Record<string, string>).message || `HTTP ${response.status}`)
    }
    
    // 204 No Content
    if (response.status === 204) {
      return {}
    }
    
    return response.json()
  }
}

// 单例
let instance: WechatPay | null = null

export function getWechatPay(): WechatPay {
  if (!instance) {
    instance = new WechatPay()
  }
  return instance
}
