/**
 * 支付宝支付服务
 * 支持电脑网站支付、手机网站支付、当面付（扫码）
 * 文档: https://opendocs.alipay.com/open/repo-0038oa
 */

import crypto from 'crypto'
import { getPaymentConfig, type AlipayConfig, type PaymentResult, type NotifyResult } from '../config/payment'

const ALIPAY_GATEWAY = 'https://openapi.alipay.com/gateway.do'
const ALIPAY_SANDBOX_GATEWAY = 'https://openapi-sandbox.dl.alipaydev.com/gateway.do'

export class Alipay {
  private config: AlipayConfig
  private gateway: string
  
  constructor() {
    this.config = getPaymentConfig().alipay
    this.gateway = this.config.sandbox ? ALIPAY_SANDBOX_GATEWAY : ALIPAY_GATEWAY
  }
  
  /**
   * 电脑网站支付
   * 返回跳转表单HTML，前端直接渲染即可
   */
  async createPagePay(params: {
    orderNo: string
    amount: number       // 单位：分
    subject: string
    body?: string
    expireMinutes?: number
  }): Promise<PaymentResult> {
    const { orderNo, amount, subject, body, expireMinutes = 30 } = params
    
    const bizContent = {
      out_trade_no: orderNo,
      total_amount: (amount / 100).toFixed(2),  // 转换为元
      subject,
      body,
      product_code: 'FAST_INSTANT_TRADE_PAY',
      time_expire: this.getExpireTime(expireMinutes),
    }
    
    try {
      const formHtml = this.buildForm('alipay.trade.page.pay', bizContent)
      return {
        success: true,
        formHtml,
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      }
    }
  }
  
  /**
   * 手机网站支付
   * 返回跳转URL
   */
  async createWapPay(params: {
    orderNo: string
    amount: number
    subject: string
    body?: string
    quitUrl?: string
    expireMinutes?: number
  }): Promise<PaymentResult> {
    const { orderNo, amount, subject, body, quitUrl, expireMinutes = 30 } = params
    
    const bizContent = {
      out_trade_no: orderNo,
      total_amount: (amount / 100).toFixed(2),
      subject,
      body,
      product_code: 'QUICK_WAP_WAY',
      quit_url: quitUrl || this.config.returnUrl,
      time_expire: this.getExpireTime(expireMinutes),
    }
    
    try {
      const redirectUrl = this.buildUrl('alipay.trade.wap.pay', bizContent)
      return {
        success: true,
        redirectUrl,
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      }
    }
  }
  
  /**
   * 当面付（扫码支付）
   * 返回二维码链接
   */
  async createQrPay(params: {
    orderNo: string
    amount: number
    subject: string
    body?: string
    expireMinutes?: number
  }): Promise<PaymentResult> {
    const { orderNo, amount, subject, body, expireMinutes = 30 } = params
    
    const bizContent = {
      out_trade_no: orderNo,
      total_amount: (amount / 100).toFixed(2),
      subject,
      body,
      time_expire: this.getExpireTime(expireMinutes),
    }
    
    try {
      const response = await this.execute('alipay.trade.precreate', bizContent)
      
      if (response.alipay_trade_precreate_response?.code === '10000') {
        return {
          success: true,
          codeUrl: response.alipay_trade_precreate_response.qr_code,
        }
      }
      
      return {
        success: false,
        error: response.alipay_trade_precreate_response?.sub_msg || 'Failed to create QR payment',
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
    status?: 'WAIT_BUYER_PAY' | 'TRADE_CLOSED' | 'TRADE_SUCCESS' | 'TRADE_FINISHED'
    transactionId?: string
    paidAmount?: number
    error?: string
  }> {
    try {
      const response = await this.execute('alipay.trade.query', {
        out_trade_no: orderNo,
      })
      
      const data = response.alipay_trade_query_response
      
      if (data?.code === '10000') {
        return {
          success: true,
          status: data.trade_status,
          transactionId: data.trade_no,
          paidAmount: Math.round(parseFloat(data.total_amount) * 100),
        }
      }
      
      return {
        success: false,
        error: data?.sub_msg || 'Query failed',
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
      const response = await this.execute('alipay.trade.close', {
        out_trade_no: orderNo,
      })
      
      return response.alipay_trade_close_response?.code === '10000'
    } catch {
      return false
    }
  }
  
  /**
   * 验证回调通知
   */
  verifyNotify(params: Record<string, string>): NotifyResult {
    try {
      // 获取签名
      const sign = params.sign
      const signType = params.sign_type || 'RSA2'
      
      if (!sign) {
        return { success: false, error: 'Missing sign' }
      }
      
      // 验证签名
      const signStr = this.buildSignString(params)
      const verified = this.verifySign(signStr, sign, signType)
      
      if (!verified) {
        return { success: false, error: 'Invalid signature' }
      }
      
      // 检查交易状态
      const tradeStatus = params.trade_status
      
      if (tradeStatus === 'TRADE_SUCCESS' || tradeStatus === 'TRADE_FINISHED') {
        return {
          success: true,
          orderNo: params.out_trade_no,
          transactionId: params.trade_no,
          paidAmount: params.total_amount ? Math.round(parseFloat(params.total_amount) * 100) : 0,
        }
      }
      
      return {
        success: false,
        orderNo: params.out_trade_no,
        error: `Trade status: ${tradeStatus}`,
      }
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Verify failed',
      }
    }
  }
  
  /**
   * 执行API请求
   */
  private async execute(method: string, bizContent: object): Promise<any> {
    const params = this.buildParams(method, bizContent)
    
    const response = await fetch(this.gateway, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/x-www-form-urlencoded',
      },
      body: new URLSearchParams(params).toString(),
    })
    
    return response.json()
  }
  
  /**
   * 构建请求参数
   */
  private buildParams(method: string, bizContent: object): Record<string, string> {
    const params: Record<string, string> = {
      app_id: this.config.appId,
      method,
      format: 'JSON',
      charset: 'utf-8',
      sign_type: 'RSA2',
      timestamp: this.formatDate(new Date()),
      version: '1.0',
      notify_url: this.config.notifyUrl,
      biz_content: JSON.stringify(bizContent),
    }
    
    // 添加签名
    const signStr = this.buildSignString(params)
    params.sign = this.sign(signStr)
    
    return params
  }
  
  /**
   * 构建跳转表单
   */
  private buildForm(method: string, bizContent: object): string {
    const params = this.buildParams(method, bizContent)
    params.return_url = this.config.returnUrl
    
    const inputs = Object.entries(params)
      .map(([key, value]) => `<input type="hidden" name="${key}" value="${this.escapeHtml(value)}" />`)
      .join('\n')
    
    return `
      <form id="alipay-form" action="${this.gateway}" method="POST">
        ${inputs}
      </form>
      <script>document.getElementById('alipay-form').submit();</script>
    `.trim()
  }
  
  /**
   * 构建跳转URL
   */
  private buildUrl(method: string, bizContent: object): string {
    const params = this.buildParams(method, bizContent)
    params.return_url = this.config.returnUrl
    
    const query = new URLSearchParams(params).toString()
    return `${this.gateway}?${query}`
  }
  
  /**
   * 构建签名字符串
   */
  private buildSignString(params: Record<string, string>): string {
    const sorted = Object.keys(params)
      .filter(key => key !== 'sign' && key !== 'sign_type' && params[key])
      .sort()
      .map(key => `${key}=${params[key]}`)
      .join('&')
    
    return sorted
  }
  
  /**
   * 签名
   */
  private sign(content: string): string {
    const signer = crypto.createSign('RSA-SHA256')
    signer.update(content, 'utf8')
    return signer.sign(this.config.privateKey, 'base64')
  }
  
  /**
   * 验证签名
   */
  private verifySign(content: string, signature: string, signType: string): boolean {
    const verifier = crypto.createVerify(signType === 'RSA2' ? 'RSA-SHA256' : 'RSA-SHA1')
    verifier.update(content, 'utf8')
    return verifier.verify(this.config.alipayPublicKey, signature, 'base64')
  }
  
  /**
   * 获取过期时间
   */
  private getExpireTime(minutes: number): string {
    const date = new Date(Date.now() + minutes * 60 * 1000)
    return this.formatDate(date)
  }
  
  /**
   * 格式化日期
   */
  private formatDate(date: Date): string {
    const pad = (n: number) => n.toString().padStart(2, '0')
    return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())} ${pad(date.getHours())}:${pad(date.getMinutes())}:${pad(date.getSeconds())}`
  }
  
  /**
   * HTML转义
   */
  private escapeHtml(str: string): string {
    return str
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#x27;')
  }
}

// 单例
let instance: Alipay | null = null

export function getAlipay(): Alipay {
  if (!instance) {
    instance = new Alipay()
  }
  return instance
}
