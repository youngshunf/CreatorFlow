import crypto from 'crypto'

// In-memory store for verification codes (use Redis in production)
const verificationCodes = new Map<string, { code: string; expiresAt: number }>()

// Code expiration time (5 minutes)
const CODE_EXPIRATION_MS = 5 * 60 * 1000

// Cooldown between sending codes (60 seconds)
const SEND_COOLDOWN_MS = 60 * 1000

// Track last send time
const lastSendTime = new Map<string, number>()

// Aliyun SMS configuration
const ALIYUN_SMS_CONFIG = {
  accessKeyId: process.env.ALIYUN_ACCESS_KEY_ID || '',
  accessKeySecret: process.env.ALIYUN_ACCESS_KEY_SECRET || '',
  signName: process.env.ALIYUN_SMS_SIGN_NAME || '',
  templateCode: process.env.ALIYUN_SMS_TEMPLATE_CODE || '',
  endpoint: 'https://dysmsapi.aliyuncs.com',
}

export class SmsService {
  /**
   * Send verification code to phone number
   */
  static async sendCode(phone: string): Promise<{ success: boolean; message?: string }> {
    // Check cooldown
    const lastSent = lastSendTime.get(phone)
    if (lastSent && Date.now() - lastSent < SEND_COOLDOWN_MS) {
      const remaining = Math.ceil((SEND_COOLDOWN_MS - (Date.now() - lastSent)) / 1000)
      return { 
        success: false, 
        message: `请等待 ${remaining} 秒后再发送` 
      }
    }
    
    // Generate 6-digit code
    const code = Math.floor(100000 + Math.random() * 900000).toString()
    
    // Store code with expiration
    verificationCodes.set(phone, {
      code,
      expiresAt: Date.now() + CODE_EXPIRATION_MS,
    })
    
    // Update last send time
    lastSendTime.set(phone, Date.now())
    
    // Development: log the code
    if (process.env.NODE_ENV !== 'production') {
      console.log(`[SMS] Verification code for ${phone}: ${code}`)
      return { success: true }
    }
    
    // Production: send via Aliyun SMS
    try {
      await this.sendAliyunSms(phone, code)
      return { success: true }
    } catch (error) {
      console.error('[SMS] Failed to send:', error)
      // Remove stored code on failure
      verificationCodes.delete(phone)
      lastSendTime.delete(phone)
      return { success: false, message: '短信发送失败，请稍后重试' }
    }
  }
  
  /**
   * Verify the code for a phone number
   */
  static async verifyCode(phone: string, code: string): Promise<boolean> {
    const stored = verificationCodes.get(phone)
    
    if (!stored) {
      return false
    }
    
    // Check expiration
    if (Date.now() > stored.expiresAt) {
      verificationCodes.delete(phone)
      return false
    }
    
    // Check code match
    if (stored.code !== code) {
      return false
    }
    
    // Code is valid - remove it (one-time use)
    verificationCodes.delete(phone)
    return true
  }
  
  /**
   * Send SMS via Aliyun SMS API
   * Docs: https://help.aliyun.com/document_detail/419273.html
   */
  private static async sendAliyunSms(phone: string, code: string): Promise<void> {
    const { accessKeyId, accessKeySecret, signName, templateCode, endpoint } = ALIYUN_SMS_CONFIG
    
    if (!accessKeyId || !accessKeySecret || !signName || !templateCode) {
      throw new Error('Aliyun SMS not configured')
    }
    
    // Build request parameters
    const params: Record<string, string> = {
      AccessKeyId: accessKeyId,
      Action: 'SendSms',
      Format: 'JSON',
      PhoneNumbers: phone,
      SignName: signName,
      SignatureMethod: 'HMAC-SHA1',
      SignatureNonce: crypto.randomUUID(),
      SignatureVersion: '1.0',
      TemplateCode: templateCode,
      TemplateParam: JSON.stringify({ code }),
      Timestamp: new Date().toISOString().replace(/\.\d{3}Z$/, 'Z'),
      Version: '2017-05-25',
    }
    
    // Sort and encode parameters
    const sortedKeys = Object.keys(params).sort()
    const canonicalizedQueryString = sortedKeys
      .map(key => `${encodeURIComponent(key)}=${encodeURIComponent(params[key]!)}`)
      .join('&')
    
    // Build string to sign
    const stringToSign = `GET&${encodeURIComponent('/')}&${encodeURIComponent(canonicalizedQueryString)}`
    
    // Calculate signature
    const signature = crypto
      .createHmac('sha1', accessKeySecret + '&')
      .update(stringToSign)
      .digest('base64')
    
    // Build final URL
    const url = `${endpoint}/?${canonicalizedQueryString}&Signature=${encodeURIComponent(signature)}`
    
    // Send request
    const response = await fetch(url)
    const result = await response.json() as { Code?: string; Message?: string }
    
    if (result.Code !== 'OK') {
      throw new Error(`Aliyun SMS error: ${result.Message || result.Code}`)
    }
  }
  
  /**
   * Clean up expired codes (call periodically)
   */
  static cleanup(): void {
    const now = Date.now()
    for (const [phone, data] of verificationCodes.entries()) {
      if (now > data.expiresAt) {
        verificationCodes.delete(phone)
      }
    }
  }
}

// Cleanup expired codes every 5 minutes
setInterval(() => SmsService.cleanup(), 5 * 60 * 1000)
