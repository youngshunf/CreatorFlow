import request from './request'

/**
 * 云端模型定义 (已转换为统一格式)
 */
export interface CloudModel {
  id: string
  name: string
  description?: string
  context_window?: number
  provider?: string
  enabled?: boolean
  visible?: boolean
}

/**
 * API 返回的原始模型格式
 * 兼容后端 GetAvailableModel (available 接口) 和 GetModelConfigList (分页列表接口)
 */
export interface RawCloudModel {
  id?: number
  model_id?: string        // available 接口返回的模型 ID
  model_name: string
  display_name: string
  description?: string
  model_type?: string
  max_context_length?: number
  max_tokens?: number
  provider?: string        // available 接口返回的供应商类型
  provider_name?: string   // 分页列表接口返回的供应商名称
  enabled?: boolean
  visible?: boolean
}

/**
 * 获取云端可用模型列表
 * 返回 unknown 类型，由调用方解析具体格式
 */
export async function getAvailableModels(): Promise<unknown> {
  return request.get('/llm/models/available')
}

/**
 * 用户 LLM 配置
 */
export interface UserLlmConfig {
  default_anthropic_model?: string
  default_openai_model?: string
}

/**
 * 获取用户的 LLM 配置（包括默认模型）
 */
export async function getUserLlmConfig(): Promise<UserLlmConfig> {
  return request.get('/llm/user/config')
}
