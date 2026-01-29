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
}

/**
 * API 返回的原始模型格式
 */
export interface RawCloudModel {
  id?: number
  model_name: string
  display_name: string
  description?: string  // 模型描述，用于前端显示
  model_type?: string
  max_context_length?: number
  max_tokens?: number
  provider_name?: string
  enabled?: boolean
}

/**
 * 获取云端可用模型列表
 * 返回 unknown 类型，由调用方解析具体格式
 */
export async function getAvailableModels(): Promise<unknown> {
  return request.get('/llm/models')
}
