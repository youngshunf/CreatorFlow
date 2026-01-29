import { useState, useEffect } from 'react'
import { getAvailableModels, type CloudModel, type RawCloudModel } from '@/api/models'
import type { ModelDefinition } from '@config/models'

/**
 * 将原始 API 模型转换为统一格式
 */
function rawModelToCloudModel(raw: RawCloudModel): CloudModel {
  return {
    id: raw.model_name,
    name: raw.display_name,
    description: raw.description,
    context_window: raw.max_context_length,
    provider: raw.provider_name,
    enabled: raw.enabled,
  }
}

/**
 * 将云端模型转换为本地模型定义格式
 */
function cloudModelToDefinition(model: CloudModel): ModelDefinition {
  // 从模型名称提取短名称 (如 "claude-sonnet-4-5" -> "sonnet-4-5")
  const shortName = model.name.replace('claude-', '').split('-').slice(0, 2).join('-') || model.name
  
  return {
    id: model.id,
    name: model.name,
    shortName,
    // 优先使用云端配置的描述，如果没有则留空
    description: model.description || '',
    contextWindow: model.context_window,
  }
}

/**
 * 从 API 响应中提取模型数组
 * 支持多种响应格式：
 * - 直接数组：[...]
 * - 对象包含数组：{ items: [...] } 或 { models: [...] } 等
 */
function extractModelsArray(response: unknown): CloudModel[] {
  // 直接是数组
  if (Array.isArray(response)) {
    return response.map(rawModelToCloudModel)
  }
  
  // 是对象，尝试从常见字段中提取
  if (response && typeof response === 'object') {
    const obj = response as Record<string, unknown>
    
    // 尝试各种可能的字段名
    for (const key of ['items', 'models', 'list', 'data']) {
      if (Array.isArray(obj[key])) {
        return (obj[key] as RawCloudModel[]).map(rawModelToCloudModel)
      }
    }
  }
  
  return []
}

/**
 * 获取云端可用模型列表的 Hook
 * 自动缓存结果，避免重复请求
 */
export function useCloudModels() {
  const [models, setModels] = useState<ModelDefinition[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)

  useEffect(() => {
    let mounted = true
    
    async function fetchModels() {
      try {
        setLoading(true)
        const response = await getAvailableModels()
        const cloudModels = extractModelsArray(response)
        
        if (mounted) {
          // 只显示启用的模型
          const enabledModels = cloudModels.filter(m => m.enabled !== false)
          
          if (enabledModels.length > 0) {
            setModels(enabledModels.map(cloudModelToDefinition))
            setError(null)
          }
          // 如果没有获取到模型，保持使用 fallback
        }
      } catch (err) {
        console.error('[useCloudModels] Failed to fetch models:', err)
        if (mounted) {
          setError(err instanceof Error ? err : new Error('Failed to fetch models'))
          // 获取失败时不清空已有模型，使用 fallback
        }
      } finally {
        if (mounted) {
          setLoading(false)
        }
      }
    }
    
    fetchModels()
    
    return () => {
      mounted = false
    }
  }, [])

  return { models, loading, error }
}
