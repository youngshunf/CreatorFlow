// 使用原生 fetch 替代 axios，避免依赖安装问题

import { getCloudApiUrl, getCurrentEnv, isDebugMode } from '@config/environments'

// 定义后端响应的标准结构
export interface ApiResponse<T = any> {
  code: number;
  msg: string;
  data: T;
}

// 从 environments.ts 获取 API URL（根据构建环境自动选择）
const BASE_URL = getCloudApiUrl();

// Debug: Log BASE_URL in staging/development
if (isDebugMode()) {
  console.log('[Request] BASE_URL:', BASE_URL);
  console.log('[Request] Environment:', getCurrentEnv());
}

interface RequestOptions extends RequestInit {
  params?: Record<string, string>;
}

class RequestClient {
  private baseURL: string;

  constructor(baseURL: string) {
    this.baseURL = baseURL;
  }

  private async request<T>(url: string, options: RequestOptions = {}): Promise<T> {
    const { params, ...init } = options;
    
    // 构建 URL
    let fullUrl = url.startsWith('http') ? url : `${this.baseURL}${url}`;
    if (params) {
      const searchParams = new URLSearchParams(params);
      fullUrl += `?${searchParams.toString()}`;
    }

    // 设置 Headers
    const headers = new Headers(init.headers);
    if (!headers.has('Content-Type') && !(init.body instanceof FormData)) {
      headers.set('Content-Type', 'application/json');
    }

    // 注入 Token
    const token = localStorage.getItem('access_token');
    if (token) {
      headers.set('Authorization', `Bearer ${token}`);
    }

    // Debug logging
    if (isDebugMode()) {
      console.log('[Request] Fetching:', fullUrl);
      console.log('[Request] Method:', init.method || 'GET');
    }

    try {
      const startTime = Date.now();
      const response = await fetch(fullUrl, {
        ...init,
        headers,
      });
      
      if (isDebugMode()) {
        console.log(`[Request] Response: ${response.status} in ${Date.now() - startTime}ms`);
      }

      // 处理 401 未授权
      if (response.status === 401) {
        localStorage.removeItem('access_token');
        localStorage.removeItem('llm_token');
        window.dispatchEvent(new Event('auth:unauthorized'));
        throw new Error('Unauthorized');
      }

      // 处理 HTTP 错误
      if (!response.ok) {
        const errorText = await response.text();
        let errorMessage = response.statusText;
        try {
          const errorJson = JSON.parse(errorText);
          errorMessage = errorJson.msg || errorJson.error || errorMessage;
        } catch {
          // JSON 解析失败，使用原始文本或状态文本
          if (errorText) errorMessage = errorText;
        }
        throw new Error(errorMessage);
      }

      // 解析 JSON
      const res = await response.json() as ApiResponse<T>;
      
      // 处理业务错误 (clound-backend 标准: code 200 为成功)
      if (res.code !== undefined && res.code !== 200) {
        throw new Error(res.msg || 'API Error');
      }

      // 如果后端返回标准结构，取 data；否则直接返回整个 body (兼容性)
      return (res.code === 200 && res.data !== undefined) ? res.data : res;

    } catch (error: any) {
      console.error('[Request] Error:', error.message || error);
      console.error('[Request] URL was:', fullUrl);
      throw error;
    }
  }

  get<T>(url: string, options?: RequestOptions) {
    return this.request<T>(url, { ...options, method: 'GET' });
  }

  post<T>(url: string, data?: any, options?: RequestOptions) {
    return this.request<T>(url, {
      ...options,
      method: 'POST',
      body: JSON.stringify(data),
    });
  }

  put<T>(url: string, data?: any, options?: RequestOptions) {
    return this.request<T>(url, {
      ...options,
      method: 'PUT',
      body: JSON.stringify(data),
    });
  }

  delete<T>(url: string, options?: RequestOptions) {
    return this.request<T>(url, { ...options, method: 'DELETE' });
  }

  upload<T>(url: string, formData: FormData, options?: RequestOptions) {
    return this.request<T>(url, {
      ...options,
      method: 'POST',
      body: formData,
    });
  }
}

const request = new RequestClient(BASE_URL);

export default request;
