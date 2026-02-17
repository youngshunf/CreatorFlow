# 媒体生成网关设计方案

> 日期: 2026-02-17
> 状态: 设计完成，待实施
> 位置: clound-backend/backend/app/llm/

## 概述

在现有 LiteLLM 文本网关基础上，扩展图像和视频生成能力。采用 OpenAI 兼容 API 风格，通过适配器模式统一各厂商接口差异。

**优先接入厂商**：DALL-E 3（图像）、通义万相（图像）、可灵 Kling（视频）

**核心决策**：
- API 风格：OpenAI 兼容（/v1/images/generations、/v1/videos/generations）
- 异步策略：轮询 + 可选 Webhook 回调
- 存储策略：生成结果自动转存到自有 OSS，返回永久链接

---

## 1. 目录结构

在现有 clound-backend/backend/app/llm/ 下扩展：

    backend/app/llm/
    ├── core/
    │   ├── gateway.py                  # 现有文本网关（不动）
    │   ├── media_gateway.py            # 新增：媒体生成网关
    │   └── media_adapters/             # 新增：厂商适配器
    │       ├── __init__.py
    │       ├── base.py                 # 抽象基类
    │       ├── registry.py             # 适配器注册表
    │       ├── openai_image.py         # DALL-E 3
    │       ├── dashscope_image.py      # 通义万相图像
    │       └── kling_video.py          # 可灵视频
    ├── api/v1/
    │   ├── proxy.py                    # 现有（不动）
    │   ├── images.py                   # 新增：图像生成 API
    │   └── videos.py                   # 新增：视频生成 API
    ├── model/
    │   ├── media_task.py               # 新增：异步任务表
    │   └── ...
    ├── schema/
    │   ├── image.py                    # 新增：图像请求/响应 Schema
    │   └── video.py                    # 新增：视频请求/响应 Schema
    └── enums.py                        # 扩展：TaskStatus、MediaErrorCode

---

## 2. 适配器抽象

### 2.1 基类

    class MediaAdapter(ABC):
        provider_type: str           # 对应 ProviderType 枚举
        media_type: str              # "image" | "video"

        @abstractmethod
        async def submit(self, request: MediaRequest) -> SubmitResult:
            """提交生成任务，返回厂商任务ID或直接结果"""

        @abstractmethod
        async def poll(self, vendor_task_id: str) -> PollResult:
            """查询任务状态，返回进度或结果URL"""

        def normalize_params(self, request: MediaRequest) -> dict:
            """将统一参数转换为厂商特定参数（可覆盖）"""

### 2.2 注册机制

    ADAPTER_REGISTRY: dict[str, type[MediaAdapter]] = {}

    def register_adapter(provider: str, media_type: str):
        def decorator(cls):
            ADAPTER_REGISTRY[f"{provider}:{media_type}"] = cls
            return cls
        return decorator

    def get_adapter(provider: str, media_type: str, config: ModelConfig) -> MediaAdapter:
        key = f"{provider}:{media_type}"
        cls = ADAPTER_REGISTRY.get(key)
        if not cls:
            raise UnsupportedModelError(f"不支持的适配器: {key}")
        return cls(config)

### 2.3 厂商差异对照

    | 特性       | DALL-E 3       | 通义万相           | 可灵视频          |
    |-----------|----------------|-------------------|------------------|
    | 接口风格   | 同步返回        | 异步轮询           | 异步轮询          |
    | submit    | 直接拿到图片URL  | 返回 task_id       | 返回 task_id      |
    | poll      | 不需要          | 需要               | 需要              |
    | 认证方式   | Bearer token    | API-Key header     | AccessKey + 签名  |
    | 尺寸参数   | "1024x1024"     | {"width":1024,...} | aspect_ratio      |

---

## 3. API 接口

### 3.1 图像生成 — POST /v1/images/generations

请求：

    {
      "model": "dall-e-3",
      "prompt": "一只在月球上弹吉他的猫",
      "n": 1,
      "size": "1024x1024",
      "quality": "hd",
      "response_format": "url",
      "webhook_url": null
    }

响应（同步返回，超时 60s 降级为异步）：

    {
      "id": "img-abc123",
      "status": "completed",
      "data": [
        {
          "url": "https://your-oss.com/media/img-abc123-0.png",
          "revised_prompt": "..."
        }
      ],
      "usage": { "credits": 0.04 },
      "created": 1708123456
    }

### 3.2 视频生成 — POST /v1/videos/generations

请求：

    {
      "model": "kling-v1",
      "prompt": "航拍城市日落延时摄影",
      "duration": 5,
      "aspect_ratio": "16:9",
      "image_url": null,
      "webhook_url": "https://your-app.com/hooks/video-done"
    }

响应（立即返回任务ID）：

    {
      "id": "vid-xyz789",
      "status": "processing",
      "progress": 0,
      "estimated_seconds": 120,
      "created": 1708123456
    }

### 3.3 视频轮询 — GET /v1/videos/generations/{task_id}

    // 完成
    {
      "id": "vid-xyz789",
      "status": "completed",
      "data": {
        "url": "https://your-oss.com/media/vid-xyz789.mp4",
        "duration": 5.0,
        "resolution": "1920x1080"
      },
      "usage": { "credits": 1.5 },
      "created": 1708123456
    }

    // 失败
    {
      "id": "vid-xyz789",
      "status": "failed",
      "error": { "code": "content_policy", "message": "..." }
    }

### 3.4 任务状态

统一四种：pending → processing → completed / failed

---

## 4. 异步任务管理

### 4.1 数据库任务表

    class MediaTask(Base):
        __tablename__ = 'media_task'

        id: Mapped[str]                  # "img-xxx" / "vid-xxx" 前缀区分
        user_id: Mapped[int]
        model_name: Mapped[str]
        media_type: Mapped[str]          # "image" / "video"
        status: Mapped[str]              # pending / processing / completed / failed
        progress: Mapped[int]            # 0-100
        prompt: Mapped[str]
        params: Mapped[dict]             # JSONB，原始请求参数
        vendor_task_id: Mapped[str|None]
        vendor_urls: Mapped[list|None]   # JSONB，厂商临时URL
        oss_urls: Mapped[list|None]      # JSONB，转存后永久URL
        error_code: Mapped[str|None]
        error_message: Mapped[str|None]
        webhook_url: Mapped[str|None]
        credits_cost: Mapped[Decimal]
        created_at: Mapped[datetime]
        completed_at: Mapped[datetime|None]

### 4.2 处理流水线

    用户请求
      ├─ 1. 创建 media_task 记录（status=pending）
      ├─ 2. 预扣积分
      ├─ 3. 调用适配器 submit()
      │     ├─ 同步型（DALL-E）→ 直接拿到 vendor_urls
      │     └─ 异步型（万相/可灵）→ 拿到 vendor_task_id
      ├─ 4. [异步型] Celery 定时轮询
      │     └─ 调用适配器 poll() → 更新 progress / vendor_urls
      ├─ 5. 拿到 vendor_urls → 下载并上传到 OSS
      │     └─ 更新 oss_urls，status=completed
      ├─ 6. 结算积分（预扣 vs 实际差额退回）
      └─ 7. [可选] 触发 webhook 回调

### 4.3 轮询策略

退避间隔：前 30 秒每 2 秒查一次，之后每 5 秒，超过 10 分钟标记超时失败。

### 4.4 OSS 转存

独立 Celery 任务，与生成逻辑解耦。OSS 上传失败时任务仍为 completed，vendor_urls 作为降级。

---

## 5. 参数映射

### 5.1 DALL-E 3

    @register_adapter("openai", "image")
    class OpenAIImageAdapter(MediaAdapter):
        def normalize_params(self, req):
            return {
                "prompt": req.prompt,
                "model": "dall-e-3",
                "n": req.n,
                "size": req.size,
                "quality": req.quality or "standard",
                "response_format": "url",
            }

### 5.2 通义万相

    @register_adapter("qwen", "image")
    class DashScopeImageAdapter(MediaAdapter):
        def normalize_params(self, req):
            w, h = req.size.split("x")
            return {
                "model": "wanx-v1",
                "input": {"prompt": req.prompt},
                "parameters": {
                    "size": f"{w}*{h}",
                    "n": req.n,
                    "style": req.style or "<auto>",
                },
            }

### 5.3 可灵视频

    @register_adapter("kling", "video")
    class KlingVideoAdapter(MediaAdapter):
        def normalize_params(self, req):
            return {
                "model_name": "kling-v1",
                "prompt": req.prompt,
                "cfg_scale": req.cfg_scale or 0.5,
                "mode": req.mode or "std",
                "aspect_ratio": req.aspect_ratio or "16:9",
                "duration": str(req.duration or 5),
                "image_url": req.image_url,
            }

### 5.4 新增厂商步骤

1. 在 media_adapters/ 下新建文件
2. 实现 MediaAdapter 基类的 submit / poll / normalize_params
3. 加上 @register_adapter("provider", "media_type") 装饰器
4. 在数据库中添加供应商和模型配置记录
5. 上层代码、API 路由、任务管理无需改动

---

## 6. 错误处理

### 6.1 统一错误码

    class MediaErrorCode(StrEnum):
        CONTENT_POLICY = "content_policy"
        RATE_LIMITED = "rate_limited"
        QUOTA_EXCEEDED = "quota_exceeded"
        INVALID_PARAMS = "invalid_params"
        MODEL_UNAVAILABLE = "model_unavailable"
        VENDOR_ERROR = "vendor_error"
        TIMEOUT = "timeout"
        OSS_UPLOAD_FAILED = "oss_upload_failed"

### 6.2 厂商错误映射

每个适配器维护自己的映射表：

    DASHSCOPE_ERROR_MAP = {
        "DataInspectionFailed": MediaErrorCode.CONTENT_POLICY,
        "Throttling": MediaErrorCode.RATE_LIMITED,
        "InvalidParameter": MediaErrorCode.INVALID_PARAMS,
    }

---

## 7. 限流与计费

### 7.1 限流

    请求进入
      ├─ 1. 检查用户积分余额（预估消耗）
      ├─ 2. 检查供应商级 RPM 限制（Redis 滑动窗口，复用现有逻辑）
      ├─ 3. 检查模型级并发限制（视频生成有并发上限）
      └─ 4. 通过 → 提交；不通过 → 429 + Retry-After

### 7.2 计费

model_config 表扩展：

    cost_per_generation: Mapped[Decimal|None]   # 每次生成费用（图像用）
    cost_per_second: Mapped[Decimal|None]        # 每秒费用（视频按时长用）

- 图像：cost_per_generation × n
- 视频：cost_per_second × duration
- 预扣制：提交时按最大值预扣，完成后按实际结算

### 7.3 故障转移

复用现有熔断器：连续 3 次失败 → 熔断 60s → 有同类型备选模型则自动降级 → 无备选则返回错误

---

## 8. 实施阶段

### Phase 1：基础框架
- 适配器基类 + 注册机制
- media_task 数据库表 + 迁移
- MediaGateway 核心逻辑
- Celery 轮询任务

### Phase 2：图像生成
- DALL-E 3 适配器
- 通义万相适配器
- POST /v1/images/generations API
- OSS 转存

### Phase 3：视频生成
- 可灵视频适配器
- POST /v1/videos/generations API
- GET /v1/videos/generations/{task_id} 轮询
- Webhook 回调

### Phase 4：运营能力
- 管理后台模型配置界面
- 用量统计与报表
- 告警（厂商异常、积分不足）
