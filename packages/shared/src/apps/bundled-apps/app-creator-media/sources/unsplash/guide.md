# Unsplash

高质量免费图片素材库，提供数百万张精美图片供创作者使用。

## Scope

- 搜索和浏览高质量图片
- 下载图片用于创作
- 获取图片详情和作者信息
- 收藏和管理图片集合

## Guidelines

### 使用规范
1. 遵守 Unsplash 使用条款，注明图片来源
2. 商业用途需确认许可证类型
3. 批量下载请控制请求频率

### API 认证
1. 访问 [Unsplash Developers](https://unsplash.com/developers) 注册应用
2. 获取 Access Key 作为 API 密钥
3. 在 Authorization 头中使用 `Client-ID {access_key}` 格式

### 常用端点
- `GET /search/photos?query={关键词}` - 搜索图片
- `GET /photos/{id}` - 获取图片详情
- `GET /photos/random` - 获取随机图片
- `GET /users/{username}/photos` - 获取用户图片

## Context

Unsplash API 提供免费的高质量图片访问。注意：
- 免费版每小时限制 50 次请求
- 生产版每小时可达 5000 次请求
- 图片使用需遵守 Unsplash 许可协议

## API Notes

响应格式为 JSON，图片 URL 在 `urls` 字段中：
- `urls.raw` - 原始尺寸
- `urls.full` - 全尺寸
- `urls.regular` - 1080px 宽
- `urls.small` - 400px 宽
- `urls.thumb` - 200px 宽
