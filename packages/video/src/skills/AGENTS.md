# Video Creation Agent Guide

本文档为 AI Agent 提供视频创作能力的使用指南。

## 核心概念

### Remotion 基础

Remotion 是一个基于 React 的程序化视频创作框架。核心概念包括：

1. **Composition** - 视频组件，定义视频的内容和动画
2. **Frame** - 帧，视频的基本时间单位
3. **FPS** - 帧率，每秒帧数
4. **Duration** - 时长，以帧数表示

### 时间计算

```typescript
// 秒转帧
const frames = seconds * fps;

// 帧转秒
const seconds = frames / fps;

// 示例：30fps 下 10 秒 = 300 帧
const durationInFrames = 10 * 30; // 300
```

## 可用工具

### 1. video_create_project

创建新的视频项目。

**参数：**
- `name` (必需): 项目名称
- `durationInSeconds` (必需): 视频时长（秒）
- `template`: 使用的模板
- `width`: 视频宽度（默认 1920）
- `height`: 视频高度（默认 1080）
- `fps`: 帧率（默认 30）

**模板选项：**
- `blank` - 空白项目
- `title-animation` - 标题动画
- `slideshow` - 图片轮播
- `product-showcase` - 产品展示
- `social-media` - 社交媒体（横屏）
- `social-media-vertical` - 社交媒体（竖屏 9:16）
- `social-media-square` - 社交媒体（方形 1:1）
- `marketing` - 营销视频
- `tutorial` - 教程视频

**示例：**
```json
{
  "name": "产品介绍视频",
  "template": "product-showcase",
  "durationInSeconds": 30,
  "width": 1920,
  "height": 1080
}
```

### 2. video_generate_composition

根据自然语言描述生成视频组件。

**参数：**
- `projectId` (必需): 项目 ID
- `description` (必需): 视频内容描述
- `style`: 视觉风格
- `colorScheme`: 自定义配色

**风格选项：**
- `modern` - 现代风格（默认）
- `minimal` - 极简风格
- `playful` - 活泼风格
- `corporate` - 企业风格
- `cinematic` - 电影风格
- `vibrant` - 鲜艳风格
- `nature` - 自然风格

**示例：**
```json
{
  "projectId": "proj_123",
  "description": "展示三个产品特性，每个特性配有图标和简短说明，使用淡入动画",
  "style": "modern",
  "colorScheme": {
    "primary": "#6366f1",
    "secondary": "#8b5cf6",
    "background": "#1a1a2e"
  }
}
```

### 3. video_update_composition

更新现有视频组件。

**参数：**
- `projectId` (必需): 项目 ID
- `compositionId` (必需): 组件 ID
- `changes` (必需): 变更描述
- `props`: 要更新的属性

**示例：**
```json
{
  "projectId": "proj_123",
  "compositionId": "comp_001",
  "changes": "将标题改为'智能创作助手'，背景色改为深蓝色",
  "props": {
    "title": "智能创作助手",
    "colors": {
      "background": "#0f172a"
    }
  }
}
```

### 4. video_preview

启动视频预览服务器。

**参数：**
- `projectId` (必需): 项目 ID
- `compositionId`: 指定组件 ID
- `startFrame`: 起始帧

**示例：**
```json
{
  "projectId": "proj_123",
  "startFrame": 0
}
```

### 5. video_render

渲染视频到文件。

**参数：**
- `projectId` (必需): 项目 ID
- `compositionId` (必需): 组件 ID
- `outputFormat`: 输出格式（mp4/webm/gif）
- `quality`: 质量预设（draft/standard/high）

**质量预设说明：**
- `draft` - 快速预览，低质量，文件小
- `standard` - 平衡质量和大小（推荐）
- `high` - 最高质量，文件大

**示例：**
```json
{
  "projectId": "proj_123",
  "compositionId": "comp_001",
  "outputFormat": "mp4",
  "quality": "standard"
}
```

### 6. video_add_asset

添加素材到项目。

**参数：**
- `projectId` (必需): 项目 ID
- `assetPath` (必需): 素材文件路径
- `assetType` (必需): 素材类型

**素材类型：**
- `image` - 图片（png, jpg, gif, webp, svg）
- `video` - 视频（mp4, webm, mov）
- `audio` - 音频（mp3, wav, ogg, m4a）
- `font` - 字体（ttf, otf, woff, woff2）

**示例：**
```json
{
  "projectId": "proj_123",
  "assetPath": "/path/to/logo.png",
  "assetType": "image"
}
```

## 动画模式

### 基础动画

使用 Remotion 的 `interpolate` 函数创建动画：

```typescript
import { interpolate, useCurrentFrame } from 'remotion';

const frame = useCurrentFrame();

// 淡入动画（0-30帧）
const opacity = interpolate(frame, [0, 30], [0, 1], {
  extrapolateRight: 'clamp',
});

// 滑入动画
const translateY = interpolate(frame, [0, 30], [50, 0], {
  extrapolateRight: 'clamp',
});

// 缩放动画
const scale = interpolate(frame, [0, 30], [0.8, 1], {
  extrapolateRight: 'clamp',
});
```

### 弹性动画

使用 `spring` 函数创建弹性效果：

```typescript
import { spring, useCurrentFrame, useVideoConfig } from 'remotion';

const frame = useCurrentFrame();
const { fps } = useVideoConfig();

const scale = spring({
  frame,
  fps,
  config: {
    damping: 10,
    stiffness: 100,
    mass: 0.5,
  },
});
```

### 序列动画

使用 `Sequence` 组件创建时间序列：

```typescript
import { Sequence } from 'remotion';

<Sequence from={0} durationInFrames={60}>
  <TitleSlide />
</Sequence>
<Sequence from={60} durationInFrames={90}>
  <ContentSlide />
</Sequence>
<Sequence from={150} durationInFrames={60}>
  <EndSlide />
</Sequence>
```

## 模板参数

所有模板都接受以下通用参数：

```typescript
interface TemplateProps {
  title?: string;           // 主标题
  subtitle?: string;        // 副标题
  items?: Array<{           // 内容项
    title: string;
    description?: string;
    image?: string;
    icon?: string;
  }>;
  colors?: {                // 配色方案
    primary: string;
    secondary: string;
    background: string;
    text: string;
  };
  logo?: string;            // Logo 图片路径
  cta?: {                   // 行动号召
    text: string;
    url?: string;
  };
  fontFamily?: string;      // 字体
  animationStyle?: string;  // 动画风格
}
```

## 最佳实践

### 1. 选择合适的分辨率

| 平台 | 推荐分辨率 | 宽高比 |
|------|-----------|--------|
| YouTube | 1920×1080 | 16:9 |
| TikTok/Reels | 1080×1920 | 9:16 |
| Instagram Feed | 1080×1080 | 1:1 |
| Instagram Story | 1080×1920 | 9:16 |

### 2. 帧率选择

- **24fps** - 电影感
- **30fps** - 标准视频（推荐）
- **60fps** - 流畅动画

### 3. 时长建议

- **社交媒体** - 15-60秒
- **产品介绍** - 30-90秒
- **教程视频** - 2-5分钟

### 4. 动画时长

- **淡入/淡出** - 0.5-1秒（15-30帧@30fps）
- **滑动动画** - 0.3-0.8秒
- **弹性动画** - 0.5-1.5秒

## 工作流程示例

### 创建社交媒体视频

```
1. 创建项目
   video_create_project({
     name: "产品发布预告",
     template: "social-media-vertical",
     durationInSeconds: 15
   })

2. 添加素材
   video_add_asset({
     projectId: "proj_xxx",
     assetPath: "/path/to/product.png",
     assetType: "image"
   })

3. 生成内容
   video_generate_composition({
     projectId: "proj_xxx",
     description: "展示产品图片，配合标题'新品上市'和副标题'限时优惠'",
     style: "vibrant"
   })

4. 预览
   video_preview({ projectId: "proj_xxx" })

5. 渲染
   video_render({
     projectId: "proj_xxx",
     compositionId: "comp_xxx",
     outputFormat: "mp4",
     quality: "high"
   })
```

## 错误处理

常见错误及解决方案：

| 错误 | 原因 | 解决方案 |
|------|------|----------|
| BUNDLE_FAILED | 组件代码错误 | 检查 React 组件语法 |
| COMPOSITION_NOT_FOUND | 组件 ID 不存在 | 确认组件 ID 正确 |
| RENDER_FAILED | 渲染过程出错 | 检查素材路径和配置 |
| ASSET_NOT_FOUND | 素材文件不存在 | 确认文件路径正确 |

## 注意事项

1. **素材路径** - 使用绝对路径或相对于项目目录的路径
2. **内存管理** - 大型视频项目可能需要较多内存
3. **渲染时间** - 高质量渲染可能需要较长时间
4. **格式兼容** - GIF 格式不支持音频
