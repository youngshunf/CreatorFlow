# Requirements Document

## Introduction

本文档定义了将 Remotion 视频创作引擎集成到 CreatorFlow 桌面应用的功能需求。该集成将实现 AI 驱动的程序化视频创作能力，允许用户通过自然语言描述生成专业视频内容。

核心价值：
- **React 生态兼容** - CreatorFlow 已使用 React，可无缝集成 Remotion
- **AI 友好** - Remotion 官方提供 Claude Code Agent Skills，专为 AI 生成视频优化
- **程序化创作** - 代码即视频，完美契合 AI 生成场景
- **灵活渲染** - 支持本地渲染、预览服务器、多种输出格式

## Glossary

- **Video_Package**: `@creator-flow/video` 包，负责视频处理核心逻辑，包含 Remotion 配置、组件、模板和渲染器
- **Composition**: Remotion 中的视频组合单元，定义视频的尺寸、帧率、时长和内容
- **Video_Project**: 视频项目数据结构，包含配置、组件代码、素材资源和渲染历史
- **Video_Service**: Electron 主进程中的视频服务模块，负责 IPC 处理、项目管理和渲染调度
- **Video_Editor**: 渲染进程中的视频编辑器 UI 组件集合
- **Render_Worker**: 后台渲染工作进程，执行 Remotion 视频渲染任务
- **Preview_Server**: Remotion 预览服务器，提供实时视频预览能力
- **Agent_Skills**: AI Agent 视频创作工具集，供 Claude SDK 调用
- **Template**: 预设视频模板，包含社交媒体、营销视频、教程视频等类型
- **Asset**: 视频素材资源，包括图片、视频、音频、字体等

## Requirements

### Requirement 1: Video Package 核心架构

**User Story:** 作为开发者，我希望有一个独立的视频处理包，以便集中管理 Remotion 相关的配置、组件和渲染逻辑。

#### Acceptance Criteria

1. THE Video_Package SHALL 在 `CreatorFlow/packages/video/` 目录下创建，遵循 monorepo 包结构规范
2. THE Video_Package SHALL 包含 Remotion 核心依赖：`remotion ^4.0.0`、`@remotion/player`、`@remotion/cli`、`@remotion/renderer`、`@remotion/bundler`
3. THE Video_Package SHALL 导出 `RemotionRoot` 根组件，注册所有可用的 Composition
4. THE Video_Package SHALL 提供 `remotion.config.ts` 配置文件，设置视频图像格式和输出覆盖选项
5. THE Video_Package SHALL 使用 TypeScript 编写，并提供完整的类型定义导出

### Requirement 2: 视频项目数据模型

**User Story:** 作为用户，我希望视频项目有清晰的数据结构，以便管理视频配置、组件代码和素材资源。

#### Acceptance Criteria

1. THE Video_Package SHALL 定义 `VideoProject` 数据模型，使用 Zod schema 进行验证
2. WHEN 创建 Video_Project 时，THE Video_Package SHALL 包含以下必需字段：id、name、createdAt、updatedAt、config（width、height、fps、durationInFrames）
3. THE Video_Project SHALL 支持存储多个 Composition，每个包含 id、name、code（React 组件代码）、props
4. THE Video_Project SHALL 支持管理 Asset 列表，每个包含 id、type（image/video/audio/font）、name、path
5. THE Video_Project SHALL 记录渲染历史，包含 id、compositionId、outputPath、status、progress、createdAt

### Requirement 3: 基础 Composition 组件

**User Story:** 作为用户，我希望有预置的视频组件，以便快速创建常见类型的视频内容。

#### Acceptance Criteria

1. THE Video_Package SHALL 提供 `TitleAnimation` 组件，支持标题和副标题的淡入动画效果
2. THE Video_Package SHALL 提供 `Slideshow` 组件，支持图片轮播展示
3. THE Video_Package SHALL 提供 `DataVisualization` 组件，支持数据图表动画展示
4. THE Video_Package SHALL 提供 `ProductShowcase` 组件，支持产品特性展示动画
5. WHEN 创建 Composition 组件时，THE Video_Package SHALL 使用 Remotion 的 `useCurrentFrame`、`interpolate`、`spring` 等 API 实现动画
6. THE Composition 组件 SHALL 接受标准化的 props 接口，包含 title、subtitle、items、colors、logo、cta 等可选字段

### Requirement 4: 可复用 UI 组件

**User Story:** 作为开发者，我希望有可复用的视频 UI 组件，以便在不同 Composition 中共享动画效果。

#### Acceptance Criteria

1. THE Video_Package SHALL 提供 `AnimatedText` 组件，支持文字淡入、滑入、缩放等动画效果
2. THE Video_Package SHALL 提供 `Transition` 组件，支持场景切换过渡效果
3. THE Video_Package SHALL 提供 `Background` 组件，支持纯色、渐变、图片背景
4. THE Video_Package SHALL 提供 `Logo` 组件，支持 Logo 图片的动画展示
5. WHEN 使用可复用组件时，THE 组件 SHALL 通过 props 接受动画配置参数

### Requirement 5: 预设模板系统

**User Story:** 作为用户，我希望有预设的视频模板，以便快速开始创建特定类型的视频。

#### Acceptance Criteria

1. THE Video_Package SHALL 提供 `social-media` 模板，适用于短视频平台（9:16 竖屏、1:1 方形）
2. THE Video_Package SHALL 提供 `marketing` 模板，适用于产品营销视频（16:9 横屏）
3. THE Video_Package SHALL 提供 `tutorial` 模板，适用于教程讲解视频
4. WHEN 用户选择模板时，THE Video_Package SHALL 返回预配置的 Composition 代码和默认 props
5. THE 模板 SHALL 支持自定义颜色方案、字体、Logo 等品牌元素

### Requirement 6: 本地渲染引擎

**User Story:** 作为用户，我希望能够将视频渲染为 MP4 文件，以便导出和分享。

#### Acceptance Criteria

1. THE Render_Worker SHALL 使用 `@remotion/bundler` 打包视频项目代码
2. THE Render_Worker SHALL 使用 `@remotion/renderer` 的 `renderMedia` API 执行视频渲染
3. WHEN 渲染视频时，THE Render_Worker SHALL 支持三种质量预设：draft（快速预览）、standard（标准质量）、high（高质量）
4. THE Render_Worker SHALL 支持输出格式：mp4（H.264）、webm（VP9）、gif
5. WHILE 渲染进行中，THE Render_Worker SHALL 通过 IPC 向渲染进程报告进度（0-100%）
6. IF 渲染失败，THEN THE Render_Worker SHALL 返回详细的错误信息

### Requirement 7: 预览服务器

**User Story:** 作为用户，我希望能够实时预览视频效果，以便在渲染前调整内容。

#### Acceptance Criteria

1. THE Preview_Server SHALL 启动 Remotion Studio 开发服务器，提供实时预览
2. WHEN 启动预览时，THE Preview_Server SHALL 返回预览 URL 供渲染进程加载
3. THE Preview_Server SHALL 支持热重载，当组件代码变更时自动刷新预览
4. WHEN 用户关闭预览时，THE Preview_Server SHALL 正确清理资源并停止服务

### Requirement 8: Electron 主进程视频服务

**User Story:** 作为开发者，我希望主进程提供完整的视频服务 API，以便渲染进程调用。

#### Acceptance Criteria

1. THE Video_Service SHALL 在 `apps/electron/src/main/video/` 目录下创建
2. THE Video_Service SHALL 注册以下 IPC 处理器：video:create-project、video:list-projects、video:get-project、video:update-project、video:delete-project
3. THE Video_Service SHALL 注册素材管理 IPC：video:add-asset、video:remove-asset
4. THE Video_Service SHALL 注册渲染相关 IPC：video:render、video:start-preview、video:stop-preview
5. WHEN 渲染视频时，THE Video_Service SHALL 显示系统保存对话框让用户选择输出路径
6. THE Video_Service SHALL 通过 IPC 事件 `video:render-progress` 向渲染进程推送渲染进度

### Requirement 9: 视频编辑器 UI 组件

**User Story:** 作为用户，我希望有直观的视频编辑界面，以便管理和编辑视频项目。

#### Acceptance Criteria

1. THE Video_Editor SHALL 提供 `VideoEditor` 主组件，采用三栏布局：左侧项目列表、中间预览区、右侧属性面板
2. THE Video_Editor SHALL 提供 `VideoPreview` 组件，使用 `@remotion/player` 嵌入视频预览播放器
3. THE Video_Editor SHALL 提供 `VideoTimeline` 组件，显示视频时间轴和当前帧位置
4. THE Video_Editor SHALL 提供 `VideoProperties` 组件，允许编辑视频配置和组件 props
5. THE Video_Editor SHALL 提供 `VideoProjectList` 组件，显示工作区内的视频项目列表
6. THE Video_Editor SHALL 提供 `VideoTemplates` 组件，展示可用模板供用户选择
7. THE Video_Editor SHALL 提供 `VideoExport` 组件，配置导出选项并触发渲染

### Requirement 10: AI Agent 视频创作工具

**User Story:** 作为用户，我希望通过自然语言与 AI 对话来创建视频，以便降低视频创作门槛。

#### Acceptance Criteria

1. THE Agent_Skills SHALL 提供 `video_create_project` 工具，支持创建新视频项目并选择模板
2. THE Agent_Skills SHALL 提供 `video_generate_composition` 工具，根据用户描述生成 Remotion 组件代码
3. THE Agent_Skills SHALL 提供 `video_update_composition` 工具，根据用户反馈修改现有组件
4. THE Agent_Skills SHALL 提供 `video_preview` 工具，启动视频预览供用户查看效果
5. THE Agent_Skills SHALL 提供 `video_render` 工具，将视频渲染为指定格式的文件
6. THE Agent_Skills SHALL 提供 `video_add_asset` 工具，添加图片、视频、音频等素材到项目
7. WHEN AI 生成组件代码时，THE Agent_Skills SHALL 遵循 Remotion 最佳实践：使用函数组件、避免 useEffect/useState、使用 interpolate 创建动画

### Requirement 11: Agent Skills 指南文档

**User Story:** 作为 AI Agent，我需要详细的 Remotion 使用指南，以便正确生成视频组件代码。

#### Acceptance Criteria

1. THE Video_Package SHALL 在 `src/skills/` 目录下提供 `AGENTS.md` 指南文档
2. THE AGENTS.md SHALL 说明 Remotion 核心概念：Composition、useCurrentFrame、useVideoConfig、interpolate
3. THE AGENTS.md SHALL 提供代码生成规则：使用 TypeScript、函数组件、AbsoluteFill 容器
4. THE AGENTS.md SHALL 包含常用动画模式示例：淡入、滑入、弹性动画、序列动画
5. THE AGENTS.md SHALL 说明可用模板及其参数接口
6. THE AGENTS.md SHALL 说明渲染质量预设和输出格式选项

### Requirement 12: 视频项目存储

**User Story:** 作为用户，我希望视频项目能够持久化存储，以便下次打开时继续编辑。

#### Acceptance Criteria

1. THE Video_Service SHALL 将视频项目存储在工作区目录下的 `.creator-flow/video-projects/` 路径
2. WHEN 保存 Video_Project 时，THE Video_Service SHALL 将项目数据序列化为 JSON 文件
3. THE Video_Service SHALL 将项目素材文件存储在项目目录下的 `assets/` 子目录
4. WHEN 加载工作区时，THE Video_Service SHALL 自动扫描并加载所有视频项目
5. IF 项目文件损坏，THEN THE Video_Service SHALL 记录错误日志并跳过该项目

### Requirement 13: 错误处理与用户反馈

**User Story:** 作为用户，我希望在操作失败时获得清晰的错误提示，以便了解问题并采取措施。

#### Acceptance Criteria

1. IF 渲染失败，THEN THE Video_Editor SHALL 显示错误对话框，包含错误类型和建议解决方案
2. IF AI 生成的代码无法编译，THEN THE Agent_Skills SHALL 返回语法错误详情并尝试自动修复
3. IF 素材文件不存在或格式不支持，THEN THE Video_Service SHALL 返回明确的错误消息
4. WHILE 执行长时间操作时，THE Video_Editor SHALL 显示进度指示器和取消按钮
5. IF 用户取消渲染，THEN THE Render_Worker SHALL 正确清理临时文件并释放资源

### Requirement 14: 国际化支持

**User Story:** 作为用户，我希望视频编辑器支持多语言，以便使用我熟悉的语言操作。

#### Acceptance Criteria

1. THE Video_Editor SHALL 使用 CreatorFlow 现有的 i18n 系统（中文 key + 翻译文件）
2. THE Video_Editor SHALL 为所有用户可见文案提供中文和英文翻译
3. THE 模板名称和描述 SHALL 支持国际化
4. THE 错误消息 SHALL 支持国际化

### Requirement 15: Preload API 暴露

**User Story:** 作为开发者，我需要通过 preload 脚本安全地暴露视频 API 给渲染进程。

#### Acceptance Criteria

1. THE Video_Service SHALL 在 `apps/electron/src/preload/` 下创建 `video-api.ts`
2. THE video-api SHALL 通过 contextBridge 暴露类型安全的 `electronAPI.video` 接口
3. THE video-api SHALL 包含所有视频相关 IPC 调用的封装方法
4. THE video-api SHALL 提供 `onRenderProgress` 事件监听器用于接收渲染进度更新
