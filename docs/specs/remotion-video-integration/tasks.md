# Implementation Plan: Remotion Video Integration

## Overview

本实施计划将 Remotion 视频创作引擎集成到 CreatorFlow 桌面应用。采用分阶段实施策略，从基础架构开始，逐步构建渲染引擎、UI 组件和 AI 集成能力。

## Tasks

- [x] 1. 创建 @creator-flow/video 包基础架构
  - [x] 1.1 初始化包结构和配置
    - 在 `CreatorFlow/packages/video/` 创建目录结构
    - 创建 `package.json`，添加 Remotion 依赖
    - 创建 `tsconfig.json`，配置 TypeScript
    - 创建 `remotion.config.ts`，配置 Remotion
    - _Requirements: 1.1, 1.2, 1.4_
  
  - [x] 1.2 实现核心类型定义
    - 创建 `src/types.ts`，定义 VideoProject、Composition、Asset、RenderHistory 等 Zod schema
    - 创建 `src/index.ts`，导出所有公共 API
    - _Requirements: 2.1, 2.2, 2.3, 2.4, 2.5_
  
  - [ ]* 1.3 编写类型定义属性测试
    - **Property 1: VideoProject Schema Validation**
    - **Validates: Requirements 2.1, 2.2, 2.3, 2.4, 2.5**

- [x] 2. 实现基础 Composition 组件
  - [x] 2.1 创建 RemotionRoot 根组件
    - 创建 `src/Root.tsx`，注册所有 Composition
    - _Requirements: 1.3_
  
  - [x] 2.2 实现 TitleAnimation 组件
    - 创建 `src/compositions/TitleAnimation.tsx`
    - 实现标题和副标题的淡入动画效果
    - _Requirements: 3.1, 3.5, 3.6_
  
  - [x] 2.3 实现 Slideshow 组件
    - 创建 `src/compositions/Slideshow.tsx`
    - 实现图片轮播展示功能
    - _Requirements: 3.2, 3.5, 3.6_
  
  - [x] 2.4 实现 DataVisualization 组件
    - 创建 `src/compositions/DataVisualization.tsx`
    - 实现数据图表动画展示
    - _Requirements: 3.3, 3.5, 3.6_
  
  - [x] 2.5 实现 ProductShowcase 组件
    - 创建 `src/compositions/ProductShowcase.tsx`
    - 实现产品特性展示动画
    - _Requirements: 3.4, 3.5, 3.6_

- [x] 3. 实现可复用 UI 组件
  - [x] 3.1 实现 AnimatedText 组件
    - 创建 `src/components/AnimatedText.tsx`
    - 支持淡入、滑入、缩放等动画效果
    - _Requirements: 4.1, 4.5_
  
  - [x] 3.2 实现 Transition 组件
    - 创建 `src/components/Transition.tsx`
    - 支持场景切换过渡效果
    - _Requirements: 4.2, 4.5_
  
  - [x] 3.3 实现 Background 组件
    - 创建 `src/components/Background.tsx`
    - 支持纯色、渐变、图片背景
    - _Requirements: 4.3, 4.5_
  
  - [x] 3.4 实现 Logo 组件
    - 创建 `src/components/Logo.tsx`
    - 支持 Logo 图片的动画展示
    - _Requirements: 4.4, 4.5_

- [x] 4. Checkpoint - 验证基础组件
  - 运行 `bun run dev` 启动 Remotion Studio
  - 验证所有 Composition 组件可正常预览
  - 确保所有测试通过，如有问题请询问用户

- [x] 5. 实现预设模板系统
  - [x] 5.1 定义模板接口和类型
    - 创建 `src/templates/types.ts`
    - 定义 VideoTemplate、TemplateProps 接口
    - _Requirements: 5.4_
  
  - [x] 5.2 实现 social-media 模板
    - 创建 `src/templates/social-media.ts`
    - 配置 9:16 竖屏和 1:1 方形预设
    - _Requirements: 5.1, 5.5_
  
  - [x] 5.3 实现 marketing 模板
    - 创建 `src/templates/marketing.ts`
    - 配置 16:9 横屏预设
    - _Requirements: 5.2, 5.5_
  
  - [x] 5.4 实现 tutorial 模板
    - 创建 `src/templates/tutorial.ts`
    - 配置教程视频预设
    - _Requirements: 5.3, 5.5_
  
  - [ ]* 5.5 编写模板属性测试
    - **Property 7: Template Return Value Validation**
    - **Validates: Requirements 5.4, 5.5**

- [x] 6. 实现 Electron 主进程视频服务
  - [x] 6.1 创建视频服务目录结构
    - 在 `apps/electron/src/main/video/` 创建目录
    - 创建 `index.ts` 入口文件
    - _Requirements: 8.1_
  
  - [x] 6.2 实现 ProjectManager
    - 创建 `project-manager.ts`
    - 实现项目的 CRUD 操作
    - 实现素材管理功能
    - 实现项目存储和加载逻辑
    - _Requirements: 12.1, 12.2, 12.3, 12.4_
  
  - [ ]* 6.3 编写项目持久化属性测试
    - **Property 2: Project Persistence Round-Trip**
    - **Validates: Requirements 12.1, 12.2, 12.3, 12.4**
  
  - [ ]* 6.4 编写损坏项目处理属性测试
    - **Property 8: Corrupted Project Handling**
    - **Validates: Requirements 12.5**

- [x] 7. 实现渲染引擎
  - [x] 7.1 实现 RenderWorker
    - 创建 `render-worker.ts`
    - 实现 bundle 和 renderMedia 调用
    - 实现质量预设配置
    - 实现进度报告机制
    - _Requirements: 6.1, 6.2, 6.3, 6.4_
  
  - [ ]* 7.2 编写渲染进度属性测试
    - **Property 5: Render Progress Reporting**
    - **Validates: Requirements 6.5**
  
  - [ ]* 7.3 编写渲染错误处理属性测试
    - **Property 6: Render Error Handling**
    - **Validates: Requirements 6.6**
  
  - [x] 7.4 实现 PreviewServer
    - 创建 `preview-server.ts`
    - 实现 Remotion Studio 启动和停止
    - 实现预览 URL 返回
    - _Requirements: 7.1, 7.2, 7.4_

- [x] 8. 实现 IPC 处理器
  - [x] 8.1 注册项目管理 IPC
    - 创建 `ipc-handlers.ts`
    - 注册 video:create-project、video:list-projects、video:get-project、video:update-project、video:delete-project
    - _Requirements: 8.2_
  
  - [x] 8.2 注册素材管理 IPC
    - 注册 video:add-asset、video:remove-asset
    - 实现素材验证逻辑
    - _Requirements: 8.3, 13.3_
  
  - [ ]* 8.3 编写无效素材错误处理属性测试
    - **Property 9: Invalid Asset Error Handling**
    - **Validates: Requirements 13.3**
  
  - [x] 8.4 注册渲染和预览 IPC
    - 注册 video:render、video:cancel-render、video:start-preview、video:stop-preview
    - 实现渲染进度事件推送
    - _Requirements: 8.4, 8.5, 8.6_
  
  - [ ]* 8.5 编写 IPC 通道完整性属性测试
    - **Property 3: IPC Channel Completeness**
    - **Validates: Requirements 8.2, 8.3, 8.4, 8.6**

- [x] 9. Checkpoint - 验证主进程服务
  - 确保所有 IPC 处理器正确注册
  - 测试项目创建、保存、加载流程
  - 确保所有测试通过，如有问题请询问用户

- [x] 10. 实现 Preload API
  - [x] 10.1 创建 video-api.ts
    - 在 `apps/electron/src/preload/` 创建 `video-api.ts`
    - 通过 contextBridge 暴露 VideoAPI 接口
    - 封装所有视频相关 IPC 调用
    - _Requirements: 15.1, 15.2, 15.3_
  
  - [x] 10.2 实现事件监听器
    - 实现 onRenderProgress 事件监听
    - _Requirements: 15.4_
  
  - [ ]* 10.3 编写 Preload API 完整性属性测试
    - **Property 10: Preload API Completeness**
    - **Validates: Requirements 15.3**

- [x] 11. 实现视频编辑器 UI 组件
  - [x] 11.1 创建 VideoEditor 主组件
    - 在 `apps/electron/src/renderer/components/video/` 创建目录
    - 创建 `VideoEditor.tsx`，实现三栏布局
    - _Requirements: 9.1_
  
  - [x] 11.2 实现 VideoPreview 组件
    - 创建 `VideoPreview.tsx`
    - 集成 @remotion/player
    - _Requirements: 9.2_
  
  - [x] 11.3 实现 VideoTimeline 组件
    - 创建 `VideoTimeline.tsx`
    - 实现时间轴和帧位置显示
    - _Requirements: 9.3_
  
  - [x] 11.4 实现 VideoProperties 组件
    - 创建 `VideoProperties.tsx`
    - 实现配置和 props 编辑
    - _Requirements: 9.4_
  
  - [x] 11.5 实现 VideoProjectList 组件
    - 创建 `VideoProjectList.tsx`
    - 实现项目列表展示和选择
    - _Requirements: 9.5_
  
  - [x] 11.6 实现 VideoTemplates 组件
    - 创建 `VideoTemplates.tsx`
    - 实现模板选择器
    - _Requirements: 9.6_
  
  - [x] 11.7 实现 VideoExport 组件
    - 创建 `VideoExport.tsx`
    - 实现导出选项配置
    - _Requirements: 9.7_

- [x] 12. Checkpoint - 验证 UI 组件
  - 在 Electron 应用中测试视频编辑器
  - 验证预览、时间轴、属性面板功能
  - 确保所有测试通过，如有问题请询问用户

- [x] 13. 实现 AI Agent Skills
  - [x] 13.1 定义视频创作工具
    - 创建 `packages/video/src/skills/video-tools.ts`
    - 定义 video_create_project、video_generate_composition、video_update_composition、video_preview、video_render、video_add_asset 工具
    - _Requirements: 10.1, 10.2, 10.3, 10.4, 10.5, 10.6_
  
  - [ ]* 13.2 编写 Agent Skills 完整性属性测试
    - **Property 4: Agent Skills Tool Completeness**
    - **Validates: Requirements 10.1, 10.2, 10.3, 10.4, 10.5, 10.6**
  
  - [x] 13.3 创建 AGENTS.md 指南文档
    - 创建 `packages/video/src/skills/AGENTS.md`
    - 说明 Remotion 核心概念和代码生成规则
    - 提供动画模式示例和模板参数说明
    - _Requirements: 11.1, 11.2, 11.3, 11.4, 11.5, 11.6_

- [x] 14. 实现国际化支持
  - [x] 14.1 添加视频编辑器翻译
    - 在 `apps/electron/src/renderer/locales/` 添加视频相关翻译
    - 使用 t() 函数包装所有用户可见文案
    - _Requirements: 14.1, 14.2_
  
  - [x] 14.2 添加模板和错误消息翻译
    - 添加模板名称和描述的翻译
    - 添加错误消息的翻译
    - _Requirements: 14.3, 14.4_

- [x] 15. 集成和连接
  - [x] 15.1 在主进程入口注册视频服务
    - 修改 `apps/electron/src/main/index.ts`
    - 导入并初始化视频服务
    - 注册 IPC 处理器
  
  - [x] 15.2 在 preload 入口暴露视频 API
    - 修改 `apps/electron/src/preload/index.ts`
    - 导入并暴露 videoAPI
  
  - [x] 15.3 添加视频页面路由
    - 在渲染进程添加视频创作页面路由
    - 集成 VideoEditor 组件

- [x] 16. Final Checkpoint - 完整功能验证
  - 运行所有单元测试和属性测试
  - 测试完整的视频创作流程：创建项目 → 选择模板 → 预览 → 渲染
  - 测试 AI Agent 视频创作工具
  - 确保所有测试通过，如有问题请询问用户

## Notes

- 标记为 `*` 的任务是可选的测试任务，可以跳过以加快 MVP 开发
- 每个任务都引用了具体的需求编号，确保可追溯性
- Checkpoint 任务用于阶段性验证，确保增量开发的正确性
- 属性测试验证通用正确性属性，单元测试验证具体示例和边界情况
