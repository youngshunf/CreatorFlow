# 实现计划：Bun 视频集成

## 概述

本实现计划将 CreatorFlow 桌面应用的视频功能从依赖外部 `npx` 命令迁移到使用内置 Bun 运行时。实现分为 6 个主要阶段：Bun 路径解析、视频服务管理器、预览服务器重构、渲染工作器重构、IPC 集成和测试。

## 任务

- [x] 1. 创建 Bun 路径解析模块
  - [x] 1.1 创建 `apps/electron/src/main/bun-path.ts` 文件
    - 从 `sessions.ts` 提取 Bun 路径解析逻辑
    - 实现 `BunPathResolver` 接口
    - 支持打包模式和开发模式
    - 支持 Windows、macOS、Linux 平台
    - _Requirements: 1.1, 1.2, 1.3, 1.4, 1.5, 1.6_
  
  - [ ]* 1.2 编写 Bun 路径解析的属性测试
    - **Property 1: Bun 路径解析的完整性**
    - **Validates: Requirements 1.1, 1.6**
  
  - [x] 1.3 重构 `sessions.ts` 使用新的 BunPathResolver
    - 移除重复的 Bun 路径解析代码
    - 导入并使用 `createBunPathResolver()`
    - _Requirements: 1.1_

- [x] 2. 创建视频服务管理器
  - [x] 2.1 创建 `apps/electron/src/main/video/service-manager.ts` 文件
    - 实现 `VideoServiceManager` 接口
    - 实现 MCP Video Server 子进程管理
    - 实现服务状态跟踪
    - 实现错误处理和重试逻辑
    - _Requirements: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6_
  
  - [ ]* 2.2 编写服务生命周期的属性测试
    - **Property 2: 服务生命周期状态一致性**
    - **Validates: Requirements 2.2, 2.5**
  
  - [x] 2.3 创建 `apps/electron/src/main/video/index.ts` 导出文件
    - 导出 VideoServiceManager
    - 导出相关类型
    - _Requirements: 2.1_

- [x] 3. 检查点 - 确保基础模块测试通过
  - 确保所有测试通过，如有问题请询问用户。

- [x] 4. 重构预览服务器
  - [x] 4.1 修改 `apps/electron/src/main/video/preview-server.ts`
    - 注入 BunPathResolver 依赖
    - 将 `npx remotion studio` 替换为 Bun 子进程
    - 更新 `startRemotionStudio` 方法使用 Bun
    - 保持现有的端口管理和进程清理逻辑
    - _Requirements: 3.1, 3.2, 3.3, 3.4, 3.5, 3.6_
  
  - [ ]* 4.2 编写预览服务器的属性测试
    - **Property 3: 预览服务器的多实例管理**
    - **Validates: Requirements 3.3, 3.4, 3.5, 3.6**

- [x] 5. 重构渲染工作器
  - [x] 5.1 创建独立的渲染脚本 `apps/mcp-video/src/render-script.ts`
    - 接收命令行参数（项目路径、组合 ID、输出路径、质量、格式）
    - 调用 Remotion bundler 和 renderer API
    - 通过 stdout 输出 JSON 格式的进度信息
    - _Requirements: 4.1, 4.2, 4.3, 4.6_
  
  - [x] 5.2 修改 `apps/electron/src/main/video/render-worker.ts`
    - 注入 BunPathResolver 依赖
    - 将动态 import 替换为 Bun 子进程调用渲染脚本
    - 解析子进程 stdout 获取进度信息
    - 实现取消渲染（发送 SIGTERM）
    - _Requirements: 4.1, 4.2, 4.3, 4.4, 4.5_
  
  - [ ]* 5.3 编写渲染进度的属性测试
    - **Property 4: 渲染进度报告的顺序性**
    - **Validates: Requirements 4.3, 4.4, 4.5**
  
  - [ ]* 5.4 编写渲染输出格式的属性测试
    - **Property 5: 渲染输出格式的正确性**
    - **Validates: Requirements 4.6**

- [x] 6. 检查点 - 确保视频功能测试通过
  - 确保所有测试通过，如有问题请询问用户。

- [x] 7. 集成 MCP Video Server
  - [x] 7.1 在 VideoServiceManager 中实现 MCP Server 启动逻辑
    - 使用 BunPathResolver 获取 Bun 路径
    - 启动 `apps/mcp-video/src/index.ts` 作为子进程
    - 配置 stdio 传输模式
    - 记录可用工具列表
    - _Requirements: 5.1, 5.2, 5.3_
  
  - [x] 7.2 实现进程间通信
    - 转发子进程 stdout 到日志系统
    - 捕获 stderr 并记录错误
    - 实现信号发送（SIGTERM、SIGKILL）
    - 实现超时强制终止
    - _Requirements: 6.1, 6.2, 6.3, 6.4, 6.5_
  
  - [ ]* 7.3 编写进程通信的属性测试
    - **Property 6: 进程通信的可靠性**
    - **Validates: Requirements 6.2, 6.3, 6.4, 6.5**

- [x] 8. 实现错误处理和恢复
  - [x] 8.1 实现错误分类和处理
    - 定义 VideoServiceError 枚举
    - 实现错误分类逻辑
    - 实现用户通知机制
    - _Requirements: 7.1, 7.2, 7.3_
  
  - [x] 8.2 实现服务恢复机制
    - 实现最大重试次数配置
    - 实现配置状态保持
    - 实现自动重启逻辑
    - _Requirements: 7.4, 7.5_
  
  - [ ]* 8.3 编写错误恢复的属性测试
    - **Property 7: 错误恢复的幂等性**
    - **Validates: Requirements 7.1, 7.2, 7.4, 7.5**

- [x] 9. 添加 IPC 通道
  - [x] 9.1 在 `apps/electron/src/main/ipc.ts` 中添加视频服务 IPC 处理
    - 添加 VIDEO_IPC_CHANNELS 定义
    - 实现 start/stop/status 处理器
    - 实现预览和渲染相关处理器
    - _Requirements: 5.1, 5.2_
  
  - [x] 9.2 在 `apps/electron/src/preload/index.ts` 中暴露视频 API
    - 添加 videoService 命名空间
    - 暴露 startService、stopService、getStatus 方法
    - 暴露 startPreview、stopPreview 方法
    - 暴露 startRender、cancelRender 方法
    - _Requirements: 5.1_

- [x] 10. 检查点 - 确保集成测试通过
  - 确保所有测试通过，如有问题请询问用户。

- [x] 11. 更新应用初始化
  - [x] 11.1 在应用启动时初始化 VideoServiceManager
    - 在 `apps/electron/src/main/index.ts` 中创建 VideoServiceManager 实例
    - 注册应用退出时的清理逻辑
    - _Requirements: 2.3, 2.4_
  
  - [x] 11.2 更新 electron-builder 配置（如需要）
    - 确保 MCP Video Server 相关文件被打包
    - 确保 Remotion 依赖被正确处理
    - _Requirements: 5.1_

- [x] 12. 最终检查点 - 确保所有测试通过
  - 确保所有测试通过，如有问题请询问用户。

## 备注

- 标记 `*` 的任务为可选任务，可以跳过以加快 MVP 开发
- 每个任务都引用了具体的需求以便追溯
- 检查点确保增量验证
- 属性测试验证通用正确性属性
- 单元测试验证具体示例和边界情况
