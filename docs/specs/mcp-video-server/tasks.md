# Implementation Plan: MCP Video Server

## Overview

本实现计划将 MCP Video Server 分解为可增量执行的任务，每个任务构建在前一个任务的基础上。使用 TypeScript + Bun 运行时，基于 FastMCP 框架构建。

## Tasks

- [x] 1. 项目初始化和基础设施
  - [x] 1.1 创建 `CreatorFlow/apps/mcp-video/` 目录结构
    - 创建 package.json，配置 Bun 运行时和依赖
    - 创建 tsconfig.json，配置 TypeScript
    - 添加依赖：fastmcp, zod, @remotion/bundler, @remotion/renderer
    - _Requirements: 8.1, 8.2_
  
  - [x] 1.2 创建类型定义 `src/types/index.ts`
    - 复用 @creator-flow/video 的类型定义
    - 定义 MCP 响应类型（SuccessResponse, ErrorResponse）
    - 定义错误代码枚举
    - _Requirements: 8.4, 9.1_

  - [x] 1.3 创建错误处理模块 `src/types/errors.ts`
    - 实现 MCPError 类
    - 实现错误工厂函数（createNotFoundError, createValidationError 等）
    - _Requirements: 9.1, 9.4_

- [x] 2. 项目存储服务
  - [x] 2.1 创建路径工具 `src/utils/paths.ts`
    - 实现 getVideoProjectsDir(workspacePath) 返回 `{workspacePath}/视频创作/`
    - 实现 getProjectPath(workspacePath, projectName) 返回项目目录
    - 实现 getAssetsPath(workspacePath, projectName) 返回素材目录
    - 实现 getOutputPath(workspacePath, projectName) 返回输出目录
    - _Requirements: 1.5, 1.7_

  - [x] 2.2 创建验证器 `src/utils/validators.ts`
    - 实现项目名称验证（非空、有效字符）
    - 实现素材类型验证（扩展名检查）
    - 实现工作区路径验证
    - _Requirements: 1.4, 4.2_

  - [x] 2.3 实现 ProjectStore 服务 `src/services/project-store.ts`
    - 实现 createProject() 方法
    - 实现 getProject() 方法
    - 实现 listProjects() 方法
    - 实现 updateProject() 方法
    - 实现 deleteProject() 方法
    - 实现 JSON 序列化/反序列化
    - _Requirements: 1.1, 1.5, 2.1, 2.2, 2.4, 3.1, 3.2, 3.3_

  - [ ]* 2.4 编写 ProjectStore 属性测试
    - **Property 1: Project Creation Produces Valid Persisted Project**
    - **Property 2: Project Query Operations Return Correct Data**
    - **Property 3: Project Mutations Work Correctly**
    - **Property 6: Project Serialization Round-Trip**
    - **Validates: Requirements 1.1, 1.5, 1.7, 2.1, 2.2, 2.4, 3.1, 3.2, 3.3, 5.7**

- [x] 3. Checkpoint - 确保项目存储服务测试通过
  - 运行所有测试，确保通过
  - 如有问题，询问用户

- [x] 4. 素材管理服务
  - [x] 4.1 扩展 ProjectStore 添加素材管理方法
    - 实现 addAsset() 方法（复制文件、注册素材）
    - 实现 removeAsset() 方法（删除注册、可选删除文件）
    - 实现 listAssets() 方法（按类型分组）
    - _Requirements: 4.1, 4.4, 4.5_

  - [ ]* 4.2 编写素材管理属性测试
    - **Property 4: Asset Management Operations Work Correctly**
    - **Validates: Requirements 4.1, 4.2, 4.4, 4.5**

- [x] 5. 组合管理服务
  - [x] 5.1 扩展 ProjectStore 添加组合管理方法
    - 实现 addComposition() 方法
    - 实现 updateComposition() 方法
    - 实现 removeComposition() 方法
    - _Requirements: 7.1, 7.3, 7.4_

  - [ ]* 5.2 编写组合管理属性测试
    - **Property 5: Composition Management Operations Work Correctly**
    - **Validates: Requirements 7.1, 7.3, 7.4**

- [x] 6. 渲染引擎服务
  - [x] 6.1 实现 RenderEngine 服务 `src/services/render-engine.ts`
    - 实现 Remotion bundle 构建
    - 实现 renderMedia 调用
    - 实现质量预设配置
    - 实现进度回调
    - _Requirements: 5.1, 5.2, 5.3, 5.4, 5.5, 5.6_

- [x] 7. 预览服务器服务
  - [x] 7.1 实现 PreviewServerManager 服务 `src/services/preview-server.ts`
    - 实现 start() 方法（启动 Remotion Studio）
    - 实现 stop() 方法（停止服务器）
    - 实现活跃服务器跟踪
    - 实现端口查找
    - _Requirements: 6.1, 6.3, 6.4, 6.5_

  - [ ]* 7.2 编写预览服务器属性测试
    - **Property 9: Preview Server Duplicate Prevention**
    - **Validates: Requirements 6.5**

- [x] 8. Checkpoint - 确保所有服务测试通过
  - 运行所有测试，确保通过
  - 如有问题，询问用户

- [x] 9. MCP 工具实现
  - [x] 9.1 实现项目管理工具 `src/tools/project.ts`
    - 实现 video_create_project 工具
    - 实现 video_list_projects 工具
    - 实现 video_get_project 工具
    - 实现 video_update_project 工具
    - 实现 video_delete_project 工具
    - _Requirements: 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 2.1, 2.2, 2.3, 2.4, 3.1, 3.2, 3.3, 3.4_

  - [x] 9.2 实现素材管理工具 `src/tools/asset.ts`
    - 实现 video_add_asset 工具
    - 实现 video_remove_asset 工具
    - 实现 video_list_assets 工具
    - _Requirements: 4.1, 4.2, 4.3, 4.4, 4.5_

  - [x] 9.3 实现组合管理工具 `src/tools/composition.ts`
    - 实现 video_add_composition 工具
    - 实现 video_update_composition 工具
    - 实现 video_remove_composition 工具
    - _Requirements: 7.1, 7.2, 7.3, 7.4, 7.5_

  - [x] 9.4 实现渲染工具 `src/tools/render.ts`
    - 实现 video_render 工具
    - _Requirements: 5.1, 5.2, 5.3, 5.4, 5.5, 5.6_

  - [x] 9.5 实现预览工具 `src/tools/preview.ts`
    - 实现 video_preview_start 工具
    - 实现 video_preview_stop 工具
    - _Requirements: 6.1, 6.3, 6.4, 6.5_

  - [x] 9.6 实现模板工具 `src/tools/template.ts`
    - 实现 video_list_templates 工具
    - 集成 @creator-flow/video 模板
    - _Requirements: 10.1, 10.2, 10.3_

- [x] 10. MCP 服务器集成
  - [x] 10.1 创建 FastMCP 服务器 `src/server.ts`
    - 配置 FastMCP 实例
    - 注册所有工具
    - 配置 stdio 和 HTTP 传输
    - _Requirements: 8.1, 8.2, 8.3, 8.5_

  - [x] 10.2 创建入口文件 `src/index.ts`
    - 解析命令行参数
    - 启动服务器
    - 打印启动信息
    - _Requirements: 8.5_

  - [ ]* 10.3 编写 MCP 协议合规性测试
    - **Property 7: Error Handling Returns Structured Errors**
    - **Property 8: MCP Protocol Compliance**
    - **Validates: Requirements 8.3, 8.4, 9.1, 9.4**

- [x] 11. 文档和配置
  - [x] 11.1 创建 README.md
    - 安装说明
    - 使用说明（stdio/HTTP 模式）
    - 工具列表和参数说明
    - 配置示例（Claude Desktop、Kiro）

  - [x] 11.2 添加 package.json 脚本
    - `dev`: 开发模式启动
    - `start`: 生产模式启动
    - `test`: 运行测试
    - `build`: 构建（如需要）

- [x] 12. Final Checkpoint - 确保所有测试通过
  - 运行完整测试套件
  - 验证 stdio 和 HTTP 模式
  - 如有问题，询问用户

## Notes

- 标记 `*` 的任务为可选测试任务，可跳过以加快 MVP 开发
- 每个任务引用具体的需求条款以确保可追溯性
- Checkpoint 任务用于验证增量进度
- 属性测试使用 fast-check 库，每个测试至少运行 100 次迭代
