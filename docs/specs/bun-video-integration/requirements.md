# 需求文档

## 简介

本功能旨在改进 CreatorFlow 桌面应用的视频功能架构，使用内置的 Bun 运行时来运行视频相关功能（MCP Video Server、预览服务器、渲染工作器），而不是依赖外部的 `npx` 命令或动态 import。这将减少对外部 Node.js/npm 的依赖，提升用户体验和应用的独立性。

## 术语表

- **Bun_Runtime**: 内置在桌面应用 `vendor/bun/` 目录下的 Bun 可执行文件，用于运行 JavaScript/TypeScript 代码
- **MCP_Video_Server**: 基于 FastMCP 的视频服务器，提供项目管理、素材管理、渲染和预览等 16 个工具
- **Preview_Server**: 视频预览服务器，用于启动 Remotion Studio 进行实时预览
- **Render_Worker**: 视频渲染工作器，使用 Remotion 的 bundler 和 renderer API 进行视频渲染
- **Video_Service_Manager**: 视频服务管理器，负责启动、停止和管理所有视频相关子进程
- **Electron_Main_Process**: Electron 主进程，负责管理窗口和子进程

## 需求

### 需求 1：Bun 路径解析服务

**用户故事：** 作为桌面应用开发者，我希望有一个统一的 Bun 路径解析服务，以便在不同模块中复用 Bun 可执行文件路径。

#### 验收标准

1. THE Bun_Path_Resolver SHALL 提供获取内置 Bun 可执行文件路径的方法
2. WHEN 应用在打包模式运行时，THE Bun_Path_Resolver SHALL 返回 `vendor/bun/` 目录下的 Bun 路径
3. WHEN 应用在开发模式运行时，THE Bun_Path_Resolver SHALL 返回系统 PATH 中的 `bun` 命令
4. WHEN 在 Windows 平台时，THE Bun_Path_Resolver SHALL 从 `process.resourcesPath` 解析 `bun.exe`
5. WHEN 在 macOS 或 Linux 平台时，THE Bun_Path_Resolver SHALL 从应用基础路径解析 `bun`
6. IF Bun 可执行文件不存在，THEN THE Bun_Path_Resolver SHALL 抛出描述性错误

### 需求 2：视频服务管理器

**用户故事：** 作为桌面应用用户，我希望视频相关服务能够自动管理，以便我可以无缝使用视频创作功能。

#### 验收标准

1. THE Video_Service_Manager SHALL 使用内置 Bun 运行时启动 MCP Video Server
2. THE Video_Service_Manager SHALL 管理 MCP Video Server 的生命周期（启动、停止、重启）
3. WHEN 桌面应用启动时，THE Video_Service_Manager SHALL 按需启动视频服务
4. WHEN 桌面应用关闭时，THE Video_Service_Manager SHALL 优雅地停止所有视频服务子进程
5. IF MCP Video Server 进程意外退出，THEN THE Video_Service_Manager SHALL 记录错误并支持重启
6. THE Video_Service_Manager SHALL 提供获取 MCP Video Server 连接信息的方法

### 需求 3：预览服务器 Bun 集成

**用户故事：** 作为视频创作者，我希望预览服务器使用内置运行时启动，以便我不需要安装额外的 Node.js 环境。

#### 验收标准

1. WHEN 启动预览服务器时，THE Preview_Server SHALL 使用内置 Bun 运行时而非 `npx` 命令
2. THE Preview_Server SHALL 通过 Bun 子进程启动 Remotion Studio
3. WHEN 预览服务器启动成功时，THE Preview_Server SHALL 返回预览 URL 和端口
4. WHEN 预览服务器启动失败时，THE Preview_Server SHALL 返回详细的错误信息
5. THE Preview_Server SHALL 支持同时运行多个项目的预览服务器
6. WHEN 停止预览服务器时，THE Preview_Server SHALL 正确清理子进程资源

### 需求 4：渲染工作器 Bun 集成

**用户故事：** 作为视频创作者，我希望视频渲染使用内置运行时执行，以便渲染过程更加稳定和可控。

#### 验收标准

1. THE Render_Worker SHALL 使用内置 Bun 运行时执行渲染任务
2. WHEN 渲染视频时，THE Render_Worker SHALL 通过 Bun 子进程调用 Remotion 渲染 API
3. THE Render_Worker SHALL 支持进度回调，报告渲染状态（bundling、preparing、rendering、completed、failed）
4. THE Render_Worker SHALL 支持取消正在进行的渲染任务
5. WHEN 渲染失败时，THE Render_Worker SHALL 返回分类的错误信息（BUNDLE_FAILED、COMPOSITION_NOT_FOUND、RENDER_FAILED、OUTPUT_WRITE_FAILED、CANCELLED）
6. THE Render_Worker SHALL 支持多种输出格式（mp4、webm、gif）和质量预设

### 需求 5：MCP Video Server 集成

**用户故事：** 作为 AI Agent 用户，我希望通过 MCP 协议访问视频创作功能，以便 AI 可以帮助我创建视频。

#### 验收标准

1. THE Electron_Main_Process SHALL 使用内置 Bun 启动 MCP Video Server 作为子进程
2. THE MCP_Video_Server SHALL 通过 stdio 传输模式与 Electron 主进程通信
3. WHEN MCP Video Server 启动时，THE Video_Service_Manager SHALL 记录可用工具列表
4. THE MCP_Video_Server SHALL 提供项目管理工具（create_project、list_projects、get_project、update_project、delete_project）
5. THE MCP_Video_Server SHALL 提供素材管理工具（add_asset、remove_asset、list_assets）
6. THE MCP_Video_Server SHALL 提供渲染管理工具（render_video、get_render_status、cancel_render）
7. THE MCP_Video_Server SHALL 提供预览管理工具（start_preview、stop_preview）

### 需求 6：进程间通信

**用户故事：** 作为桌面应用开发者，我希望主进程与视频服务子进程之间有可靠的通信机制。

#### 验收标准

1. THE Video_Service_Manager SHALL 通过 stdio 管道与 MCP Video Server 通信
2. WHEN 子进程输出日志时，THE Video_Service_Manager SHALL 将日志转发到主进程日志系统
3. WHEN 子进程发生错误时，THE Video_Service_Manager SHALL 捕获 stderr 输出并记录
4. THE Video_Service_Manager SHALL 支持向子进程发送信号（SIGTERM、SIGKILL）
5. IF 子进程在指定时间内未响应 SIGTERM，THEN THE Video_Service_Manager SHALL 发送 SIGKILL 强制终止

### 需求 7：错误处理与恢复

**用户故事：** 作为桌面应用用户，我希望视频服务在出错时能够自动恢复，以便我的工作不会中断。

#### 验收标准

1. IF MCP Video Server 启动失败，THEN THE Video_Service_Manager SHALL 记录详细错误并通知用户
2. IF 预览服务器启动超时，THEN THE Preview_Server SHALL 返回超时错误
3. IF 渲染过程中发生错误，THEN THE Render_Worker SHALL 返回分类的错误信息
4. THE Video_Service_Manager SHALL 支持配置最大重试次数
5. WHEN 服务重启时，THE Video_Service_Manager SHALL 保持之前的配置状态
