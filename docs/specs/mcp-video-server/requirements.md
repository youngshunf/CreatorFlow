# Requirements Document

## Introduction

本文档定义了 CreatorFlow 视频创作 MCP Server 的需求规格。该 MCP Server 将使外部 AI Agent（如 Claude Desktop、Kiro 等）能够通过 MCP 协议调用 CreatorFlow 的视频创作功能，包括项目管理、素材管理、视频渲染和实时预览。

MCP Server 作为独立进程运行，不依赖 Electron，直接调用 Remotion API 进行视频处理。

## Glossary

- **MCP_Server**: Model Context Protocol 服务器，提供工具接口供 AI Agent 调用
- **Video_Project**: 视频项目，包含配置、组合、素材和渲染历史的完整数据结构
- **Composition**: 视频组合，定义视频内容和动画的 React 组件
- **Asset**: 素材资源，包括图片、视频、音频和字体文件
- **Render_Engine**: 渲染引擎，基于 Remotion 的视频渲染核心
- **Preview_Server**: 预览服务器，提供实时视频预览的 HTTP 服务
- **Project_Store**: 项目存储，基于文件系统的 JSON 格式项目数据持久化层，存储在工作区的用户可见目录下
- **Workspace**: 工作区，用户创建的项目工作目录
- **Workspace_Root**: 工作区根路径，用户指定的工作区目录绝对路径
- **Video_Projects_Dir**: 视频项目目录，位于 `{Workspace_Root}/视频创作/`，用户可直接查看和编辑
- **Template**: 视频模板，预定义的视频组合配置

## Requirements

### Requirement 1: 视频项目创建

**User Story:** As an AI Agent, I want to create video projects with configurable parameters, so that I can start video creation workflows programmatically.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_create_project with valid parameters, THE MCP_Server SHALL create a new Video_Project with unique ID and return the project details
2. WHEN creating a project with a template parameter, THE MCP_Server SHALL apply the template's default configuration and compositions
3. WHEN creating a project without specifying dimensions, THE MCP_Server SHALL use default values (1920x1080, 30fps)
4. IF the project name is empty or invalid, THEN THE MCP_Server SHALL return an error with descriptive message
5. WHEN a project is created, THE Project_Store SHALL persist the project data to `{Workspace_Root}/视频创作/{projectName}/project.json`
6. WHEN creating a project, THE MCP_Server SHALL require a workspacePath parameter (absolute path to Workspace_Root) to determine the storage location
7. WHEN a project is created, THE MCP_Server SHALL create a user-friendly directory structure with Chinese naming for easy browsing

### Requirement 2: 视频项目查询

**User Story:** As an AI Agent, I want to list and retrieve video projects, so that I can manage and continue working on existing projects.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_list_projects with a workspacePath, THE MCP_Server SHALL return a list of all projects in that Workspace with basic metadata (id, name, createdAt, updatedAt)
2. WHEN an AI Agent calls video_get_project with a valid project ID and workspacePath, THE MCP_Server SHALL return the complete project details including compositions and assets
3. IF the requested project ID does not exist in the Workspace, THEN THE MCP_Server SHALL return a not-found error with the requested ID
4. WHEN listing projects, THE MCP_Server SHALL sort results by updatedAt in descending order

### Requirement 3: 视频项目更新与删除

**User Story:** As an AI Agent, I want to update and delete video projects, so that I can modify project settings and clean up unused projects.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_update_project with valid parameters, THE MCP_Server SHALL update the specified fields and return the updated project
2. WHEN a project is updated, THE Project_Store SHALL update the updatedAt timestamp automatically
3. WHEN an AI Agent calls video_delete_project with a valid project ID, THE MCP_Server SHALL remove the project and all associated files
4. IF attempting to update or delete a non-existent project, THEN THE MCP_Server SHALL return a not-found error

### Requirement 4: 素材管理

**User Story:** As an AI Agent, I want to add and manage assets in video projects, so that I can include images, videos, audio, and fonts in compositions.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_add_asset with a valid file path, THE MCP_Server SHALL copy the file to `{Workspace_Root}/视频创作/{projectName}/素材/` and register it in the project
2. WHEN adding an asset, THE MCP_Server SHALL validate the file type against supported extensions (image: png/jpg/jpeg/gif/webp/svg, video: mp4/webm/mov, audio: mp3/wav/ogg/m4a, font: ttf/otf/woff/woff2)
3. IF the asset file does not exist or has unsupported format, THEN THE MCP_Server SHALL return an error with details
4. WHEN an AI Agent calls video_remove_asset, THE MCP_Server SHALL remove the asset from the project and optionally delete the file from the project 素材 directory
5. WHEN an AI Agent calls video_list_assets, THE MCP_Server SHALL return all assets in the project grouped by type

### Requirement 5: 视频渲染

**User Story:** As an AI Agent, I want to render video compositions to files, so that I can produce final video outputs in various formats.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_render with valid parameters, THE Render_Engine SHALL start rendering the specified composition
2. THE Render_Engine SHALL support output formats: MP4, WebM, and GIF
3. THE Render_Engine SHALL support quality presets: draft (fast, lower quality), standard (balanced), high (slow, best quality)
4. WHILE rendering is in progress, THE MCP_Server SHALL provide progress updates (percentage, current status)
5. WHEN rendering completes successfully, THE MCP_Server SHALL return the output file path and render statistics
6. IF rendering fails, THEN THE MCP_Server SHALL return an error with the failure reason and partial progress
7. FOR ALL valid Video_Project objects, serializing to JSON then deserializing SHALL produce an equivalent object (round-trip property)

### Requirement 6: 实时预览

**User Story:** As an AI Agent, I want to start and stop preview servers, so that users can view video compositions in real-time.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_preview_start with a valid project ID, THE Preview_Server SHALL start on an available port and return the preview URL
2. WHILE the Preview_Server is running, THE MCP_Server SHALL serve the composition with hot-reload capability
3. WHEN an AI Agent calls video_preview_stop, THE Preview_Server SHALL gracefully shutdown and release the port
4. IF the preview server fails to start, THEN THE MCP_Server SHALL return an error with the failure reason
5. THE MCP_Server SHALL track active preview servers and prevent duplicate servers for the same project

### Requirement 7: 组合管理

**User Story:** As an AI Agent, I want to add and manage compositions in video projects, so that I can define video content and animations.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_add_composition with valid parameters, THE MCP_Server SHALL add a new composition to the project
2. WHEN adding a composition, THE MCP_Server SHALL validate that the composition code references valid components
3. WHEN an AI Agent calls video_update_composition, THE MCP_Server SHALL update the composition props and code
4. WHEN an AI Agent calls video_remove_composition, THE MCP_Server SHALL remove the composition from the project
5. IF the composition ID does not exist in the project, THEN THE MCP_Server SHALL return a not-found error

### Requirement 8: MCP 协议兼容性

**User Story:** As an AI Agent developer, I want the MCP Server to follow standard MCP protocol, so that it can integrate with any MCP-compatible client.

#### Acceptance Criteria

1. THE MCP_Server SHALL support stdio transport mode for desktop MCP clients
2. THE MCP_Server SHALL support HTTP transport mode for remote/cloud deployments
3. THE MCP_Server SHALL expose tools with proper JSON Schema definitions for input validation
4. THE MCP_Server SHALL return structured JSON responses with consistent error format
5. WHEN the MCP_Server starts, THE MCP_Server SHALL log available tools and transport configuration

### Requirement 9: 错误处理与日志

**User Story:** As an AI Agent, I want clear error messages and logging, so that I can diagnose issues and handle failures gracefully.

#### Acceptance Criteria

1. WHEN any operation fails, THE MCP_Server SHALL return a structured error with code, message, and context
2. THE MCP_Server SHALL log all tool invocations with timestamps and parameters
3. THE MCP_Server SHALL log render progress and completion status
4. IF an unexpected error occurs, THEN THE MCP_Server SHALL catch it and return a safe error response without exposing internal details

### Requirement 10: 模板支持

**User Story:** As an AI Agent, I want to use predefined templates, so that I can quickly create professional-looking videos.

#### Acceptance Criteria

1. WHEN an AI Agent calls video_list_templates, THE MCP_Server SHALL return all available templates with metadata
2. WHEN creating a project with a template, THE MCP_Server SHALL apply the template's default config, compositions, and props
3. THE MCP_Server SHALL support templates from @creator-flow/video package (social-media, marketing, tutorial categories)
