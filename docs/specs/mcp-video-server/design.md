# Design Document: MCP Video Server

## Overview

MCP Video Server 是一个独立的 MCP (Model Context Protocol) 服务器，为 AI Agent 提供视频创作能力。它基于 Bun 运行时和 FastMCP 框架构建，直接调用 Remotion API 进行视频渲染和预览。

### 设计目标

1. **独立运行**: 不依赖 Electron，可作为独立进程运行
2. **MCP 兼容**: 支持 stdio 和 HTTP 两种传输模式
3. **用户友好**: 项目文件使用中文命名，用户可直接查看和编辑
4. **类型安全**: 使用 TypeScript 和 Zod 进行类型验证
5. **复用现有代码**: 复用 @sprouty-ai/video 包的类型定义和模板

## Architecture

```mermaid
graph TB
    subgraph "External"
        Agent[AI Agent<br/>Claude Desktop / Kiro]
    end
    
    subgraph "MCP Video Server"
        MCP[FastMCP Server]
        Tools[Tool Handlers]
        ProjectStore[Project Store]
        RenderEngine[Render Engine]
        PreviewServer[Preview Server]
    end
    
    subgraph "Storage"
        FS[(File System<br/>视频创作/)]
    end
    
    subgraph "Dependencies"
        Remotion[Remotion API]
        VideoPackage[@sprouty-ai/video]
    end
    
    Agent -->|MCP Protocol| MCP
    MCP --> Tools
    Tools --> ProjectStore
    Tools --> RenderEngine
    Tools --> PreviewServer
    ProjectStore --> FS
    RenderEngine --> Remotion
    PreviewServer --> Remotion
    Tools --> VideoPackage
```

### 目录结构

```
Sprouty AI/apps/mcp-video/
├── src/
│   ├── index.ts              # 入口文件
│   ├── server.ts             # FastMCP 服务器配置
│   ├── tools/                # MCP 工具实现
│   │   ├── index.ts
│   │   ├── project.ts        # 项目管理工具
│   │   ├── asset.ts          # 素材管理工具
│   │   ├── render.ts         # 渲染工具
│   │   ├── preview.ts        # 预览工具
│   │   └── composition.ts    # 组合管理工具
│   ├── services/             # 核心服务
│   │   ├── project-store.ts  # 项目存储服务
│   │   ├── render-engine.ts  # 渲染引擎
│   │   └── preview-server.ts # 预览服务器
│   ├── types/                # 类型定义
│   │   ├── index.ts
│   │   └── errors.ts         # 错误类型
│   └── utils/                # 工具函数
│       ├── paths.ts          # 路径处理
│       └── validators.ts     # 验证器
├── package.json
├── tsconfig.json
└── README.md
```

## Components and Interfaces

### 1. FastMCP Server (server.ts)

```typescript
import { FastMCP } from 'fastmcp';

interface ServerConfig {
  transport: 'stdio' | 'http';
  host?: string;
  port?: number;
}

// FastMCP 服务器实例
const mcp = FastMCP('sprouty-ai-video');

// 启动服务器
function runServer(config: ServerConfig): void;
```

### 2. Project Store Service (services/project-store.ts)

```typescript
interface ProjectStoreConfig {
  workspacePath: string;
}

class ProjectStore {
  constructor(config: ProjectStoreConfig);
  
  // 项目 CRUD
  createProject(input: CreateProjectInput): Promise<VideoProject>;
  getProject(projectId: string): Promise<VideoProject | null>;
  listProjects(): Promise<ProjectSummary[]>;
  updateProject(projectId: string, updates: Partial<VideoProject>): Promise<VideoProject>;
  deleteProject(projectId: string): Promise<boolean>;
  
  // 路径工具
  getProjectPath(projectName: string): string;
  getAssetsPath(projectName: string): string;
  
  // 私有方法
  private loadProject(projectPath: string): VideoProject | null;
  private saveProject(project: VideoProject): void;
}

interface CreateProjectInput {
  name: string;
  template?: string;
  config?: Partial<VideoConfig>;
  description?: string;
}

interface ProjectSummary {
  id: string;
  name: string;
  createdAt: string;
  updatedAt: string;
  compositionCount: number;
  assetCount: number;
}
```

### 3. Render Engine Service (services/render-engine.ts)

```typescript
import { bundle, renderMedia } from '@remotion/bundler';
import type { RenderMediaOnProgress } from '@remotion/renderer';

interface RenderConfig {
  projectPath: string;
  compositionId: string;
  outputFormat: 'mp4' | 'webm' | 'gif';
  quality: 'draft' | 'standard' | 'high';
  outputPath?: string;
}

interface RenderResult {
  success: boolean;
  outputPath: string;
  duration: number;
  fileSize: number;
}

interface RenderProgress {
  status: 'bundling' | 'preparing' | 'rendering' | 'completed' | 'failed';
  progress: number;
  error?: string;
}

class RenderEngine {
  // 渲染视频
  async render(
    config: RenderConfig,
    onProgress?: (progress: RenderProgress) => void
  ): Promise<RenderResult>;
  
  // 获取质量预设配置
  private getQualityConfig(quality: QualityPreset): QualityPresetConfig;
  
  // 构建 Remotion bundle
  private async bundleProject(projectPath: string): Promise<string>;
}
```

### 4. Preview Server Service (services/preview-server.ts)

```typescript
interface PreviewConfig {
  projectPath: string;
  compositionId?: string;
  port?: number;
}

interface PreviewInstance {
  projectId: string;
  port: number;
  url: string;
  pid: number;
}

class PreviewServerManager {
  private activeServers: Map<string, PreviewInstance>;
  
  // 启动预览服务器
  async start(config: PreviewConfig): Promise<PreviewInstance>;
  
  // 停止预览服务器
  async stop(projectId: string): Promise<boolean>;
  
  // 获取活跃的预览服务器
  getActiveServer(projectId: string): PreviewInstance | null;
  
  // 列出所有活跃服务器
  listActiveServers(): PreviewInstance[];
  
  // 查找可用端口
  private findAvailablePort(startPort: number): Promise<number>;
}
```

### 5. MCP Tool Handlers

#### 5.1 Project Tools (tools/project.ts)

```typescript
// video_create_project
const createProjectTool = {
  name: 'video_create_project',
  description: '创建新的视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    name: z.string().describe('项目名称'),
    template: z.string().optional().describe('模板ID'),
    width: z.number().optional().default(1920),
    height: z.number().optional().default(1080),
    fps: z.number().optional().default(30),
    durationInSeconds: z.number().describe('视频时长（秒）'),
    description: z.string().optional(),
  }),
};

// video_list_projects
const listProjectsTool = {
  name: 'video_list_projects',
  description: '列出工作区中的所有视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
  }),
};

// video_get_project
const getProjectTool = {
  name: 'video_get_project',
  description: '获取视频项目详情',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
  }),
};

// video_update_project
const updateProjectTool = {
  name: 'video_update_project',
  description: '更新视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    name: z.string().optional(),
    description: z.string().optional(),
    config: z.object({
      width: z.number().optional(),
      height: z.number().optional(),
      fps: z.number().optional(),
      durationInFrames: z.number().optional(),
    }).optional(),
  }),
};

// video_delete_project
const deleteProjectTool = {
  name: 'video_delete_project',
  description: '删除视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
  }),
};
```

#### 5.2 Asset Tools (tools/asset.ts)

```typescript
// video_add_asset
const addAssetTool = {
  name: 'video_add_asset',
  description: '添加素材到视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    assetPath: z.string().describe('素材文件路径'),
    assetType: z.enum(['image', 'video', 'audio', 'font']).describe('素材类型'),
    name: z.string().optional().describe('素材名称'),
  }),
};

// video_remove_asset
const removeAssetTool = {
  name: 'video_remove_asset',
  description: '从视频项目移除素材',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    assetId: z.string().describe('素材ID'),
    deleteFile: z.boolean().optional().default(false).describe('是否删除文件'),
  }),
};

// video_list_assets
const listAssetsTool = {
  name: 'video_list_assets',
  description: '列出项目中的所有素材',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
  }),
};
```

#### 5.3 Render Tools (tools/render.ts)

```typescript
// video_render
const renderTool = {
  name: 'video_render',
  description: '渲染视频到文件',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    compositionId: z.string().describe('组合ID'),
    outputFormat: z.enum(['mp4', 'webm', 'gif']).optional().default('mp4'),
    quality: z.enum(['draft', 'standard', 'high']).optional().default('standard'),
    outputPath: z.string().optional().describe('输出文件路径'),
  }),
};
```

#### 5.4 Preview Tools (tools/preview.ts)

```typescript
// video_preview_start
const previewStartTool = {
  name: 'video_preview_start',
  description: '启动视频预览服务器',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    compositionId: z.string().optional().describe('组合ID'),
  }),
};

// video_preview_stop
const previewStopTool = {
  name: 'video_preview_stop',
  description: '停止视频预览服务器',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
  }),
};
```

#### 5.5 Composition Tools (tools/composition.ts)

```typescript
// video_add_composition
const addCompositionTool = {
  name: 'video_add_composition',
  description: '添加组合到视频项目',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    name: z.string().describe('组合名称'),
    code: z.string().describe('组合代码（React 组件）'),
    props: z.record(z.any()).optional().describe('组合属性'),
  }),
};

// video_update_composition
const updateCompositionTool = {
  name: 'video_update_composition',
  description: '更新视频组合',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    compositionId: z.string().describe('组合ID'),
    name: z.string().optional(),
    code: z.string().optional(),
    props: z.record(z.any()).optional(),
  }),
};

// video_remove_composition
const removeCompositionTool = {
  name: 'video_remove_composition',
  description: '从视频项目移除组合',
  inputSchema: z.object({
    workspacePath: z.string().describe('工作区根路径'),
    projectId: z.string().describe('项目ID'),
    compositionId: z.string().describe('组合ID'),
  }),
};

// video_list_templates
const listTemplatesTool = {
  name: 'video_list_templates',
  description: '列出所有可用的视频模板',
  inputSchema: z.object({}),
};
```

## Data Models

### VideoProject (复用 @sprouty-ai/video)

```typescript
interface VideoProject {
  id: string;
  name: string;
  description?: string;
  createdAt: string;
  updatedAt: string;
  config: VideoConfig;
  compositions: Composition[];
  assets: Asset[];
  renders: RenderHistory[];
}

interface VideoConfig {
  width: number;      // 默认 1920
  height: number;     // 默认 1080
  fps: number;        // 默认 30
  durationInFrames: number;
}

interface Composition {
  id: string;
  name: string;
  code: string;       // React 组件代码或文件路径
  props: Record<string, any>;
}

interface Asset {
  id: string;
  type: 'image' | 'video' | 'audio' | 'font';
  name: string;
  path: string;       // 相对于项目目录的路径
}

interface RenderHistory {
  id: string;
  compositionId: string;
  outputPath: string;
  status: 'pending' | 'rendering' | 'completed' | 'failed';
  progress: number;
  createdAt: string;
  error?: string;
}
```

### 文件系统结构

```
{Workspace_Root}/
└── 视频创作/
    └── {项目名称}/
        ├── project.json          # 项目配置
        ├── 素材/                  # 素材文件
        │   ├── images/
        │   ├── videos/
        │   ├── audio/
        │   └── fonts/
        ├── 组合/                  # 组合代码
        │   └── {compositionId}.tsx
        └── 输出/                  # 渲染输出
            └── {renderName}.mp4
```

### MCP 响应格式

```typescript
// 成功响应
interface SuccessResponse<T> {
  success: true;
  data: T;
}

// 错误响应
interface ErrorResponse {
  success: false;
  error: {
    code: string;
    message: string;
    details?: Record<string, any>;
  };
}

// 错误代码
type ErrorCode =
  | 'PROJECT_NOT_FOUND'
  | 'ASSET_NOT_FOUND'
  | 'COMPOSITION_NOT_FOUND'
  | 'INVALID_INPUT'
  | 'FILE_NOT_FOUND'
  | 'UNSUPPORTED_FORMAT'
  | 'RENDER_FAILED'
  | 'PREVIEW_FAILED'
  | 'INTERNAL_ERROR';
```



## Correctness Properties

*A property is a characteristic or behavior that should hold true across all valid executions of a system—essentially, a formal statement about what the system should do. Properties serve as the bridge between human-readable specifications and machine-verifiable correctness guarantees.*

### Property 1: Project Creation Produces Valid Persisted Project

*For any* valid project creation input (name, workspacePath, optional template and config), creating a project SHALL result in:
- A project with a unique ID being returned
- The project data being persisted to `{workspacePath}/视频创作/{projectName}/project.json`
- The directory structure including `素材/` subdirectory being created

**Validates: Requirements 1.1, 1.5, 1.7**

### Property 2: Project Query Operations Return Correct Data

*For any* workspace with N projects, listing projects SHALL return exactly N projects sorted by updatedAt descending, and getting any project by ID SHALL return the complete project data matching what was created/updated.

**Validates: Requirements 2.1, 2.2, 2.4**

### Property 3: Project Mutations Work Correctly

*For any* existing project, updating it SHALL modify only the specified fields and automatically update the updatedAt timestamp, and deleting it SHALL remove the project from both memory and file system.

**Validates: Requirements 3.1, 3.2, 3.3**

### Property 4: Asset Management Operations Work Correctly

*For any* valid asset file with supported extension, adding it to a project SHALL copy the file to the project's `素材/` directory and register it in the project. Removing an asset SHALL unregister it and optionally delete the file. Listing assets SHALL return all assets grouped by type.

**Validates: Requirements 4.1, 4.2, 4.4, 4.5**

### Property 5: Composition Management Operations Work Correctly

*For any* valid composition definition, adding it to a project SHALL register it with a unique ID. Updating a composition SHALL modify the specified fields. Removing a composition SHALL unregister it from the project.

**Validates: Requirements 7.1, 7.3, 7.4**

### Property 6: Project Serialization Round-Trip

*For any* valid VideoProject object, serializing it to JSON and then deserializing SHALL produce an equivalent object with all fields preserved.

**Validates: Requirements 5.7**

### Property 7: Error Handling Returns Structured Errors

*For any* invalid input (empty project name, non-existent file, unsupported format, unexpected error), the MCP_Server SHALL return a structured error response with code, message, and optional details, never exposing internal implementation details.

**Validates: Requirements 1.4, 4.3, 9.1, 9.4**

### Property 8: MCP Protocol Compliance

*For any* registered MCP tool, the tool SHALL have a valid JSON Schema definition for its input, and all responses SHALL follow the consistent success/error response format.

**Validates: Requirements 8.3, 8.4**

### Property 9: Preview Server Duplicate Prevention

*For any* project with an active preview server, attempting to start another preview server for the same project SHALL either return the existing server's URL or prevent the duplicate.

**Validates: Requirements 6.5**

## Error Handling

### Error Categories

1. **Validation Errors** (4xx equivalent)
   - `INVALID_INPUT`: Missing or malformed parameters
   - `PROJECT_NOT_FOUND`: Requested project does not exist
   - `ASSET_NOT_FOUND`: Requested asset does not exist
   - `COMPOSITION_NOT_FOUND`: Requested composition does not exist
   - `FILE_NOT_FOUND`: Source file does not exist
   - `UNSUPPORTED_FORMAT`: File format not supported

2. **Operation Errors** (5xx equivalent)
   - `RENDER_FAILED`: Video rendering failed
   - `PREVIEW_FAILED`: Preview server failed to start
   - `INTERNAL_ERROR`: Unexpected internal error

### Error Response Format

```typescript
interface MCPError {
  code: ErrorCode;
  message: string;        // Human-readable message
  details?: {
    field?: string;       // Field that caused the error
    expected?: string;    // Expected value/format
    received?: string;    // Actual value received
    path?: string;        // File path if relevant
  };
}
```

### Error Handling Strategy

1. **Input Validation**: Use Zod schemas to validate all inputs before processing
2. **File Operations**: Wrap all file operations in try-catch, return specific error codes
3. **Remotion Operations**: Catch Remotion errors and translate to MCP error format
4. **Unexpected Errors**: Catch all unhandled errors, log details, return safe INTERNAL_ERROR

## Testing Strategy

### Dual Testing Approach

本项目采用单元测试和属性测试相结合的方式：

- **单元测试**: 验证具体示例、边界情况和错误条件
- **属性测试**: 验证跨所有输入的通用属性

### Property-Based Testing Configuration

- **测试框架**: fast-check (TypeScript 属性测试库)
- **最小迭代次数**: 每个属性测试 100 次
- **标签格式**: `Feature: mcp-video-server, Property {number}: {property_text}`

### Test Categories

1. **Unit Tests**
   - Project CRUD operations with specific examples
   - Asset type validation with known extensions
   - Error response format verification
   - Path generation with edge cases (special characters, long names)

2. **Property Tests**
   - Property 1: Project creation invariants
   - Property 2: Query operation correctness
   - Property 3: Mutation operation correctness
   - Property 4: Asset management invariants
   - Property 5: Composition management invariants
   - Property 6: Serialization round-trip
   - Property 7: Error handling consistency
   - Property 8: MCP schema validity
   - Property 9: Preview server uniqueness

3. **Integration Tests**
   - Remotion bundle and render (manual/CI)
   - Preview server lifecycle (manual/CI)
   - MCP transport modes (stdio/HTTP)

### Test File Structure

```
Sprouty AI/apps/mcp-video/
└── tests/
    ├── unit/
    │   ├── project-store.test.ts
    │   ├── asset-validator.test.ts
    │   └── path-utils.test.ts
    ├── property/
    │   ├── project.property.test.ts
    │   ├── asset.property.test.ts
    │   ├── composition.property.test.ts
    │   └── serialization.property.test.ts
    └── integration/
        ├── render.integration.test.ts
        └── preview.integration.test.ts
```
