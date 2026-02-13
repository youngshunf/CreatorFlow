# Video MCP Source 部署完成

## 完成时间
2026-02-12

## 部署位置

**仓库**：`skill-marketplace`
**路径**：`apps/app-creator-media/sources/video-mcp/`

## 文件清单

### 1. config.json
MCP source 配置文件

```json
{
  "id": "video-mcp_preset",
  "name": "视频创作引擎",
  "slug": "video-mcp",
  "enabled": true,
  "provider": "sprouty-video",
  "type": "mcp",
  "mcp": {
    "command": "bun",
    "args": ["run", "src/mcp-server/index.ts"],
    "cwd": "app:resources/video",
    "transport": "stdio",
    "env": {}
  },
  "icon": "🎬",
  "tagline": "Remotion 视频创作 MCP 服务 — 模板管理、合成渲染、素材发现、项目管理",
  "connectionStatus": "untested"
}
```

**关键配置**：
- `cwd: "app:resources/video"` - 使用 app: 前缀指向打包后的位置
- `command: "bun"` - 系统会自动探测 bun 的安装路径
- `transport: "stdio"` - 使用标准输入输出通信

### 2. guide.md
完整的使用指南，包含：

- **Scope**：服务能力范围
- **可用工具 (18)**：按功能分类的工具列表
  - 项目管理（5 个）
  - 素材管理（4 个）
  - 组合管理（3 个）
  - 模板（2 个）
  - 渲染与预览（3 个）
  - 代码校验（1 个）
- **Guidelines**：基本用法和注意事项
- **支持的模板类型**：4 种模板详细说明
- **Context**：服务架构和技术栈
- **配置说明**：cwd 路径解析规则和命令路径探测
- **故障排除**：常见问题和解决方案

## 技术实现

### 路径解析机制

**开发环境**：
```json
"cwd": "../../packages/video"
```
解析为：`/Users/mac/saas/creator-flow/Sprouty AI/packages/video`

**生产环境（打包后）**：
```json
"cwd": "app:resources/video"
```
解析为：`/Applications/Sprouty AI.app/Contents/Resources/resources/video`

### 命令路径探测

系统自动探测 bun 命令，搜索顺序：
1. 使用 `which` 命令（扩展 PATH）
2. 检查 well-known 安装位置：
   - `~/.bun/bin/bun`
   - `/opt/homebrew/bin/bun`
   - `/usr/local/bin/bun`
3. 回退到原始命令名

### 打包资源复制

`apps/electron/scripts/copy-assets.ts` 会在打包时：
1. 复制 `packages/video/` 到 `dist/resources/video/`
2. 自动过滤不必要的文件：
   - node_modules
   - 测试文件（*.test.*）
   - 构建产物（dist/、.remotion/）

## 验证清单

- [x] config.json 配置正确
- [x] guide.md 文档完整
- [x] cwd 使用 app: 前缀
- [x] 命令路径自动探测已实现
- [x] 打包脚本包含 video 包复制
- [x] 所有 192 个测试通过
- [ ] 在打包应用中测试连接（待验证）
- [ ] 验证 18 个工具可用性（待验证）

## 下一步

1. **构建应用**
   ```bash
   cd apps/electron
   bun run build
   ```

2. **测试连接**
   - 启动打包后的应用
   - 在 app-creator-media 中启用 video-mcp source
   - 验证连接状态
   - 测试工具调用

3. **跨平台测试**
   - macOS（已配置）
   - Windows（需要测试）
   - Linux（需要测试）

## 相关文档

- [MCP 打包支持总结](./MCP_PACKAGING_SUMMARY.md)
- [Video MCP 使用指南](./VIDEO_MCP_GUIDE.md)
- [优化总结](./OPTIMIZATION_SUMMARY.md)
- [变更日志](./CHANGELOG.md)

## 注意事项

1. **开发环境 vs 生产环境**
   - 开发时使用相对路径（`../../packages/video`）
   - 打包后使用 app: 前缀（`app:resources/video`）

2. **命令路径**
   - 优先使用命令名（`bun`），让系统自动探测
   - 仅在自动探测失败时使用绝对路径

3. **资源复制**
   - 确保 `copy-assets.ts` 在打包前执行
   - 验证 `dist/resources/video/` 目录存在

4. **故障排除**
   - 检查应用日志中的 MCP 连接错误
   - 使用 `which bun` 验证命令可用性
   - 确认 video 包已正确复制到应用资源目录
