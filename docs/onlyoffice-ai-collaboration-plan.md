# OnlyOffice 云端集成 + AI 协作技术方案

## 概述

本方案采用云端 OnlyOffice Document Server + 桌面端 Electron 客户端架构，实现 Office 文档的高保真编辑与 AI 协作功能。

### 核心优势

- **零本地依赖**：用户无需安装 LibreOffice 或 OnlyOffice
- **高保真编辑**：支持复杂样式、批注、修订等高级功能
- **统一体验**：跨平台一致的编辑体验
- **易于维护**：服务端统一升级，客户端轻量化

### 架构图

```
┌─────────────────────────────────────────────────────────────┐
│                    CreatorFlow Electron                      │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │ File Manager│  │   webview   │  │    AI Service       │  │
│  │  (Chonky)   │  │ OnlyOffice  │  │  (CraftAgent)       │  │
│  │             │  │   Editor    │  │                     │  │
│  └──────┬──────┘  └──────┬──────┘  └──────────┬──────────┘  │
│         │                │                     │             │
│         └────────────────┼─────────────────────┘             │
│                          │ postMessage                       │
│                    ┌─────┴─────┐                             │
│                    │  IPC Bus  │                             │
│                    └─────┬─────┘                             │
└──────────────────────────┼──────────────────────────────────┘
                           │ HTTPS
          ┌────────────────┼────────────────┐
          │                │                │
          ▼                ▼                ▼
┌─────────────────┐ ┌─────────────┐ ┌─────────────────┐
│  DMS API        │ │ OnlyOffice  │ │  File Storage   │
│  (文件回调服务)  │ │ Doc Server  │ │  (S3/OSS/本地)  │
└─────────────────┘ └─────────────┘ └─────────────────┘
```

---

## Phase 1: 云端 OnlyOffice 部署（2 天）

### 1.1 Docker Compose 部署

```yaml
# docker-compose.yml
version: '3.8'

services:
  onlyoffice:
    image: onlyoffice/documentserver:latest
    container_name: onlyoffice-docs
    ports:
      - "8080:80"
      - "8443:443"
    environment:
      - JWT_ENABLED=true
      - JWT_SECRET=${ONLYOFFICE_JWT_SECRET}
      - JWT_HEADER=Authorization
      - JWT_IN_BODY=false
    volumes:
      - oo_data:/var/www/onlyoffice/Data
      - oo_logs:/var/log/onlyoffice
      - oo_fonts:/usr/share/fonts/truetype/custom
    restart: unless-stopped

volumes:
  oo_data:
  oo_logs:
  oo_fonts:
```

### 1.2 Nginx 反向代理配置

```nginx
# /etc/nginx/conf.d/onlyoffice.conf
upstream onlyoffice {
    server 127.0.0.1:8080;
}

server {
    listen 443 ssl http2;
    server_name docs.creatorflow.app;

    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;

    location / {
        proxy_pass http://onlyoffice;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # WebSocket 支持
        proxy_read_timeout 3600s;
        proxy_send_timeout 3600s;
    }
}
```

### 1.3 健康检查端点

```bash
# 验证部署成功
curl https://docs.creatorflow.app/healthcheck
# 预期返回: true
```

---

## Phase 2: DMS 文件服务适配（3 天）

### 2.1 DMS API 接口设计

OnlyOffice 需要通过回调 URL 读取和保存文件，我们需要实现以下接口：

```typescript
// 接口定义
interface DMSApi {
  // 获取文件内容 (OnlyOffice 调用)
  'GET /api/dms/files/:fileId': () => FileStream;
  
  // 保存文件内容 (OnlyOffice 回调)
  'POST /api/dms/callback': (body: OnlyOfficeCallback) => { error: 0 };
  
  // 获取文件编辑配置 (Electron 调用)
  'POST /api/dms/editor-config': (body: EditorConfigRequest) => EditorConfig;
}
```

### 2.2 文件上传与存储

```typescript
// server/routes/dms.ts
import { Router } from 'express';
import { S3Client, PutObjectCommand, GetObjectCommand } from '@aws-sdk/client-s3';
import { v4 as uuidv4 } from 'uuid';

const router = Router();
const s3 = new S3Client({ region: process.env.AWS_REGION });

// 上传文件
router.post('/files/upload', upload.single('file'), async (req, res) => {
  const fileId = uuidv4();
  const key = `documents/${fileId}/${req.file.originalname}`;
  
  await s3.send(new PutObjectCommand({
    Bucket: process.env.S3_BUCKET,
    Key: key,
    Body: req.file.buffer,
    ContentType: req.file.mimetype
  }));
  
  // 记录文件元数据
  await db.files.create({
    id: fileId,
    name: req.file.originalname,
    key,
    size: req.file.size,
    workspaceId: req.body.workspaceId,
    createdBy: req.user.id
  });
  
  res.json({ fileId, name: req.file.originalname });
});

// 获取文件内容 (供 OnlyOffice 拉取)
router.get('/files/:fileId', async (req, res) => {
  const file = await db.files.findById(req.params.fileId);
  if (!file) return res.status(404).send('File not found');
  
  const object = await s3.send(new GetObjectCommand({
    Bucket: process.env.S3_BUCKET,
    Key: file.key
  }));
  
  res.setHeader('Content-Type', file.mimeType);
  res.setHeader('Content-Disposition', `attachment; filename="${file.name}"`);
  object.Body.pipe(res);
});
```

### 2.3 OnlyOffice 回调处理

```typescript
// server/routes/dms.ts
interface OnlyOfficeCallback {
  key: string;           // 文档唯一标识
  status: number;        // 2=ready to save, 6=force save
  url?: string;          // 修改后文档的下载 URL
  changesurl?: string;   // 变更历史 URL
  history?: object;      // 版本历史
  users?: string[];      // 当前编辑用户
}

router.post('/callback', async (req, res) => {
  const { key, status, url } = req.body as OnlyOfficeCallback;
  
  // 验证 JWT
  const token = req.headers.authorization?.replace('Bearer ', '');
  if (!verifyJWT(token)) {
    return res.status(403).json({ error: 1 });
  }
  
  // status 2 或 6 表示需要保存
  if (status === 2 || status === 6) {
    try {
      // 从 OnlyOffice 下载修改后的文件
      const response = await fetch(url);
      const buffer = await response.arrayBuffer();
      
      // 保存到存储
      const file = await db.files.findByKey(key);
      await s3.send(new PutObjectCommand({
        Bucket: process.env.S3_BUCKET,
        Key: file.key,
        Body: Buffer.from(buffer)
      }));
      
      // 更新元数据
      await db.files.update(file.id, {
        updatedAt: new Date(),
        version: file.version + 1
      });
      
      res.json({ error: 0 });
    } catch (err) {
      console.error('Save failed:', err);
      res.json({ error: 1 });
    }
  } else {
    res.json({ error: 0 });
  }
});
```

### 2.4 JWT 签名服务

```typescript
// server/services/jwt.ts
import jwt from 'jsonwebtoken';

const JWT_SECRET = process.env.ONLYOFFICE_JWT_SECRET;

export function signEditorConfig(config: object): string {
  return jwt.sign(config, JWT_SECRET, { expiresIn: '1h' });
}

export function verifyJWT(token: string): boolean {
  try {
    jwt.verify(token, JWT_SECRET);
    return true;
  } catch {
    return false;
  }
}
```

---

## Phase 3: Electron 客户端集成（3 天）

### 3.1 OnlyOffice 配置生成

```typescript
// apps/electron/src/main/onlyoffice/config.ts
interface EditorConfig {
  document: {
    fileType: string;
    key: string;
    title: string;
    url: string;
    permissions: {
      edit: boolean;
      download: boolean;
      print: boolean;
    };
  };
  editorConfig: {
    callbackUrl: string;
    lang: string;
    mode: 'edit' | 'view';
    user: {
      id: string;
      name: string;
    };
    customization: {
      autosave: boolean;
      forcesave: boolean;
      plugins: boolean;
    };
  };
  token?: string;
}

export async function generateEditorConfig(
  filePath: string,
  userId: string,
  userName: string
): Promise<EditorConfig> {
  // 1. 上传本地文件到 DMS
  const fileInfo = await uploadFileToDMS(filePath);
  
  // 2. 生成唯一文档 key
  const documentKey = `${fileInfo.fileId}_${Date.now()}`;
  
  // 3. 构建配置
  const config: EditorConfig = {
    document: {
      fileType: getFileExtension(filePath),
      key: documentKey,
      title: path.basename(filePath),
      url: `${DMS_BASE_URL}/api/dms/files/${fileInfo.fileId}`,
      permissions: {
        edit: true,
        download: true,
        print: true
      }
    },
    editorConfig: {
      callbackUrl: `${DMS_BASE_URL}/api/dms/callback`,
      lang: 'zh',
      mode: 'edit',
      user: {
        id: userId,
        name: userName
      },
      customization: {
        autosave: true,
        forcesave: true,
        plugins: true
      }
    }
  };
  
  // 4. JWT 签名
  const response = await fetch(`${DMS_BASE_URL}/api/dms/editor-config`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(config)
  });
  
  return response.json();
}
```

### 3.2 IPC 通道定义

```typescript
// apps/electron/src/shared/types.ts
export type IpcChannels = {
  // ... 现有 channels
  
  // OnlyOffice 相关
  'onlyoffice:getConfig': (filePath: string) => Promise<EditorConfig>;
  'onlyoffice:processAI': (params: AIProcessParams) => Promise<string>;
  'onlyoffice:downloadFile': (fileId: string, savePath: string) => Promise<void>;
};

interface AIProcessParams {
  sessionId: string;
  action: string;
  text: string;
  customPrompt?: string;
}
```

### 3.3 主进程 IPC 处理

```typescript
// apps/electron/src/main/onlyoffice/ipc.ts
import { ipcMain } from 'electron';
import { generateEditorConfig } from './config';
import { handleOnlyOfficeAI } from './ai-handler';

export function registerOnlyOfficeIpc() {
  ipcMain.handle('onlyoffice:getConfig', async (_, filePath: string) => {
    const user = await getCurrentUser();
    return generateEditorConfig(filePath, user.id, user.name);
  });
  
  ipcMain.handle('onlyoffice:processAI', async (_, params: AIProcessParams) => {
    const agent = await getSessionAgent(params.sessionId);
    return handleOnlyOfficeAI(agent, params.action, params.text, params.customPrompt);
  });
}
```

### 3.4 Renderer 编辑器组件

```typescript
// apps/electron/src/renderer/components/onlyoffice/OnlyOfficeEditor.tsx
import { useEffect, useRef, useState } from 'react';
import { Loader2 } from 'lucide-react';

interface OnlyOfficeEditorProps {
  filePath: string;
  sessionId: string;
}

export function OnlyOfficeEditor({ filePath, sessionId }: OnlyOfficeEditorProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let docEditor: any = null;

    const initEditor = async () => {
      try {
        // 1. 获取编辑器配置
        const config = await window.electronAPI.onlyoffice.getConfig(filePath);
        
        // 2. 添加事件处理
        config.events = {
          onDocumentReady: () => setLoading(false),
          onError: (e: any) => setError(e.data?.errorDescription || '编辑器错误'),
          onRequestClose: () => window.close()
        };
        
        // 3. 初始化编辑器
        docEditor = new (window as any).DocsAPI.DocEditor(
          containerRef.current!.id,
          config
        );
      } catch (err) {
        setError(err instanceof Error ? err.message : '初始化失败');
      }
    };

    initEditor();

    return () => {
      docEditor?.destroyEditor();
    };
  }, [filePath]);

  if (error) {
    return (
      <div className="flex items-center justify-center h-full">
        <p className="text-destructive">{error}</p>
      </div>
    );
  }

  return (
    <div className="relative w-full h-full">
      {loading && (
        <div className="absolute inset-0 flex items-center justify-center bg-background/80">
          <Loader2 className="w-8 h-8 animate-spin" />
        </div>
      )}
      <div
        id="onlyoffice-editor"
        ref={containerRef}
        className="w-full h-full"
      />
    </div>
  );
}
```

### 3.5 加载 OnlyOffice API 脚本

```typescript
// apps/electron/src/renderer/hooks/useOnlyOfficeScript.ts
import { useEffect, useState } from 'react';

const ONLYOFFICE_URL = 'https://docs.creatorflow.app';

export function useOnlyOfficeScript() {
  const [loaded, setLoaded] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  useEffect(() => {
    if ((window as any).DocsAPI) {
      setLoaded(true);
      return;
    }

    const script = document.createElement('script');
    script.src = `${ONLYOFFICE_URL}/web-apps/apps/api/documents/api.js`;
    script.async = true;
    script.onload = () => setLoaded(true);
    script.onerror = () => setError(new Error('Failed to load OnlyOffice API'));
    
    document.body.appendChild(script);

    return () => {
      document.body.removeChild(script);
    };
  }, []);

  return { loaded, error };
}
```

---

## Phase 4: AI 协作插件开发（5 天）

### 4.1 插件架构设计

```
onlyoffice-ai-plugin/
├── config.json              # 插件描述与入口配置
├── index.html               # 插件 UI 面板
├── scripts/
│   ├── plugin.js            # 主逻辑：与编辑器 API 交互
│   ├── ai-bridge.js         # 与 Electron postMessage 桥接
│   └── ui-controller.js     # 面板 UI 交互逻辑
├── styles/
│   └── plugin.css
└── translations/
    ├── en.json
    └── zh.json
```

### 4.2 config.json 示例

```json
{
  "name": "CreatorFlow AI",
  "guid": "asc.{CF-AI-PLUGIN-GUID}",
  "version": "1.0.0",
  "variations": [
    {
      "description": "AI 协作助手",
      "url": "index.html",
      "icons": ["icon.png"],
      "isViewer": false,
      "EditorsSupport": ["word", "cell", "slide"],
      "isVisual": true,
      "isModal": false,
      "isInsideMode": true,
      "initDataType": "text",
      "initData": "",
      "buttons": [
        { "text": "AI 助手", "primary": true }
      ]
    }
  ]
}
```

### 4.3 核心 AI 功能实现

#### 4.3.1 获取选区文本

```javascript
// scripts/plugin.js
window.Asc.plugin.init = function() {
  // 插件初始化
};

window.Asc.plugin.button = function(id) {
  if (id === 0) { // 主按钮点击
    this.callCommand(function() {
      var selection = Api.GetDocument().GetSelectedText();
      return selection;
    }, false, true, function(result) {
      // 发送到 AI 服务
      sendToElectron('ai:process', { text: result, action: 'current' });
    });
  }
};
```

#### 4.3.2 Electron 桥接通信

```javascript
// scripts/ai-bridge.js
function sendToElectron(channel, data) {
  window.parent.postMessage({
    type: 'onlyoffice-ai',
    channel: channel,
    payload: data
  }, '*');
}

window.addEventListener('message', function(event) {
  if (event.data.type === 'ai-response') {
    handleAIResponse(event.data.payload);
  }
});

function handleAIResponse(payload) {
  switch (payload.action) {
    case 'replace':
      replaceSelection(payload.text);
      break;
    case 'insert':
      insertAtCursor(payload.text);
      break;
    case 'comment':
      addComment(payload.text);
      break;
  }
}
```

#### 4.3.3 替换选区内容

```javascript
function replaceSelection(newText) {
  window.Asc.plugin.callCommand(function() {
    var doc = Api.GetDocument();
    var range = doc.GetRangeBySelect();
    if (range) {
      range.Delete();
      range.AddText(newText);
    }
  }, true, true);
}
```

### 4.4 AI 功能菜单

```javascript
// scripts/ui-controller.js
const AI_ACTIONS = [
  { id: 'translate', label: '翻译', icon: '🌐', prompt: 'translate' },
  { id: 'rewrite', label: '改写润色', icon: '✨', prompt: 'rewrite' },
  { id: 'expand', label: '扩展内容', icon: '📝', prompt: 'expand' },
  { id: 'summarize', label: '总结摘要', icon: '📋', prompt: 'summarize' },
  { id: 'simplify', label: '简化表达', icon: '🎯', prompt: 'simplify' },
  { id: 'formal', label: '正式化', icon: '👔', prompt: 'formal' },
  { id: 'casual', label: '口语化', icon: '💬', prompt: 'casual' },
  { id: 'fix', label: '修正语法', icon: '🔧', prompt: 'fix_grammar' },
  { id: 'custom', label: '自定义指令', icon: '⚡', prompt: null }
];

function renderActionMenu() {
  const container = document.getElementById('ai-actions');
  AI_ACTIONS.forEach(action => {
    const btn = document.createElement('button');
    btn.className = 'ai-action-btn';
    btn.innerHTML = `${action.icon} ${action.label}`;
    btn.onclick = () => executeAction(action);
    container.appendChild(btn);
  });
}
```

### 4.5 Electron 端 AI 处理

```typescript
// apps/electron/src/main/onlyoffice/ai-handler.ts
import { CraftAgent } from '@craft-agent/shared';

export async function handleOnlyOfficeAI(
  agent: CraftAgent,
  action: string,
  text: string,
  customPrompt?: string
): Promise<string> {
  const prompts: Record<string, string> = {
    translate: `将以下文本翻译为英文（如果是英文则翻译为中文），保持原有格式：\n\n${text}`,
    rewrite: `改写润色以下文本，使其更加流畅专业，保持原意：\n\n${text}`,
    expand: `扩展以下内容，添加更多细节和论述，保持风格一致：\n\n${text}`,
    summarize: `总结以下内容的要点，简明扼要：\n\n${text}`,
    simplify: `简化以下文本，使用更简单易懂的表达：\n\n${text}`,
    formal: `将以下文本改写为正式的书面语风格：\n\n${text}`,
    casual: `将以下文本改写为轻松的口语化风格：\n\n${text}`,
    fix_grammar: `修正以下文本中的语法和拼写错误：\n\n${text}`
  };

  const prompt = customPrompt || prompts[action] || text;
  
  const response = await agent.chat(prompt);
  return response.content;
}
```

### 4.6 webview 消息监听

```typescript
// apps/electron/src/renderer/components/onlyoffice/OnlyOfficeEditor.tsx
useEffect(() => {
  const handleMessage = async (event: MessageEvent) => {
    if (event.data.type !== 'onlyoffice-ai') return;
    
    const { channel, payload } = event.data;
    
    if (channel === 'ai:process') {
      setProcessing(true);
      try {
        const result = await window.electronAPI.onlyoffice.processAI({
          sessionId: currentSessionId,
          action: payload.action,
          text: payload.text
        });
        
        // 发回给 webview 中的插件
        webviewRef.current?.contentWindow?.postMessage({
          type: 'ai-response',
          payload: { action: 'replace', text: result }
        }, '*');
      } catch (error) {
        console.error('AI processing failed:', error);
      } finally {
        setProcessing(false);
      }
    }
  };
  
  window.addEventListener('message', handleMessage);
  return () => window.removeEventListener('message', handleMessage);
}, [currentSessionId]);
```

---

## Phase 5: 文件管理器集成（4 天）

### 5.1 安装 Chonky

```bash
bun add chonky chonky-icon-fontawesome
```

### 5.2 文件管理器 IPC

```typescript
// apps/electron/src/shared/types.ts - 新增 channels
export type IpcChannels = {
  // ... 现有 channels
  
  // 文件管理器
  'fm:listDirectory': (path: string) => Promise<FileEntry[]>;
  'fm:createFolder': (path: string, name: string) => Promise<void>;
  'fm:delete': (paths: string[]) => Promise<void>;
  'fm:rename': (oldPath: string, newPath: string) => Promise<void>;
  'fm:copy': (src: string[], dest: string) => Promise<void>;
  'fm:move': (src: string[], dest: string) => Promise<void>;
  'fm:openInEditor': (path: string) => Promise<void>;
};

export interface FileEntry {
  id: string;
  name: string;
  isDir: boolean;
  size?: number;
  modDate?: Date;
  ext?: string;
}
```

### 5.3 主进程实现

```typescript
// apps/electron/src/main/file-manager.ts
import { ipcMain } from 'electron';
import * as fs from 'fs/promises';
import * as path from 'path';

export function registerFileManagerIpc() {
  ipcMain.handle('fm:listDirectory', async (_, dirPath: string) => {
    const entries = await fs.readdir(dirPath, { withFileTypes: true });
    return entries
      .filter(e => !e.name.startsWith('.')) // 隐藏文件
      .map(e => ({
        id: path.join(dirPath, e.name),
        name: e.name,
        isDir: e.isDirectory(),
        ext: e.isDirectory() ? undefined : path.extname(e.name).toLowerCase()
      }));
  });
  
  ipcMain.handle('fm:delete', async (_, paths: string[]) => {
    for (const p of paths) {
      await fs.rm(p, { recursive: true });
    }
  });
  
  ipcMain.handle('fm:rename', async (_, oldPath: string, newPath: string) => {
    await fs.rename(oldPath, newPath);
  });
  
  // ... 其他操作
}
```

### 5.4 Chonky 文件浏览器组件

```typescript
// apps/electron/src/renderer/components/file-manager/FileManager.tsx
import { FileBrowser, FileNavbar, FileToolbar, FileList, FileContextMenu, ChonkyActions } from 'chonky';
import { ChonkyIconFA } from 'chonky-icon-fontawesome';

interface FileManagerProps {
  workspacePath: string;
  onFileOpen: (file: FileEntry) => void;
}

export function FileManager({ workspacePath, onFileOpen }: FileManagerProps) {
  const [currentPath, setCurrentPath] = useState(workspacePath);
  const [files, setFiles] = useState<FileArray>([]);
  
  useEffect(() => {
    loadDirectory(currentPath);
  }, [currentPath]);
  
  const loadDirectory = async (path: string) => {
    const entries = await window.electronAPI.fm.listDirectory(path);
    setFiles(entries.map(e => ({
      id: e.id,
      name: e.name,
      isDir: e.isDir,
      ext: e.ext
    })));
  };
  
  const handleFileAction = useCallback((data: ChonkyFileActionData) => {
    switch (data.id) {
      case ChonkyActions.OpenFiles.id:
        const file = data.payload.targetFile;
        if (file?.isDir) {
          setCurrentPath(file.id);
        } else if (file) {
          onFileOpen(file);
        }
        break;
      case ChonkyActions.DeleteFiles.id:
        handleDelete(data.state.selectedFiles);
        break;
      // ... 其他操作
    }
  }, []);
  
  // 自定义 AI 协作 Action
  const aiCollabAction: FileAction = {
    id: 'ai_collaborate',
    button: {
      name: 'AI 协作',
      toolbar: true,
      contextMenu: true,
      icon: ChonkyIconName.flash
    }
  };
  
  return (
    <FileBrowser
      files={files}
      folderChain={buildFolderChain(currentPath, workspacePath)}
      onFileAction={handleFileAction}
      fileActions={[aiCollabAction, ...ChonkyActions.DefaultFileActions]}
      iconComponent={ChonkyIconFA}
    >
      <FileNavbar />
      <FileToolbar />
      <FileList />
      <FileContextMenu />
    </FileBrowser>
  );
}
```

### 5.5 集成到侧边栏

```typescript
// apps/electron/src/renderer/components/app-shell/sidebar-types.ts
export type SidebarMode = 
  | { type: 'chat' }
  | { type: 'files' }  // 新增
  | { type: 'settings' };
```

```typescript
// apps/electron/src/renderer/components/app-shell/LeftSidebar.tsx
{sidebarMode.type === 'files' && (
  <FileManager
    workspacePath={currentWorkspace.path}
    onFileOpen={handleFileOpen}
  />
)}
```

### 5.6 文件打开路由

```typescript
const handleFileOpen = (file: FileEntry) => {
  const officeExtensions = ['.docx', '.xlsx', '.pptx', '.doc', '.xls', '.ppt', '.pdf'];
  
  if (officeExtensions.includes(file.ext || '')) {
    // 打开 OnlyOffice 编辑器
    navigate({
      type: 'document-editor',
      filePath: file.id
    });
  } else {
    // 打开代码/文本查看器
    navigate({
      type: 'file-viewer',
      filePath: file.id
    });
  }
};
```

---

## Phase 6: 测试与优化（3 天）

### 6.1 测试矩阵

**功能测试**
- [ ] OnlyOffice 连接与认证
- [ ] 文档打开（DOCX/XLSX/PPTX/PDF）
- [ ] 文档编辑与保存
- [ ] AI 功能（翻译/改写/扩展/总结）
- [ ] 文件管理器基本操作

**边界条件**
- [ ] 大文件处理（>50MB）
- [ ] 网络断开恢复
- [ ] 并发编辑冲突
- [ ] 特殊字符文件名

**性能测试**
- [ ] 文档加载时间 < 3s
- [ ] AI 响应时间 < 5s
- [ ] 内存占用监控

### 6.2 错误处理增强

```typescript
// apps/electron/src/renderer/components/onlyoffice/ErrorBoundary.tsx
export function OnlyOfficeErrorBoundary({ children }: PropsWithChildren) {
  const [error, setError] = useState<Error | null>(null);
  
  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-full">
        <AlertCircle className="w-12 h-12 text-destructive" />
        <h3 className="mt-4 text-lg font-medium">文档服务连接失败</h3>
        <p className="mt-2 text-muted-foreground">{error.message}</p>
        <Button className="mt-4" onClick={() => setError(null)}>
          重试
        </Button>
      </div>
    );
  }
  
  return <ErrorBoundary onError={setError}>{children}</ErrorBoundary>;
}
```

### 6.3 性能优化

**文件预加载**
```typescript
// 预先获取 JWT 和 documentKey，减少打开延迟
const prefetchDocument = async (filePath: string) => {
  const config = await window.electronAPI.onlyoffice.prepareDocument(filePath);
  documentCache.set(filePath, config);
};
```

**AI 响应流式输出**
```typescript
// 流式显示 AI 结果，提升感知速度
const streamAIResponse = async (text: string, action: string) => {
  const stream = await window.electronAPI.onlyoffice.processAIStream({
    text,
    action
  });
  
  let accumulated = '';
  for await (const chunk of stream) {
    accumulated += chunk;
    updatePreview(accumulated);
  }
  
  return accumulated;
};
```

---

## 附录

### A. 技术栈汇总

| 组件 | 技术选型 | 用途 |
|------|----------|------|
| 文档服务 | OnlyOffice Document Server | 云端文档编辑 |
| 文件存储 | 自建 DMS API | 文件读写回调 |
| 前端集成 | webview + postMessage | 编辑器嵌入 |
| AI 协作 | OnlyOffice Plugin | 插件式 AI 功能 |
| 文件管理 | Chonky | 文件浏览器 UI |
| 通信安全 | JWT | API 认证 |

### B. 工期汇总

| 阶段 | 工作日 |
|------|--------|
| Phase 1: 云端部署 | 2 天 |
| Phase 2: DMS 适配 | 3 天 |
| Phase 3: Electron 集成 | 3 天 |
| Phase 4: AI 插件开发 | 5 天 |
| Phase 5: 文件管理器 | 4 天 |
| Phase 6: 测试优化 | 3 天 |
| **总计** | **20 天** |

### C. 风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| OnlyOffice 服务不稳定 | 高 | 健康检查 + 自动重连 + 降级提示 |
| JWT 安全漏洞 | 高 | 短有效期 + HTTPS + 定期轮换密钥 |
| 大文件上传慢 | 中 | 分片上传 + 进度提示 |
| AI 响应慢 | 中 | 流式输出 + 超时取消 |
| 浏览器兼容性 | 低 | Electron webview 版本固定 |

### D. 后续扩展方向

1. **协同编辑**：多用户同时编辑同一文档
2. **版本历史**：文档版本管理与回滚
3. **批注与修订**：AI 辅助审阅建议
4. **模板系统**：预设文档模板快速创建
5. **批量处理**：多文件 AI 批量翻译/转换
