# 资源管理器与 Office 文档 AI 协作功能实现方案

## 概述

本方案采用 **Chonky 文件管理器 + HTML 中间层 + Tiptap 编辑器** 的轻量级架构，实现文件浏览、Office 文档预览编辑与 AI 协作功能。

### 方案对比

| 方案 | 优势 | 劣势 | 适用场景 |
|------|------|------|----------|
| **本方案（轻量）** | 零服务端依赖、快速上线、离线可用 | 复杂样式保真度有限 | MVP 快速验证 |
| OnlyOffice 云端 | 高保真、完整 Office 功能 | 需部署服务、网络依赖 | 企业级协作 |

### 核心架构

```
┌─────────────────────────────────────────────────────────────┐
│                    Sprouty AI Electron                      │
│                                                              │
│  ┌─────────────────┐    ┌─────────────────────────────────┐ │
│  │   左侧边栏       │    │          主内容区                │ │
│  │  ┌───────────┐  │    │  ┌─────────────────────────────┐│ │
│  │  │ Chonky    │  │    │  │  Tiptap 富文本编辑器        ││ │
│  │  │ 文件浏览器 │  │───▶│  │  + AI 协作工具栏            ││ │
│  │  │           │  │    │  │  + 实时预览                  ││ │
│  │  └───────────┘  │    │  └─────────────────────────────┘│ │
│  └─────────────────┘    └─────────────────────────────────┘ │
│                                                              │
│  ┌──────────────────────────────────────────────────────────┐│
│  │                      Main Process                         ││
│  │  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐ ││
│  │  │ 文件操作  │  │ 格式转换  │  │ AI 处理  │  │ 导出服务  │ ││
│  │  │ IPC      │  │ Service  │  │ Agent    │  │ Service  │ ││
│  │  └──────────┘  └──────────┘  └──────────┘  └──────────┘ ││
│  └──────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────┘
```

### 技术栈

| 组件 | 技术选型 | 说明 |
|------|----------|------|
| 文件浏览器 | Chonky + chonky-icon-fontawesome | 成熟的 React 文件管理器 |
| 富文本编辑 | Tiptap (ProseMirror) | 可扩展、支持自定义 AI 扩展 |
| DOCX 读取 | mammoth.js | DOCX → HTML |
| DOCX 写入 | @anthropic/html-to-docx 或 docx | HTML → DOCX |
| Excel 处理 | SheetJS (xlsx) | 读写 XLSX/XLS |
| PDF 预览 | pdfjs-dist | PDF 渲染 |
| PDF 导出 | Puppeteer (主进程) | HTML → PDF |
| 代码高亮 | Shiki | 代码文件预览 |

---

## Phase 1: 文件管理器基础（5 天）

### 1.1 安装依赖

```bash
bun add chonky chonky-icon-fontawesome
```

### 1.2 IPC 通道定义

```typescript
// apps/electron/src/shared/types.ts
export type FileManagerChannels = {
  // 目录操作
  'fm:listDirectory': (path: string) => Promise<FileEntry[]>;
  'fm:createFolder': (parentPath: string, name: string) => Promise<string>;
  'fm:delete': (paths: string[]) => Promise<void>;
  'fm:rename': (oldPath: string, newName: string) => Promise<string>;
  'fm:move': (srcPaths: string[], destDir: string) => Promise<void>;
  'fm:copy': (srcPaths: string[], destDir: string) => Promise<void>;
  
  // 文件操作
  'fm:readFile': (path: string) => Promise<FileContent>;
  'fm:writeFile': (path: string, content: string | Buffer) => Promise<void>;
  'fm:getFileInfo': (path: string) => Promise<FileInfo>;
  
  // 监听
  'fm:watchDirectory': (path: string) => void;
  'fm:unwatchDirectory': (path: string) => void;
  'fm:onDirectoryChanged': (callback: (event: DirectoryChangeEvent) => void) => void;
};

export interface FileEntry {
  id: string;           // 完整路径作为 ID
  name: string;
  isDir: boolean;
  size?: number;
  modDate?: Date;
  ext?: string;
  thumbnailUrl?: string;
}

export interface FileInfo {
  path: string;
  name: string;
  size: number;
  isDir: boolean;
  createdAt: Date;
  modifiedAt: Date;
  mimeType?: string;
}

export interface FileContent {
  content: string | Buffer;
  encoding: 'utf-8' | 'base64' | 'binary';
  mimeType: string;
}

export interface DirectoryChangeEvent {
  type: 'add' | 'change' | 'unlink' | 'addDir' | 'unlinkDir';
  path: string;
}
```

### 1.3 主进程文件操作实现

```typescript
// apps/electron/src/main/file-manager/index.ts
import { ipcMain } from 'electron';
import * as fs from 'fs/promises';
import * as path from 'path';
import { watch, FSWatcher } from 'chokidar';
import mime from 'mime-types';

const watchers = new Map<string, FSWatcher>();

export function registerFileManagerIpc(mainWindow: BrowserWindow) {
  // 列出目录内容
  ipcMain.handle('fm:listDirectory', async (_, dirPath: string) => {
    const entries = await fs.readdir(dirPath, { withFileTypes: true });
    const result: FileEntry[] = [];
    
    for (const entry of entries) {
      // 跳过隐藏文件
      if (entry.name.startsWith('.')) continue;
      
      const fullPath = path.join(dirPath, entry.name);
      const stats = await fs.stat(fullPath);
      
      result.push({
        id: fullPath,
        name: entry.name,
        isDir: entry.isDirectory(),
        size: entry.isDirectory() ? undefined : stats.size,
        modDate: stats.mtime,
        ext: entry.isDirectory() ? undefined : path.extname(entry.name).toLowerCase()
      });
    }
    
    // 排序：文件夹在前，按名称字母序
    return result.sort((a, b) => {
      if (a.isDir !== b.isDir) return a.isDir ? -1 : 1;
      return a.name.localeCompare(b.name);
    });
  });
  
  // 创建文件夹
  ipcMain.handle('fm:createFolder', async (_, parentPath: string, name: string) => {
    const newPath = path.join(parentPath, name);
    await fs.mkdir(newPath, { recursive: true });
    return newPath;
  });
  
  // 删除文件/文件夹
  ipcMain.handle('fm:delete', async (_, paths: string[]) => {
    for (const p of paths) {
      await fs.rm(p, { recursive: true, force: true });
    }
  });
  
  // 重命名
  ipcMain.handle('fm:rename', async (_, oldPath: string, newName: string) => {
    const dir = path.dirname(oldPath);
    const newPath = path.join(dir, newName);
    await fs.rename(oldPath, newPath);
    return newPath;
  });
  
  // 移动
  ipcMain.handle('fm:move', async (_, srcPaths: string[], destDir: string) => {
    for (const src of srcPaths) {
      const dest = path.join(destDir, path.basename(src));
      await fs.rename(src, dest);
    }
  });
  
  // 复制
  ipcMain.handle('fm:copy', async (_, srcPaths: string[], destDir: string) => {
    for (const src of srcPaths) {
      const dest = path.join(destDir, path.basename(src));
      await fs.cp(src, dest, { recursive: true });
    }
  });
  
  // 读取文件
  ipcMain.handle('fm:readFile', async (_, filePath: string) => {
    const mimeType = mime.lookup(filePath) || 'application/octet-stream';
    const isText = mimeType.startsWith('text/') || 
                   ['application/json', 'application/javascript', 'application/xml'].includes(mimeType);
    
    if (isText) {
      const content = await fs.readFile(filePath, 'utf-8');
      return { content, encoding: 'utf-8', mimeType };
    } else {
      const content = await fs.readFile(filePath);
      return { content: content.toString('base64'), encoding: 'base64', mimeType };
    }
  });
  
  // 写入文件
  ipcMain.handle('fm:writeFile', async (_, filePath: string, content: string | Buffer) => {
    await fs.writeFile(filePath, content);
  });
  
  // 获取文件信息
  ipcMain.handle('fm:getFileInfo', async (_, filePath: string) => {
    const stats = await fs.stat(filePath);
    return {
      path: filePath,
      name: path.basename(filePath),
      size: stats.size,
      isDir: stats.isDirectory(),
      createdAt: stats.birthtime,
      modifiedAt: stats.mtime,
      mimeType: mime.lookup(filePath) || undefined
    };
  });
  
  // 监听目录变化
  ipcMain.on('fm:watchDirectory', (_, dirPath: string) => {
    if (watchers.has(dirPath)) return;
    
    const watcher = watch(dirPath, {
      depth: 1,
      ignoreInitial: true,
      ignored: /(^|[\/\\])\../  // 忽略隐藏文件
    });
    
    watcher.on('all', (eventType, changedPath) => {
      mainWindow.webContents.send('fm:directoryChanged', {
        type: eventType,
        path: changedPath
      });
    });
    
    watchers.set(dirPath, watcher);
  });
  
  ipcMain.on('fm:unwatchDirectory', (_, dirPath: string) => {
    const watcher = watchers.get(dirPath);
    if (watcher) {
      watcher.close();
      watchers.delete(dirPath);
    }
  });
}
```

### 1.4 Chonky 文件浏览器组件

```typescript
// apps/electron/src/renderer/components/file-manager/FileManager.tsx
import { useState, useEffect, useCallback, useMemo } from 'react';
import {
  FileBrowser,
  FileNavbar,
  FileToolbar,
  FileList,
  FileContextMenu,
  ChonkyActions,
  defineFileAction,
  ChonkyIconName,
  FileData,
  FileAction,
  ChonkyFileActionData
} from 'chonky';
import { ChonkyIconFA } from 'chonky-icon-fontawesome';
import { useT } from '../../hooks/useT';

interface FileManagerProps {
  rootPath: string;
  onFileOpen: (file: FileData) => void;
  onAICollaborate?: (files: FileData[]) => void;
}

export function FileManager({ rootPath, onFileOpen, onAICollaborate }: FileManagerProps) {
  const t = useT();
  const [currentPath, setCurrentPath] = useState(rootPath);
  const [files, setFiles] = useState<FileData[]>([]);
  const [loading, setLoading] = useState(false);

  // 加载目录
  const loadDirectory = useCallback(async (dirPath: string) => {
    setLoading(true);
    try {
      const entries = await window.electronAPI.fm.listDirectory(dirPath);
      setFiles(entries.map(e => ({
        id: e.id,
        name: e.name,
        isDir: e.isDir,
        size: e.size,
        modDate: e.modDate,
        ext: e.ext
      })));
    } catch (error) {
      console.error('Failed to load directory:', error);
    } finally {
      setLoading(false);
    }
  }, []);

  // 初始加载 & 路径变化时重新加载
  useEffect(() => {
    loadDirectory(currentPath);
  }, [currentPath, loadDirectory]);

  // 目录变化监听
  useEffect(() => {
    window.electronAPI.fm.watchDirectory(currentPath);
    
    const unsubscribe = window.electronAPI.fm.onDirectoryChanged(() => {
      loadDirectory(currentPath);
    });
    
    return () => {
      window.electronAPI.fm.unwatchDirectory(currentPath);
      unsubscribe();
    };
  }, [currentPath, loadDirectory]);

  // 构建面包屑导航
  const folderChain = useMemo(() => {
    const chain: FileData[] = [];
    let p = currentPath;
    
    while (p !== rootPath && p !== '/') {
      chain.unshift({
        id: p,
        name: p.split('/').pop() || p,
        isDir: true
      });
      p = p.substring(0, p.lastIndexOf('/')) || '/';
    }
    
    chain.unshift({
      id: rootPath,
      name: t('工作区'),
      isDir: true
    });
    
    return chain;
  }, [currentPath, rootPath, t]);

  // 自定义 AI 协作动作
  const aiCollaborateAction = defineFileAction({
    id: 'ai_collaborate',
    button: {
      name: t('AI 协作'),
      toolbar: true,
      contextMenu: true,
      icon: ChonkyIconName.flash
    }
  } as const);

  // 处理文件操作
  const handleFileAction = useCallback(async (data: ChonkyFileActionData) => {
    const { id, payload, state } = data;
    
    switch (id) {
      case ChonkyActions.OpenFiles.id: {
        const target = payload.targetFile;
        if (!target) return;
        
        if (target.isDir) {
          setCurrentPath(target.id);
        } else {
          onFileOpen(target);
        }
        break;
      }
      
      case ChonkyActions.DeleteFiles.id: {
        const paths = state.selectedFiles.map(f => f.id);
        if (paths.length === 0) return;
        
        const confirmed = await window.electronAPI.dialog.confirm({
          title: t('确认删除'),
          message: t('确定要删除选中的 {{count}} 个项目吗？', { count: paths.length })
        });
        
        if (confirmed) {
          await window.electronAPI.fm.delete(paths);
          loadDirectory(currentPath);
        }
        break;
      }
      
      case ChonkyActions.CreateFolder.id: {
        const name = await window.electronAPI.dialog.prompt({
          title: t('新建文件夹'),
          message: t('请输入文件夹名称')
        });
        
        if (name) {
          await window.electronAPI.fm.createFolder(currentPath, name);
          loadDirectory(currentPath);
        }
        break;
      }
      
      case ChonkyActions.MoveFiles.id: {
        const { files: filesToMove, destination } = payload;
        if (!destination?.isDir) return;
        
        await window.electronAPI.fm.move(
          filesToMove.map((f: FileData) => f.id),
          destination.id
        );
        loadDirectory(currentPath);
        break;
      }
      
      case 'ai_collaborate': {
        const selectedFiles = state.selectedFiles.filter(f => !f.isDir);
        if (selectedFiles.length > 0 && onAICollaborate) {
          onAICollaborate(selectedFiles);
        }
        break;
      }
    }
  }, [currentPath, loadDirectory, onFileOpen, onAICollaborate, t]);

  // 文件操作列表
  const fileActions: FileAction[] = useMemo(() => [
    aiCollaborateAction,
    ChonkyActions.CreateFolder,
    ChonkyActions.DeleteFiles,
    ChonkyActions.DownloadFiles,
    ChonkyActions.CopyFiles,
    ChonkyActions.MoveFiles
  ], [aiCollaborateAction]);

  return (
    <div className="h-full flex flex-col">
      <FileBrowser
        files={files}
        folderChain={folderChain}
        fileActions={fileActions}
        onFileAction={handleFileAction}
        iconComponent={ChonkyIconFA}
        disableDragAndDrop={false}
        disableDefaultFileActions={false}
      >
        <FileNavbar />
        <FileToolbar />
        <FileList />
        <FileContextMenu />
      </FileBrowser>
      
      {loading && (
        <div className="absolute inset-0 bg-background/50 flex items-center justify-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-primary" />
        </div>
      )}
    </div>
  );
}
```

### 1.5 侧边栏集成

```typescript
// apps/electron/src/renderer/components/app-shell/sidebar-types.ts
export type SidebarMode = 
  | { type: 'chat' }
  | { type: 'files' }      // 新增：文件管理器
  | { type: 'sources' }
  | { type: 'settings' };
```

```typescript
// apps/electron/src/renderer/components/app-shell/LeftSidebar.tsx
import { FileManager } from '../file-manager/FileManager';

// 在侧边栏渲染中添加
{sidebarMode.type === 'files' && (
  <FileManager
    rootPath={currentWorkspace?.path || '/'}
    onFileOpen={handleFileOpen}
    onAICollaborate={handleAICollaborate}
  />
)}
```

### 1.6 导航菜单项

```typescript
// 在侧边栏导航中添加文件管理器入口
<SidebarNavItem
  icon={<FolderOpen className="w-4 h-4" />}
  label={t('文件')}
  active={sidebarMode.type === 'files'}
  onClick={() => setSidebarMode({ type: 'files' })}
/>
```

---

## Phase 2: 文档预览功能（4 天）

### 2.1 安装依赖

```bash
bun add mammoth xlsx pdfjs-dist shiki
```

### 2.2 文档预览路由

```typescript
// apps/electron/src/renderer/components/document-viewer/DocumentViewer.tsx
import { useState, useEffect } from 'react';
import { DocxViewer } from './DocxViewer';
import { XlsxViewer } from './XlsxViewer';
import { PdfViewer } from './PdfViewer';
import { CodeViewer } from './CodeViewer';
import { ImageViewer } from './ImageViewer';
import { TextViewer } from './TextViewer';

interface DocumentViewerProps {
  filePath: string;
  onEdit?: () => void;
}

export function DocumentViewer({ filePath, onEdit }: DocumentViewerProps) {
  const ext = filePath.split('.').pop()?.toLowerCase();
  
  // 根据扩展名选择查看器
  switch (ext) {
    case 'docx':
    case 'doc':
      return <DocxViewer filePath={filePath} onEdit={onEdit} />;
    
    case 'xlsx':
    case 'xls':
      return <XlsxViewer filePath={filePath} onEdit={onEdit} />;
    
    case 'pdf':
      return <PdfViewer filePath={filePath} />;
    
    case 'png':
    case 'jpg':
    case 'jpeg':
    case 'gif':
    case 'webp':
    case 'svg':
      return <ImageViewer filePath={filePath} />;
    
    case 'js':
    case 'ts':
    case 'tsx':
    case 'jsx':
    case 'py':
    case 'go':
    case 'rs':
    case 'java':
    case 'c':
    case 'cpp':
    case 'h':
    case 'css':
    case 'scss':
    case 'html':
    case 'json':
    case 'yaml':
    case 'yml':
    case 'toml':
    case 'xml':
    case 'sql':
    case 'sh':
    case 'bash':
    case 'zsh':
      return <CodeViewer filePath={filePath} />;
    
    default:
      return <TextViewer filePath={filePath} />;
  }
}
```

### 2.3 DOCX 预览器

```typescript
// apps/electron/src/renderer/components/document-viewer/DocxViewer.tsx
import { useState, useEffect } from 'react';
import mammoth from 'mammoth';
import { Button } from '../ui/button';
import { Loader2, Edit } from 'lucide-react';
import { useT } from '../../hooks/useT';

interface DocxViewerProps {
  filePath: string;
  onEdit?: () => void;
}

export function DocxViewer({ filePath, onEdit }: DocxViewerProps) {
  const t = useT();
  const [html, setHtml] = useState<string>('');
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadDocx = async () => {
      setLoading(true);
      setError(null);
      
      try {
        // 从主进程读取文件
        const { content } = await window.electronAPI.fm.readFile(filePath);
        const arrayBuffer = Uint8Array.from(atob(content), c => c.charCodeAt(0)).buffer;
        
        // 转换为 HTML
        const result = await mammoth.convertToHtml(
          { arrayBuffer },
          {
            styleMap: [
              "p[style-name='Heading 1'] => h1:fresh",
              "p[style-name='Heading 2'] => h2:fresh",
              "p[style-name='Heading 3'] => h3:fresh",
              "b => strong",
              "i => em",
              "u => u"
            ]
          }
        );
        
        setHtml(result.value);
        
        // 记录警告
        if (result.messages.length > 0) {
          console.warn('Mammoth warnings:', result.messages);
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : t('文档加载失败'));
      } finally {
        setLoading(false);
      }
    };
    
    loadDocx();
  }, [filePath, t]);

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-destructive">
        {error}
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col">
      {/* 工具栏 */}
      <div className="flex items-center justify-between p-2 border-b">
        <span className="text-sm text-muted-foreground">
          {filePath.split('/').pop()}
        </span>
        {onEdit && (
          <Button size="sm" variant="outline" onClick={onEdit}>
            <Edit className="w-4 h-4 mr-2" />
            {t('编辑')}
          </Button>
        )}
      </div>
      
      {/* 文档内容 */}
      <div className="flex-1 overflow-auto p-8">
        <div 
          className="prose prose-sm dark:prose-invert max-w-none"
          dangerouslySetInnerHTML={{ __html: html }}
        />
      </div>
    </div>
  );
}
```

### 2.4 Excel 预览器

```typescript
// apps/electron/src/renderer/components/document-viewer/XlsxViewer.tsx
import { useState, useEffect } from 'react';
import * as XLSX from 'xlsx';
import { Button } from '../ui/button';
import { Loader2, Edit } from 'lucide-react';
import { useT } from '../../hooks/useT';

interface XlsxViewerProps {
  filePath: string;
  onEdit?: () => void;
}

export function XlsxViewer({ filePath, onEdit }: XlsxViewerProps) {
  const t = useT();
  const [workbook, setWorkbook] = useState<XLSX.WorkBook | null>(null);
  const [activeSheet, setActiveSheet] = useState<string>('');
  const [tableHtml, setTableHtml] = useState<string>('');
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadXlsx = async () => {
      setLoading(true);
      setError(null);
      
      try {
        const { content } = await window.electronAPI.fm.readFile(filePath);
        const arrayBuffer = Uint8Array.from(atob(content), c => c.charCodeAt(0)).buffer;
        
        const wb = XLSX.read(arrayBuffer, { type: 'array' });
        setWorkbook(wb);
        setActiveSheet(wb.SheetNames[0]);
      } catch (err) {
        setError(err instanceof Error ? err.message : t('文档加载失败'));
      } finally {
        setLoading(false);
      }
    };
    
    loadXlsx();
  }, [filePath, t]);

  // 切换工作表时更新 HTML
  useEffect(() => {
    if (!workbook || !activeSheet) return;
    
    const sheet = workbook.Sheets[activeSheet];
    const html = XLSX.utils.sheet_to_html(sheet, { editable: false });
    setTableHtml(html);
  }, [workbook, activeSheet]);

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-destructive">
        {error}
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col">
      {/* 工具栏 */}
      <div className="flex items-center justify-between p-2 border-b">
        <div className="flex items-center gap-2">
          {/* 工作表标签 */}
          {workbook?.SheetNames.map(name => (
            <Button
              key={name}
              size="sm"
              variant={name === activeSheet ? 'default' : 'ghost'}
              onClick={() => setActiveSheet(name)}
            >
              {name}
            </Button>
          ))}
        </div>
        {onEdit && (
          <Button size="sm" variant="outline" onClick={onEdit}>
            <Edit className="w-4 h-4 mr-2" />
            {t('编辑')}
          </Button>
        )}
      </div>
      
      {/* 表格内容 */}
      <div className="flex-1 overflow-auto">
        <div 
          className="xlsx-table"
          dangerouslySetInnerHTML={{ __html: tableHtml }}
        />
      </div>
    </div>
  );
}
```

### 2.5 PDF 预览器

```typescript
// apps/electron/src/renderer/components/document-viewer/PdfViewer.tsx
import { useState, useEffect, useRef } from 'react';
import * as pdfjsLib from 'pdfjs-dist';
import { Loader2, ZoomIn, ZoomOut, ChevronLeft, ChevronRight } from 'lucide-react';
import { Button } from '../ui/button';
import { useT } from '../../hooks/useT';

// 设置 worker
pdfjsLib.GlobalWorkerOptions.workerSrc = new URL(
  'pdfjs-dist/build/pdf.worker.mjs',
  import.meta.url
).toString();

interface PdfViewerProps {
  filePath: string;
}

export function PdfViewer({ filePath }: PdfViewerProps) {
  const t = useT();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [pdf, setPdf] = useState<pdfjsLib.PDFDocumentProxy | null>(null);
  const [currentPage, setCurrentPage] = useState(1);
  const [totalPages, setTotalPages] = useState(0);
  const [scale, setScale] = useState(1.5);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // 加载 PDF
  useEffect(() => {
    const loadPdf = async () => {
      setLoading(true);
      setError(null);
      
      try {
        const { content } = await window.electronAPI.fm.readFile(filePath);
        const data = Uint8Array.from(atob(content), c => c.charCodeAt(0));
        
        const pdfDoc = await pdfjsLib.getDocument({ data }).promise;
        setPdf(pdfDoc);
        setTotalPages(pdfDoc.numPages);
        setCurrentPage(1);
      } catch (err) {
        setError(err instanceof Error ? err.message : t('PDF 加载失败'));
      } finally {
        setLoading(false);
      }
    };
    
    loadPdf();
  }, [filePath, t]);

  // 渲染当前页
  useEffect(() => {
    if (!pdf || !canvasRef.current) return;
    
    const renderPage = async () => {
      const page = await pdf.getPage(currentPage);
      const viewport = page.getViewport({ scale });
      
      const canvas = canvasRef.current!;
      const context = canvas.getContext('2d')!;
      
      canvas.height = viewport.height;
      canvas.width = viewport.width;
      
      await page.render({
        canvasContext: context,
        viewport
      }).promise;
    };
    
    renderPage();
  }, [pdf, currentPage, scale]);

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-destructive">
        {error}
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col">
      {/* 工具栏 */}
      <div className="flex items-center justify-center gap-4 p-2 border-b">
        <Button
          size="icon"
          variant="ghost"
          onClick={() => setCurrentPage(p => Math.max(1, p - 1))}
          disabled={currentPage <= 1}
        >
          <ChevronLeft className="w-4 h-4" />
        </Button>
        
        <span className="text-sm">
          {currentPage} / {totalPages}
        </span>
        
        <Button
          size="icon"
          variant="ghost"
          onClick={() => setCurrentPage(p => Math.min(totalPages, p + 1))}
          disabled={currentPage >= totalPages}
        >
          <ChevronRight className="w-4 h-4" />
        </Button>
        
        <div className="w-px h-4 bg-border" />
        
        <Button
          size="icon"
          variant="ghost"
          onClick={() => setScale(s => Math.max(0.5, s - 0.25))}
        >
          <ZoomOut className="w-4 h-4" />
        </Button>
        
        <span className="text-sm w-12 text-center">
          {Math.round(scale * 100)}%
        </span>
        
        <Button
          size="icon"
          variant="ghost"
          onClick={() => setScale(s => Math.min(3, s + 0.25))}
        >
          <ZoomIn className="w-4 h-4" />
        </Button>
      </div>
      
      {/* PDF 画布 */}
      <div className="flex-1 overflow-auto flex justify-center p-4 bg-muted/30">
        <canvas ref={canvasRef} className="shadow-lg" />
      </div>
    </div>
  );
}
```

### 2.6 代码预览器

```typescript
// apps/electron/src/renderer/components/document-viewer/CodeViewer.tsx
import { useState, useEffect } from 'react';
import { codeToHtml } from 'shiki';
import { Loader2 } from 'lucide-react';
import { useT } from '../../hooks/useT';

interface CodeViewerProps {
  filePath: string;
}

const EXT_TO_LANG: Record<string, string> = {
  js: 'javascript',
  ts: 'typescript',
  tsx: 'tsx',
  jsx: 'jsx',
  py: 'python',
  go: 'go',
  rs: 'rust',
  java: 'java',
  c: 'c',
  cpp: 'cpp',
  h: 'c',
  css: 'css',
  scss: 'scss',
  html: 'html',
  json: 'json',
  yaml: 'yaml',
  yml: 'yaml',
  toml: 'toml',
  xml: 'xml',
  sql: 'sql',
  sh: 'bash',
  bash: 'bash',
  zsh: 'bash'
};

export function CodeViewer({ filePath }: CodeViewerProps) {
  const t = useT();
  const [html, setHtml] = useState<string>('');
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadCode = async () => {
      setLoading(true);
      setError(null);
      
      try {
        const { content } = await window.electronAPI.fm.readFile(filePath);
        const ext = filePath.split('.').pop()?.toLowerCase() || '';
        const lang = EXT_TO_LANG[ext] || 'text';
        
        const highlighted = await codeToHtml(content, {
          lang,
          theme: 'github-dark'
        });
        
        setHtml(highlighted);
      } catch (err) {
        setError(err instanceof Error ? err.message : t('文件加载失败'));
      } finally {
        setLoading(false);
      }
    };
    
    loadCode();
  }, [filePath, t]);

  if (loading) {
    return (
      <div className="flex items-center justify-center h-full">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-destructive">
        {error}
      </div>
    );
  }

  return (
    <div className="h-full overflow-auto">
      <div 
        className="shiki-code p-4 text-sm"
        dangerouslySetInnerHTML={{ __html: html }}
      />
    </div>
  );
}
```

---

## Phase 3: Tiptap 文档编辑器（5 天）

### 3.1 安装依赖

```bash
bun add @tiptap/react @tiptap/starter-kit @tiptap/extension-image @tiptap/extension-table @tiptap/extension-table-row @tiptap/extension-table-cell @tiptap/extension-table-header @tiptap/extension-text-align @tiptap/extension-underline @tiptap/extension-highlight @tiptap/extension-placeholder
```

### 3.2 Tiptap 编辑器组件

```typescript
// apps/electron/src/renderer/components/document-editor/TiptapEditor.tsx
import { useEditor, EditorContent, Editor } from '@tiptap/react';
import StarterKit from '@tiptap/starter-kit';
import Image from '@tiptap/extension-image';
import Table from '@tiptap/extension-table';
import TableRow from '@tiptap/extension-table-row';
import TableCell from '@tiptap/extension-table-cell';
import TableHeader from '@tiptap/extension-table-header';
import TextAlign from '@tiptap/extension-text-align';
import Underline from '@tiptap/extension-underline';
import Highlight from '@tiptap/extension-highlight';
import Placeholder from '@tiptap/extension-placeholder';
import { EditorToolbar } from './EditorToolbar';
import { AIToolbar } from './AIToolbar';
import { useT } from '../../hooks/useT';

interface TiptapEditorProps {
  initialContent: string;
  onChange?: (html: string) => void;
  onSave?: (html: string) => void;
  sessionId?: string;  // AI 协作需要
}

export function TiptapEditor({ 
  initialContent, 
  onChange, 
  onSave,
  sessionId 
}: TiptapEditorProps) {
  const t = useT();
  
  const editor = useEditor({
    extensions: [
      StarterKit.configure({
        heading: {
          levels: [1, 2, 3, 4]
        }
      }),
      Image,
      Table.configure({
        resizable: true
      }),
      TableRow,
      TableCell,
      TableHeader,
      TextAlign.configure({
        types: ['heading', 'paragraph']
      }),
      Underline,
      Highlight.configure({
        multicolor: true
      }),
      Placeholder.configure({
        placeholder: t('开始输入...')
      })
    ],
    content: initialContent,
    onUpdate: ({ editor }) => {
      onChange?.(editor.getHTML());
    }
  });

  if (!editor) {
    return null;
  }

  return (
    <div className="h-full flex flex-col">
      {/* 格式工具栏 */}
      <EditorToolbar editor={editor} onSave={() => onSave?.(editor.getHTML())} />
      
      {/* AI 工具栏 */}
      {sessionId && (
        <AIToolbar editor={editor} sessionId={sessionId} />
      )}
      
      {/* 编辑区域 */}
      <div className="flex-1 overflow-auto p-8">
        <EditorContent 
          editor={editor} 
          className="prose prose-sm dark:prose-invert max-w-none min-h-full focus:outline-none"
        />
      </div>
    </div>
  );
}
```

### 3.3 格式工具栏

```typescript
// apps/electron/src/renderer/components/document-editor/EditorToolbar.tsx
import { Editor } from '@tiptap/react';
import { 
  Bold, Italic, Underline, Strikethrough,
  Heading1, Heading2, Heading3,
  List, ListOrdered,
  AlignLeft, AlignCenter, AlignRight, AlignJustify,
  Undo, Redo, Save,
  Table, Image
} from 'lucide-react';
import { Button } from '../ui/button';
import { Separator } from '../ui/separator';
import { useT } from '../../hooks/useT';

interface EditorToolbarProps {
  editor: Editor;
  onSave?: () => void;
}

export function EditorToolbar({ editor, onSave }: EditorToolbarProps) {
  const t = useT();
  
  return (
    <div className="flex items-center gap-1 p-2 border-b flex-wrap">
      {/* 撤销/重做 */}
      <Button
        size="icon"
        variant="ghost"
        onClick={() => editor.chain().focus().undo().run()}
        disabled={!editor.can().undo()}
      >
        <Undo className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant="ghost"
        onClick={() => editor.chain().focus().redo().run()}
        disabled={!editor.can().redo()}
      >
        <Redo className="w-4 h-4" />
      </Button>
      
      <Separator orientation="vertical" className="h-6 mx-2" />
      
      {/* 标题 */}
      <Button
        size="icon"
        variant={editor.isActive('heading', { level: 1 }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleHeading({ level: 1 }).run()}
      >
        <Heading1 className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('heading', { level: 2 }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleHeading({ level: 2 }).run()}
      >
        <Heading2 className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('heading', { level: 3 }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleHeading({ level: 3 }).run()}
      >
        <Heading3 className="w-4 h-4" />
      </Button>
      
      <Separator orientation="vertical" className="h-6 mx-2" />
      
      {/* 文本格式 */}
      <Button
        size="icon"
        variant={editor.isActive('bold') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleBold().run()}
      >
        <Bold className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('italic') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleItalic().run()}
      >
        <Italic className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('underline') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleUnderline().run()}
      >
        <Underline className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('strike') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleStrike().run()}
      >
        <Strikethrough className="w-4 h-4" />
      </Button>
      
      <Separator orientation="vertical" className="h-6 mx-2" />
      
      {/* 对齐 */}
      <Button
        size="icon"
        variant={editor.isActive({ textAlign: 'left' }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().setTextAlign('left').run()}
      >
        <AlignLeft className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive({ textAlign: 'center' }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().setTextAlign('center').run()}
      >
        <AlignCenter className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive({ textAlign: 'right' }) ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().setTextAlign('right').run()}
      >
        <AlignRight className="w-4 h-4" />
      </Button>
      
      <Separator orientation="vertical" className="h-6 mx-2" />
      
      {/* 列表 */}
      <Button
        size="icon"
        variant={editor.isActive('bulletList') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleBulletList().run()}
      >
        <List className="w-4 h-4" />
      </Button>
      <Button
        size="icon"
        variant={editor.isActive('orderedList') ? 'secondary' : 'ghost'}
        onClick={() => editor.chain().focus().toggleOrderedList().run()}
      >
        <ListOrdered className="w-4 h-4" />
      </Button>
      
      <div className="flex-1" />
      
      {/* 保存 */}
      {onSave && (
        <Button size="sm" onClick={onSave}>
          <Save className="w-4 h-4 mr-2" />
          {t('保存')}
        </Button>
      )}
    </div>
  );
}
```

---

## Phase 4: AI 协作功能（4 天）

### 4.1 AI 工具栏组件

```typescript
// apps/electron/src/renderer/components/document-editor/AIToolbar.tsx
import { useState } from 'react';
import { Editor } from '@tiptap/react';
import { 
  Languages, Sparkles, Expand, FileText, 
  Wand2, Briefcase, MessageCircle, CheckCheck,
  Loader2, ChevronDown
} from 'lucide-react';
import { Button } from '../ui/button';
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '../ui/dropdown-menu';
import { useT } from '../../hooks/useT';

interface AIToolbarProps {
  editor: Editor;
  sessionId: string;
}

const AI_ACTIONS = [
  { id: 'translate', label: '翻译', icon: Languages },
  { id: 'rewrite', label: '改写润色', icon: Sparkles },
  { id: 'expand', label: '扩展内容', icon: Expand },
  { id: 'summarize', label: '总结摘要', icon: FileText },
  { id: 'simplify', label: '简化表达', icon: Wand2 },
  { id: 'formal', label: '正式化', icon: Briefcase },
  { id: 'casual', label: '口语化', icon: MessageCircle },
  { id: 'fix', label: '修正语法', icon: CheckCheck },
];

export function AIToolbar({ editor, sessionId }: AIToolbarProps) {
  const t = useT();
  const [processing, setProcessing] = useState(false);
  const [currentAction, setCurrentAction] = useState<string | null>(null);

  const getSelectedText = (): string => {
    const { from, to } = editor.state.selection;
    return editor.state.doc.textBetween(from, to, ' ');
  };

  const handleAIAction = async (actionId: string) => {
    const selectedText = getSelectedText();
    
    if (!selectedText.trim()) {
      // 没有选中文本时，使用全文
      const fullText = editor.getText();
      if (!fullText.trim()) return;
    }
    
    setProcessing(true);
    setCurrentAction(actionId);
    
    try {
      const text = selectedText.trim() || editor.getText();
      
      const result = await window.electronAPI.ai.processDocument({
        sessionId,
        action: actionId,
        text
      });
      
      if (selectedText.trim()) {
        // 替换选中文本
        editor.chain().focus().deleteSelection().insertContent(result).run();
      } else {
        // 替换全文
        editor.commands.setContent(result);
      }
    } catch (error) {
      console.error('AI processing failed:', error);
    } finally {
      setProcessing(false);
      setCurrentAction(null);
    }
  };

  const hasSelection = !editor.state.selection.empty;

  return (
    <div className="flex items-center gap-2 p-2 border-b bg-muted/30">
      <span className="text-sm text-muted-foreground mr-2">
        {t('AI 协作')}:
      </span>
      
      {/* 快捷按钮 */}
      <Button
        size="sm"
        variant="outline"
        onClick={() => handleAIAction('translate')}
        disabled={processing}
      >
        {processing && currentAction === 'translate' ? (
          <Loader2 className="w-4 h-4 mr-2 animate-spin" />
        ) : (
          <Languages className="w-4 h-4 mr-2" />
        )}
        {t('翻译')}
      </Button>
      
      <Button
        size="sm"
        variant="outline"
        onClick={() => handleAIAction('rewrite')}
        disabled={processing}
      >
        {processing && currentAction === 'rewrite' ? (
          <Loader2 className="w-4 h-4 mr-2 animate-spin" />
        ) : (
          <Sparkles className="w-4 h-4 mr-2" />
        )}
        {t('润色')}
      </Button>
      
      {/* 更多操作 */}
      <DropdownMenu>
        <DropdownMenuTrigger asChild>
          <Button size="sm" variant="outline" disabled={processing}>
            {t('更多')}
            <ChevronDown className="w-4 h-4 ml-2" />
          </Button>
        </DropdownMenuTrigger>
        <DropdownMenuContent>
          {AI_ACTIONS.slice(2).map(action => (
            <DropdownMenuItem
              key={action.id}
              onClick={() => handleAIAction(action.id)}
            >
              <action.icon className="w-4 h-4 mr-2" />
              {t(action.label)}
            </DropdownMenuItem>
          ))}
        </DropdownMenuContent>
      </DropdownMenu>
      
      <div className="flex-1" />
      
      <span className="text-xs text-muted-foreground">
        {hasSelection ? t('将处理选中文本') : t('将处理全文')}
      </span>
    </div>
  );
}
```

### 4.2 AI 处理 IPC

```typescript
// apps/electron/src/shared/types.ts
export type AIDocumentChannels = {
  'ai:processDocument': (params: AIDocumentParams) => Promise<string>;
  'ai:processDocumentStream': (params: AIDocumentParams) => AsyncIterable<string>;
};

export interface AIDocumentParams {
  sessionId: string;
  action: string;
  text: string;
  customPrompt?: string;
}
```

### 4.3 主进程 AI 处理

```typescript
// apps/electron/src/main/ai/document-handler.ts
import { ipcMain } from 'electron';
import { getSessionAgent } from '../sessions';

const PROMPTS: Record<string, (text: string) => string> = {
  translate: (text) => 
    `将以下文本翻译为英文（如果是英文则翻译为中文），保持原有格式，只返回翻译结果：\n\n${text}`,
  
  rewrite: (text) => 
    `改写润色以下文本，使其更加流畅专业，保持原意，只返回改写结果：\n\n${text}`,
  
  expand: (text) => 
    `扩展以下内容，添加更多细节和论述，保持风格一致，只返回扩展后的内容：\n\n${text}`,
  
  summarize: (text) => 
    `总结以下内容的要点，简明扼要，只返回总结内容：\n\n${text}`,
  
  simplify: (text) => 
    `简化以下文本，使用更简单易懂的表达，只返回简化结果：\n\n${text}`,
  
  formal: (text) => 
    `将以下文本改写为正式的书面语风格，只返回改写结果：\n\n${text}`,
  
  casual: (text) => 
    `将以下文本改写为轻松的口语化风格，只返回改写结果：\n\n${text}`,
  
  fix: (text) => 
    `修正以下文本中的语法和拼写错误，只返回修正后的文本：\n\n${text}`
};

export function registerAIDocumentIpc() {
  ipcMain.handle('ai:processDocument', async (_, params: AIDocumentParams) => {
    const { sessionId, action, text, customPrompt } = params;
    
    const agent = await getSessionAgent(sessionId);
    if (!agent) {
      throw new Error('Session not found');
    }
    
    const promptFn = PROMPTS[action];
    const prompt = customPrompt || (promptFn ? promptFn(text) : text);
    
    const response = await agent.chat(prompt);
    return response.content;
  });
}
```

---

## Phase 5: 文档导出功能（4 天）

### 5.1 安装依赖

```bash
bun add html-to-docx puppeteer-core
```

### 5.2 导出服务

```typescript
// apps/electron/src/main/export/document-export.ts
import { ipcMain } from 'electron';
import * as fs from 'fs/promises';
import * as path from 'path';
import HTMLtoDOCX from 'html-to-docx';
import * as XLSX from 'xlsx';

export function registerExportIpc() {
  // 导出为 DOCX
  ipcMain.handle('export:toDocx', async (_, html: string, savePath: string) => {
    const docxBuffer = await HTMLtoDOCX(html, null, {
      table: { row: { cantSplit: true } },
      footer: true,
      pageNumber: true
    });
    
    await fs.writeFile(savePath, Buffer.from(docxBuffer));
    return savePath;
  });
  
  // 导出为 PDF
  ipcMain.handle('export:toPdf', async (_, html: string, savePath: string) => {
    // 使用 Puppeteer 生成 PDF
    const puppeteer = await import('puppeteer-core');
    
    // 查找 Chromium
    const executablePath = await findChromium();
    
    const browser = await puppeteer.launch({
      executablePath,
      headless: true
    });
    
    try {
      const page = await browser.newPage();
      await page.setContent(wrapHtmlForPdf(html), {
        waitUntil: 'networkidle0'
      });
      
      await page.pdf({
        path: savePath,
        format: 'A4',
        margin: {
          top: '2cm',
          right: '2cm',
          bottom: '2cm',
          left: '2cm'
        },
        printBackground: true
      });
      
      return savePath;
    } finally {
      await browser.close();
    }
  });
  
  // 导出为 XLSX
  ipcMain.handle('export:toXlsx', async (_, data: any[][], savePath: string) => {
    const worksheet = XLSX.utils.aoa_to_sheet(data);
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(workbook, worksheet, 'Sheet1');
    
    const buffer = XLSX.write(workbook, { type: 'buffer', bookType: 'xlsx' });
    await fs.writeFile(savePath, buffer);
    
    return savePath;
  });
}

function wrapHtmlForPdf(html: string): string {
  return `
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="UTF-8">
      <style>
        body {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
          font-size: 12pt;
          line-height: 1.6;
          color: #333;
        }
        h1 { font-size: 24pt; margin-top: 24pt; }
        h2 { font-size: 18pt; margin-top: 18pt; }
        h3 { font-size: 14pt; margin-top: 14pt; }
        table { border-collapse: collapse; width: 100%; }
        td, th { border: 1px solid #ddd; padding: 8px; }
        img { max-width: 100%; }
      </style>
    </head>
    <body>
      ${html}
    </body>
    </html>
  `;
}

async function findChromium(): Promise<string> {
  // macOS
  const macPaths = [
    '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
    '/Applications/Chromium.app/Contents/MacOS/Chromium'
  ];
  
  for (const p of macPaths) {
    try {
      await fs.access(p);
      return p;
    } catch {}
  }
  
  throw new Error('Chromium not found');
}
```

### 5.3 导出对话框

```typescript
// apps/electron/src/renderer/components/document-editor/ExportDialog.tsx
import { useState } from 'react';
import { FileDown, FileText, FileSpreadsheet } from 'lucide-react';
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from '../ui/dialog';
import { Button } from '../ui/button';
import { useT } from '../../hooks/useT';

interface ExportDialogProps {
  open: boolean;
  onClose: () => void;
  html: string;
  originalPath?: string;
}

export function ExportDialog({ open, onClose, html, originalPath }: ExportDialogProps) {
  const t = useT();
  const [exporting, setExporting] = useState(false);

  const handleExport = async (format: 'docx' | 'pdf') => {
    setExporting(true);
    
    try {
      const defaultName = originalPath 
        ? originalPath.replace(/\.[^.]+$/, `.${format}`)
        : `document.${format}`;
      
      const savePath = await window.electronAPI.dialog.saveFile({
        defaultPath: defaultName,
        filters: [
          format === 'docx' 
            ? { name: 'Word Document', extensions: ['docx'] }
            : { name: 'PDF Document', extensions: ['pdf'] }
        ]
      });
      
      if (!savePath) return;
      
      if (format === 'docx') {
        await window.electronAPI.export.toDocx(html, savePath);
      } else {
        await window.electronAPI.export.toPdf(html, savePath);
      }
      
      // 打开文件所在位置
      await window.electronAPI.shell.showItemInFolder(savePath);
      onClose();
    } catch (error) {
      console.error('Export failed:', error);
    } finally {
      setExporting(false);
    }
  };

  return (
    <Dialog open={open} onOpenChange={onClose}>
      <DialogContent>
        <DialogHeader>
          <DialogTitle>{t('导出文档')}</DialogTitle>
        </DialogHeader>
        
        <div className="grid grid-cols-2 gap-4 py-4">
          <Button
            variant="outline"
            className="h-24 flex-col gap-2"
            onClick={() => handleExport('docx')}
            disabled={exporting}
          >
            <FileText className="w-8 h-8" />
            <span>{t('导出为 DOCX')}</span>
          </Button>
          
          <Button
            variant="outline"
            className="h-24 flex-col gap-2"
            onClick={() => handleExport('pdf')}
            disabled={exporting}
          >
            <FileDown className="w-8 h-8" />
            <span>{t('导出为 PDF')}</span>
          </Button>
        </div>
      </DialogContent>
    </Dialog>
  );
}
```

---

## Phase 6: 测试与优化（3 天）

### 6.1 测试清单

**文件管理器**
- [ ] 目录浏览与导航
- [ ] 文件创建/删除/重命名
- [ ] 拖拽移动/复制
- [ ] 右键菜单操作
- [ ] 目录变化实时刷新

**文档预览**
- [ ] DOCX 预览
- [ ] XLSX 预览（多工作表）
- [ ] PDF 预览（翻页/缩放）
- [ ] 代码文件语法高亮
- [ ] 图片预览

**文档编辑**
- [ ] 富文本编辑（标题/粗体/列表等）
- [ ] 撤销/重做
- [ ] 复制/粘贴格式保留

**AI 协作**
- [ ] 翻译功能
- [ ] 改写润色
- [ ] 扩展/总结/简化
- [ ] 选区处理 vs 全文处理

**导出功能**
- [ ] 导出 DOCX
- [ ] 导出 PDF
- [ ] 格式保真度

### 6.2 性能优化

**大文件处理**
```typescript
// 分块加载大文件
const CHUNK_SIZE = 1024 * 1024; // 1MB

async function loadLargeFile(path: string) {
  const info = await window.electronAPI.fm.getFileInfo(path);
  
  if (info.size > CHUNK_SIZE * 10) {
    // 超过 10MB 提示用户
    const confirmed = await window.electronAPI.dialog.confirm({
      message: '文件较大，加载可能需要一些时间，是否继续？'
    });
    if (!confirmed) return null;
  }
  
  return window.electronAPI.fm.readFile(path);
}
```

**编辑器防抖保存**
```typescript
// 自动保存（防抖）
const debouncedSave = useMemo(
  () => debounce((html: string) => {
    onSave?.(html);
  }, 2000),
  [onSave]
);

// 在 onChange 中调用
const handleChange = (html: string) => {
  onChange?.(html);
  debouncedSave(html);
};
```

---

## 附录

### A. 工期汇总

| 阶段 | 内容 | 工作日 |
|------|------|--------|
| Phase 1 | 文件管理器基础 (Chonky) | 5 天 |
| Phase 2 | 文档预览功能 | 4 天 |
| Phase 3 | Tiptap 文档编辑器 | 5 天 |
| Phase 4 | AI 协作功能 | 4 天 |
| Phase 5 | 文档导出功能 | 4 天 |
| Phase 6 | 测试与优化 | 3 天 |
| **总计** | | **25 天** |

### B. 格式支持矩阵

| 格式 | 预览 | 编辑 | AI 协作 | 导出 |
|------|------|------|---------|------|
| DOCX | ✅ | ✅ | ✅ | ✅ |
| XLSX | ✅ | ⚠️ (简单编辑) | ✅ (文本) | ✅ |
| PDF | ✅ | ❌ | ✅ (OCR) | ✅ |
| TXT/MD | ✅ | ✅ | ✅ | ✅ |
| 代码文件 | ✅ | ✅ | ✅ | ❌ |
| 图片 | ✅ | ❌ | ❌ | ❌ |

### C. 与 OnlyOffice 方案对比

| 维度 | 本方案 (轻量) | OnlyOffice 方案 |
|------|---------------|-----------------|
| 部署复杂度 | 低（纯前端） | 中（需云端服务） |
| 格式保真度 | 中（复杂样式可能丢失） | 高 |
| 离线支持 | ✅ 完全支持 | ❌ 需网络 |
| 协同编辑 | ❌ | ✅ |
| 修订/批注 | ❌ | ✅ |
| 开发周期 | 25 天 | 20 天 |
| 运维成本 | 低 | 中 |

### D. 后续迭代方向

1. **Markdown 支持**：增加 Markdown 编辑与预览
2. **版本历史**：本地文件版本管理
3. **模板系统**：常用文档模板
4. **批量处理**：多文件 AI 批量操作
5. **云端同步**：可选云存储集成
