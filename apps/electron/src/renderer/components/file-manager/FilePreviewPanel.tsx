/**
 * FilePreviewPanel - 文件预览面板
 *
 * 功能:
 * - 文本/代码文件预览与编辑
 * - Markdown 渲染预览 + 转 Word/PDF
 * - 图片预览（缩放/旋转）
 * - PDF 内嵌渲染（翻页）
 * - Office 文档预览（docx/xlsx）
 * - 视频/音频播放
 * - 文件元数据显示
 */

import * as React from 'react'
import { useState, useEffect, useCallback, useRef, useMemo } from 'react'
import {
  File,
  FileText,
  FileCode,
  Image as ImageIcon,
  ExternalLink,
  RefreshCw,
  FileQuestion,
  Info,
  FileVideo,
  FileAudio,
  FileType,
  FileSpreadsheet,
  ZoomIn,
  ZoomOut,
  RotateCw,
  Maximize2,
  Pencil,
  Save,
  X,
  ChevronLeft,
  ChevronRight,
  Eye,
  Code,
  FileDown,
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { useT } from '@/context/LocaleContext'
import { Markdown } from '@/components/markdown'
import { Document, Page, pdfjs } from 'react-pdf'
import 'react-pdf/dist/Page/AnnotationLayer.css'
import 'react-pdf/dist/Page/TextLayer.css'
import { toast } from 'sonner'
import type { FMFileEntry } from '../../../shared/types'

// 配置 pdfjs worker
pdfjs.GlobalWorkerOptions.workerSrc = new URL('pdfjs-dist/build/pdf.worker.min.mjs', import.meta.url).toString()

export interface FilePreviewPanelProps {
  /** Selected file to preview */
  file?: FMFileEntry | null
  className?: string
  /** Markdown 转换请求回调 */
  onConvertRequest?: (filePath: string, targetFormat: 'word' | 'pdf') => void
  /** 关闭预览回调（嵌入模式使用） */
  onClose?: () => void
}

/**
 * Format file size in human-readable format
 */
function formatFileSize(bytes?: number): string {
  if (bytes === undefined || bytes < 0) return '-'
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(1)} GB`
}

/**
 * Format date in locale-friendly format
 */
function formatDate(timestamp?: number): string {
  if (!timestamp) return '-'
  return new Date(timestamp).toLocaleString()
}

/**
 * Get file type category
 */
function getFileType(file: FMFileEntry): 'text' | 'markdown' | 'image' | 'code' | 'video' | 'audio' | 'pdf' | 'office' | 'unknown' {
  const ext = file.name.split('.').pop()?.toLowerCase()

  // Markdown types (从 text 中拆出)
  if (['md', 'markdown', 'mdx'].includes(ext || '')) {
    return 'markdown'
  }

  // Text types
  if (['txt', 'log', 'csv', 'ini', 'conf', 'cfg'].includes(ext || '')) {
    return 'text'
  }

  // Image types
  if (['png', 'jpg', 'jpeg', 'gif', 'webp', 'svg', 'ico', 'bmp', 'tiff', 'tif'].includes(ext || '')) {
    return 'image'
  }

  // Code types
  if (['ts', 'tsx', 'js', 'jsx', 'json', 'yaml', 'yml', 'py', 'rb', 'go', 'rs', 'java', 'c', 'cpp', 'h', 'css', 'scss', 'less', 'html', 'xml', 'sh', 'bash', 'zsh', 'vue', 'svelte', 'swift', 'kt', 'php', 'sql', 'graphql', 'toml', 'env', 'dockerfile', 'makefile'].includes(ext || '')) {
    return 'code'
  }

  // Video types
  if (['mp4', 'webm', 'mov', 'avi', 'mkv', 'm4v', 'ogv'].includes(ext || '')) {
    return 'video'
  }

  // Audio types
  if (['mp3', 'wav', 'ogg', 'flac', 'aac', 'm4a', 'wma'].includes(ext || '')) {
    return 'audio'
  }

  // PDF
  if (ext === 'pdf') {
    return 'pdf'
  }

  // Office types
  if (['doc', 'docx', 'xls', 'xlsx', 'ppt', 'pptx'].includes(ext || '')) {
    return 'office'
  }

  return 'unknown'
}

/**
 * 判断文件类型是否可编辑
 */
function isEditable(fileType: string): boolean {
  return ['text', 'markdown', 'code'].includes(fileType)
}

/**
 * Get language for syntax highlighting
 */
function getLanguage(file: FMFileEntry): string {
  const ext = file.name.split('.').pop()?.toLowerCase()
  
  const langMap: Record<string, string> = {
    ts: 'typescript',
    tsx: 'tsx',
    js: 'javascript',
    jsx: 'jsx',
    json: 'json',
    yaml: 'yaml',
    yml: 'yaml',
    py: 'python',
    rb: 'ruby',
    go: 'go',
    rs: 'rust',
    java: 'java',
    c: 'c',
    cpp: 'cpp',
    h: 'c',
    css: 'css',
    scss: 'scss',
    html: 'html',
    xml: 'xml',
    sh: 'bash',
    bash: 'bash',
    zsh: 'bash',
    md: 'markdown',
    markdown: 'markdown',
    vue: 'vue',
    svelte: 'svelte',
  }
  
  return langMap[ext || ''] || 'text'
}

/**
 * 文本/代码预览与编辑组件
 */
function TextPreviewEdit({ file, isEditing, onEditingChange }: { file: FMFileEntry; isEditing: boolean; onEditingChange: (editing: boolean) => void }) {
  const t = useT()
  const [content, setContent] = useState<string | null>(null)
  const [editContent, setEditContent] = useState<string>('')
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [isTruncated, setIsTruncated] = useState(false)
  const [isSaving, setIsSaving] = useState(false)
  const [hasUnsavedChanges, setHasUnsavedChanges] = useState(false)
  const textareaRef = useRef<HTMLTextAreaElement>(null)

  const MAX_PREVIEW_SIZE = 100 * 1024 // 100KB limit for preview

  useEffect(() => {
    const loadContent = async () => {
      setIsLoading(true)
      setError(null)
      setIsTruncated(false)

      try {
        if (file.size && file.size > MAX_PREVIEW_SIZE) {
          setIsTruncated(true)
        }

        const text = await window.electronAPI.readFile(file.path)

        if (text.length > MAX_PREVIEW_SIZE) {
          setContent(text.slice(0, MAX_PREVIEW_SIZE))
          setIsTruncated(true)
        } else {
          setContent(text)
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }

    loadContent()
  }, [file.path, file.size])

  // 进入编辑模式时加载完整内容
  const prevIsEditingRef = useRef(false)
  useEffect(() => {
    // 仅在 isEditing 从 false → true 时初始化编辑内容
    if (isEditing && !prevIsEditingRef.current && content !== null) {
      setEditContent(content)
      setHasUnsavedChanges(false)
      setError(null)
      setTimeout(() => textareaRef.current?.focus(), 50)
    }
    prevIsEditingRef.current = isEditing
  }, [isEditing, content])

  const handleSave = useCallback(async () => {
    setIsSaving(true)
    try {
      await window.electronAPI.fm.writeFile(file.path, editContent)
      setContent(editContent)
      setHasUnsavedChanges(false)
      onEditingChange(false)
      toast.success(t('保存成功'))
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err)
      setError(msg)
      toast.error(t('保存失败'))
    } finally {
      setIsSaving(false)
    }
  }, [file.path, editContent, onEditingChange, t])

  const handleCancel = useCallback(() => {
    setHasUnsavedChanges(false)
    onEditingChange(false)
  }, [onEditingChange])

  // Ctrl/Cmd+S 保存, Escape 取消
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if ((e.metaKey || e.ctrlKey) && e.key === 's') {
      e.preventDefault()
      handleSave()
    } else if (e.key === 'Escape') {
      e.preventDefault()
      handleCancel()
    }
  }, [handleSave, handleCancel])

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  if (error && !isEditing) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileQuestion className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此文件')}</p>
          <p className="text-xs text-destructive mt-1">{error}</p>
        </div>
      </div>
    )
  }

  const fileType = getFileType(file)
  const language = getLanguage(file)

  // 编辑模式
  if (isEditing) {
    return (
      <div className="h-full flex flex-col titlebar-no-drag" onKeyDown={handleKeyDown}>
        {/* 编辑工具栏 */}
        <div className="flex items-center justify-between px-3 py-1.5 bg-muted/50 border-b shrink-0 relative z-panel">
          <div className="flex items-center gap-2">
            <span className="text-xs font-medium text-muted-foreground">{t('编辑模式')}</span>
            {hasUnsavedChanges && (
              <span className="text-xs text-orange-500">{t('未保存')}</span>
            )}
            {error && (
              <span className="text-xs text-destructive">{error}</span>
            )}
          </div>
          <div className="flex items-center gap-1">
            <Button variant="ghost" size="sm" className="h-6 text-xs" onClick={handleCancel}>
              <X className="h-3 w-3 mr-1" />
              {t('取消')}
            </Button>
            <Button variant="default" size="sm" className="h-6 text-xs" onClick={handleSave} disabled={isSaving}>
              <Save className="h-3 w-3 mr-1" />
              {isSaving ? t('保存中...') : t('保存')}
            </Button>
          </div>
        </div>
        <textarea
          ref={textareaRef}
          value={editContent}
          onChange={(e) => {
            setEditContent(e.target.value)
            setHasUnsavedChanges(true)
          }}
          className="flex-1 w-full p-4 text-sm font-mono bg-background resize-none outline-none border-none"
          spellCheck={false}
        />
      </div>
    )
  }

  // 预览模式
  return (
    <div className="h-full flex flex-col">
      {isTruncated && (
        <div className="px-3 py-1.5 bg-muted/50 text-xs text-muted-foreground border-b flex items-center gap-1">
          <Info className="h-3 w-3" />
          {t('文件过大，仅显示部分内容')}
        </div>
      )}
      <ScrollArea className="flex-1">
        <pre className={cn(
          "p-4 text-sm font-mono whitespace-pre-wrap break-all",
          fileType === 'code' && "bg-muted/30"
        )}>
          <code className={`language-${language}`}>
            {content}
          </code>
        </pre>
      </ScrollArea>
    </div>
  )
}

/**
 * Markdown 预览组件（渲染预览 + 源码 + 编辑 + 转换）
 */
function MarkdownPreview({ file, isEditing, onEditingChange, onConvertRequest }: {
  file: FMFileEntry
  isEditing: boolean
  onEditingChange: (editing: boolean) => void
  onConvertRequest?: (filePath: string, targetFormat: 'word' | 'pdf') => void
}) {
  const t = useT()
  const [content, setContent] = useState<string | null>(null)
  const [editContent, setEditContent] = useState<string>('')
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [viewMode, setViewMode] = useState<'preview' | 'source'>('preview')
  const [isSaving, setIsSaving] = useState(false)
  const [hasUnsavedChanges, setHasUnsavedChanges] = useState(false)
  const textareaRef = useRef<HTMLTextAreaElement>(null)

  useEffect(() => {
    const loadContent = async () => {
      setIsLoading(true)
      setError(null)
      try {
        const text = await window.electronAPI.readFile(file.path)
        setContent(text)
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }
    loadContent()
  }, [file.path])

  const prevIsEditingRef = useRef(false)
  useEffect(() => {
    if (isEditing && !prevIsEditingRef.current && content !== null) {
      setEditContent(content)
      setHasUnsavedChanges(false)
      setError(null)
      setTimeout(() => textareaRef.current?.focus(), 50)
    }
    prevIsEditingRef.current = isEditing
  }, [isEditing, content])

  const handleSave = useCallback(async () => {
    setIsSaving(true)
    try {
      await window.electronAPI.fm.writeFile(file.path, editContent)
      setContent(editContent)
      setHasUnsavedChanges(false)
      onEditingChange(false)
      toast.success(t('保存成功'))
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err)
      setError(msg)
      toast.error(t('保存失败'))
    } finally {
      setIsSaving(false)
    }
  }, [file.path, editContent, onEditingChange, t])

  const handleCancel = useCallback(() => {
    setHasUnsavedChanges(false)
    onEditingChange(false)
  }, [onEditingChange])

  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if ((e.metaKey || e.ctrlKey) && e.key === 's') {
      e.preventDefault()
      handleSave()
    } else if (e.key === 'Escape') {
      e.preventDefault()
      handleCancel()
    }
  }, [handleSave, handleCancel])

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  if (error && !isEditing) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileQuestion className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此文件')}</p>
          <p className="text-xs text-destructive mt-1">{error}</p>
        </div>
      </div>
    )
  }

  // 编辑模式
  if (isEditing) {
    return (
      <div className="h-full flex flex-col titlebar-no-drag" onKeyDown={handleKeyDown}>
        <div className="flex items-center justify-between px-3 py-1.5 bg-muted/50 border-b shrink-0 relative z-panel">
          <div className="flex items-center gap-2">
            <span className="text-xs font-medium text-muted-foreground">{t('编辑模式')}</span>
            {hasUnsavedChanges && (
              <span className="text-xs text-orange-500">{t('未保存')}</span>
            )}
            {error && (
              <span className="text-xs text-destructive">{error}</span>
            )}
          </div>
          <div className="flex items-center gap-1">
            <Button variant="ghost" size="sm" className="h-6 text-xs" onClick={handleCancel}>
              <X className="h-3 w-3 mr-1" />
              {t('取消')}
            </Button>
            <Button variant="default" size="sm" className="h-6 text-xs" onClick={handleSave} disabled={isSaving}>
              <Save className="h-3 w-3 mr-1" />
              {isSaving ? t('保存中...') : t('保存')}
            </Button>
          </div>
        </div>
        <textarea
          ref={textareaRef}
          value={editContent}
          onChange={(e) => {
            setEditContent(e.target.value)
            setHasUnsavedChanges(true)
          }}
          className="flex-1 w-full p-4 text-sm font-mono bg-background resize-none outline-none border-none"
          spellCheck={false}
        />
      </div>
    )
  }

  // 预览模式
  return (
    <div className="h-full flex flex-col">
      {/* 工具栏 */}
      <div className="flex items-center justify-between px-2 py-1.5 border-b bg-muted/30 shrink-0 relative z-panel">
        <div className="flex items-center gap-1 titlebar-no-drag">
          <Button
            variant={viewMode === 'preview' ? 'secondary' : 'ghost'}
            size="sm"
            className="h-6 text-xs"
            onClick={() => setViewMode('preview')}
          >
            <Eye className="h-3 w-3 mr-1" />
            {t('预览')}
          </Button>
          <Button
            variant={viewMode === 'source' ? 'secondary' : 'ghost'}
            size="sm"
            className="h-6 text-xs"
            onClick={() => setViewMode('source')}
          >
            <Code className="h-3 w-3 mr-1" />
            {t('源码')}
          </Button>
        </div>
        {onConvertRequest && (
          <div className="flex items-center gap-1 titlebar-no-drag">
            <Button
              variant="ghost"
              size="sm"
              className="h-6 text-xs"
              onClick={() => onConvertRequest(file.path, 'word')}
            >
              <FileDown className="h-3 w-3 mr-1" />
              {t('转 Word')}
            </Button>
            <Button
              variant="ghost"
              size="sm"
              className="h-6 text-xs"
              onClick={() => onConvertRequest(file.path, 'pdf')}
            >
              <FileDown className="h-3 w-3 mr-1" />
              {t('转 PDF')}
            </Button>
          </div>
        )}
      </div>

      {/* 内容区 */}
      {viewMode === 'preview' ? (
        <ScrollArea className="flex-1">
          <div className="p-4 prose prose-sm dark:prose-invert max-w-none">
            <Markdown>{content || ''}</Markdown>
          </div>
        </ScrollArea>
      ) : (
        <ScrollArea className="flex-1">
          <pre className="p-4 text-sm font-mono whitespace-pre-wrap break-all">
            <code className="language-markdown">{content}</code>
          </pre>
        </ScrollArea>
      )}
    </div>
  )
}

/**
 * Image preview component with zoom and rotate controls
 */
function ImagePreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [imageSrc, setImageSrc] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [scale, setScale] = useState(1)
  const [rotation, setRotation] = useState(0)
  
  // Reset zoom and rotation when file changes
  useEffect(() => {
    setScale(1)
    setRotation(0)
  }, [file.path])
  
  useEffect(() => {
    const loadImage = async () => {
      setIsLoading(true)
      setError(null)
      
      try {
        // Load image via IPC as base64 (file:// protocol blocked in Electron renderer)
        const { data, mimeType } = await window.electronAPI.fm.readFileBase64(file.path)
        const src = `data:${mimeType};base64,${data}`
        setImageSrc(src)
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }
    
    loadImage()
  }, [file.path])
  
  const handleZoomIn = useCallback(() => {
    setScale(prev => Math.min(prev * 1.25, 5))
  }, [])
  
  const handleZoomOut = useCallback(() => {
    setScale(prev => Math.max(prev / 1.25, 0.1))
  }, [])
  
  const handleRotate = useCallback(() => {
    setRotation(prev => (prev + 90) % 360)
  }, [])
  
  const handleReset = useCallback(() => {
    setScale(1)
    setRotation(0)
  }, [])
  
  // Handle mouse wheel for zoom
  const handleWheel = useCallback((e: React.WheelEvent) => {
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault()
      if (e.deltaY < 0) {
        setScale(prev => Math.min(prev * 1.1, 5))
      } else {
        setScale(prev => Math.max(prev / 1.1, 0.1))
      }
    }
  }, [])
  
  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }
  
  if (error || !imageSrc) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <ImageIcon className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此图片')}</p>
          {error && <p className="text-xs text-destructive mt-1">{error}</p>}
        </div>
      </div>
    )
  }
  
  return (
    <div className="h-full flex flex-col">
      {/* Image toolbar */}
      <div className="flex items-center justify-center gap-1 px-2 py-1.5 border-b bg-muted/30 shrink-0 relative z-panel titlebar-no-drag">
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={handleZoomOut}
          title={t('缩小')}
        >
          <ZoomOut className="h-4 w-4" />
        </Button>
        <span className="text-xs text-muted-foreground w-12 text-center">
          {Math.round(scale * 100)}%
        </span>
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={handleZoomIn}
          title={t('放大')}
        >
          <ZoomIn className="h-4 w-4" />
        </Button>
        <div className="w-px h-4 bg-border mx-1" />
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={handleRotate}
          title={t('旋转 90°')}
        >
          <RotateCw className="h-4 w-4" />
        </Button>
        <div className="w-px h-4 bg-border mx-1" />
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={handleReset}
          disabled={scale === 1 && rotation === 0}
          title={t('重置')}
        >
          <Maximize2 className="h-4 w-4" />
        </Button>
      </div>
      
      {/* Image container */}
      <div 
        className="flex-1 overflow-auto flex items-center justify-center p-4 bg-[repeating-conic-gradient(#80808010_0%_25%,transparent_0%_50%)] bg-[length:16px_16px]"
        onWheel={handleWheel}
      >
        <img
          src={imageSrc}
          alt={file.name}
          className="max-w-none rounded shadow-sm transition-transform duration-200"
          style={{
            transform: `scale(${scale}) rotate(${rotation}deg)`,
            maxWidth: scale <= 1 ? '100%' : 'none',
            maxHeight: scale <= 1 ? '100%' : 'none',
          }}
          onError={() => setError('Failed to load image')}
          draggable={false}
        />
      </div>
    </div>
  )
}

/**
 * Video preview component
 */
function VideoPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [error, setError] = useState<string | null>(null)
  const [videoSrc, setVideoSrc] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)

  useEffect(() => {
    setError(null)
    setVideoSrc(null)
    setIsLoading(true)
    let blobUrl: string | null = null
    let cancelled = false
    const loadVideo = async () => {
      try {
        const url = `localfile://file/${encodeURIComponent(file.path)}`
        const resp = await fetch(url)
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`)
        const blob = await resp.blob()
        if (cancelled) return
        blobUrl = URL.createObjectURL(blob)
        setVideoSrc(blobUrl)
      } catch (err) {
        if (cancelled) return
        console.error('Video load error:', err)
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        if (!cancelled) setIsLoading(false)
      }
    }
    loadVideo()
    return () => {
      cancelled = true
      if (blobUrl) URL.revokeObjectURL(blobUrl)
    }
  }, [file.path])

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileVideo className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此视频')}</p>
          <p className="text-xs text-destructive mt-1">{error}</p>
        </div>
      </div>
    )
  }

  if (isLoading || !videoSrc) {
    return (
      <div className="flex items-center justify-center h-full text-muted-foreground">
        {t('加载中...')}
      </div>
    )
  }

  return (
    <div className="h-full flex items-center justify-center p-4 bg-black/5">
      <video
        key={file.path}
        src={videoSrc}
        controls
        className="max-w-full max-h-full rounded shadow-sm"
        onError={() => setError(t('无法加载此视频文件'))}
      />
    </div>
  )
}

/**
 * Audio preview component
 */
function AudioPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [error, setError] = useState<string | null>(null)
  const [audioSrc, setAudioSrc] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)

  useEffect(() => {
    setError(null)
    setAudioSrc(null)
    setIsLoading(true)
    let blobUrl: string | null = null
    let cancelled = false
    const loadAudio = async () => {
      try {
        const url = `localfile://file/${encodeURIComponent(file.path)}`
        const resp = await fetch(url)
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`)
        const blob = await resp.blob()
        if (cancelled) return
        blobUrl = URL.createObjectURL(blob)
        setAudioSrc(blobUrl)
      } catch (err) {
        if (cancelled) return
        console.error('Audio load error:', err)
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        if (!cancelled) setIsLoading(false)
      }
    }
    loadAudio()
    return () => {
      cancelled = true
      if (blobUrl) URL.revokeObjectURL(blobUrl)
    }
  }, [file.path])

  if (error) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileAudio className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此音频')}</p>
          <p className="text-xs text-destructive mt-1">{error}</p>
        </div>
      </div>
    )
  }

  if (isLoading || !audioSrc) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <FileAudio className="h-16 w-16 text-muted-foreground mb-4" />
        <p className="text-sm text-muted-foreground">{t('加载中...')}</p>
      </div>
    )
  }

  return (
    <div className="h-full flex flex-col items-center justify-center p-4">
      <FileAudio className="h-16 w-16 text-muted-foreground mb-4" />
      <audio
        key={file.path}
        src={audioSrc}
        controls
        className="w-full max-w-md"
        onError={() => setError(t('无法加载此音频文件'))}
      />
    </div>
  )
}

/**
 * PDF 内嵌预览组件（react-pdf 渲染 + 翻页）
 */
function PdfPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [pdfData, setPdfData] = useState<Uint8Array | null>(null)
  const [numPages, setNumPages] = useState(0)
  const [currentPage, setCurrentPage] = useState(1)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const containerRef = useRef<HTMLDivElement>(null)
  const [containerWidth, setContainerWidth] = useState(0)

  useEffect(() => {
    const loadPdf = async () => {
      setIsLoading(true)
      setError(null)
      setCurrentPage(1)
      try {
        const data = await window.electronAPI.readFileBinary(file.path)
        setPdfData(data)
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }
    loadPdf()
  }, [file.path])

  // 监听容器宽度以自适应 PDF 页面宽度
  useEffect(() => {
    if (!containerRef.current) return
    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setContainerWidth(entry.contentRect.width)
      }
    })
    observer.observe(containerRef.current)
    return () => observer.disconnect()
  }, [])

  // 复制一份 ArrayBuffer，避免 pdf.js worker postMessage 后原 buffer 被 detached
  const pdfFile = useMemo(() => {
    if (!pdfData) return null
    return { data: pdfData.slice(0) }
  }, [pdfData])

  const handleOpen = useCallback(() => {
    window.electronAPI.openFile(file.path)
  }, [file.path])

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  // 加载失败时降级为"使用默认应用打开"
  if (error || !pdfData) {
    return (
      <div className="h-full flex flex-col items-center justify-center p-4 text-center">
        <FileType className="h-16 w-16 text-muted-foreground mb-4" />
        <p className="text-sm font-medium mb-1">{file.name}</p>
        <p className="text-xs text-muted-foreground mb-4">{t('PDF 文档')}</p>
        {error && <p className="text-xs text-destructive mb-2">{error}</p>}
        <Button variant="outline" size="sm" onClick={handleOpen}>
          <ExternalLink className="h-3.5 w-3.5 mr-1.5" />
          {t('使用默认应用打开')}
        </Button>
      </div>
    )
  }

  return (
    <div className="h-full flex flex-col">
      {/* 翻页工具栏 */}
      <div className="flex items-center justify-center gap-2 px-2 py-1.5 border-b bg-muted/30 shrink-0 relative z-panel titlebar-no-drag">
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={() => setCurrentPage(p => Math.max(1, p - 1))}
          disabled={currentPage <= 1}
        >
          <ChevronLeft className="h-4 w-4" />
        </Button>
        <span className="text-xs text-muted-foreground">
          {currentPage} / {numPages}
        </span>
        <Button
          variant="ghost"
          size="icon"
          className="h-7 w-7"
          onClick={() => setCurrentPage(p => Math.min(numPages, p + 1))}
          disabled={currentPage >= numPages}
        >
          <ChevronRight className="h-4 w-4" />
        </Button>
        <div className="w-px h-4 bg-border mx-1" />
        <Button variant="ghost" size="sm" className="h-7 text-xs" onClick={handleOpen}>
          <ExternalLink className="h-3 w-3 mr-1" />
          {t('打开')}
        </Button>
      </div>

      {/* PDF 渲染区 */}
      <div ref={containerRef} className="flex-1 overflow-auto flex justify-center p-4 bg-muted/20">
        <Document
          file={pdfFile}
          onLoadSuccess={({ numPages: n }) => setNumPages(n)}
          onLoadError={(err) => setError(err.message)}
          loading={
            <div className="flex items-center justify-center py-8">
              <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
            </div>
          }
        >
          <Page
            pageNumber={currentPage}
            width={containerWidth > 40 ? containerWidth - 40 : undefined}
            renderTextLayer={true}
            renderAnnotationLayer={true}
          />
        </Document>
      </div>
    </div>
  )
}

/**
 * Office 文档预览组件（docx/xlsx 转 HTML 渲染）
 */
function OfficePreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [htmlContent, setHtmlContent] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [canPreview, setCanPreview] = useState(false)

  const ext = file.name.split('.').pop()?.toLowerCase()

  useEffect(() => {
    const loadOffice = async () => {
      setIsLoading(true)
      setError(null)
      setHtmlContent(null)
      setCanPreview(false)

      try {
        if (ext === 'docx') {
          const { data } = await window.electronAPI.fm.readFileBase64(file.path)
          const buffer = Uint8Array.from(atob(data), c => c.charCodeAt(0))
          const mammoth = await import('mammoth')
          const result = await mammoth.convertToHtml({ arrayBuffer: buffer.buffer })
          setHtmlContent(result.value)
          setCanPreview(true)
        } else if (ext === 'xlsx' || ext === 'xls') {
          const { data } = await window.electronAPI.fm.readFileBase64(file.path)
          const buffer = Uint8Array.from(atob(data), c => c.charCodeAt(0))
          const XLSX = await import('xlsx')
          const workbook = XLSX.read(buffer, { type: 'array' })
          // 渲染第一个 sheet
          const firstSheet = workbook.Sheets[workbook.SheetNames[0]]
          if (firstSheet) {
            const html = XLSX.utils.sheet_to_html(firstSheet)
            setHtmlContent(html)
            setCanPreview(true)
          } else {
            setError(t('工作表为空'))
          }
        }
        // doc/ppt/pptx 不支持预览
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }

    if (ext === 'docx' || ext === 'xlsx' || ext === 'xls') {
      loadOffice()
    } else {
      setIsLoading(false)
    }
  }, [file.path, ext, t])

  const handleOpen = useCallback(() => {
    window.electronAPI.openFile(file.path)
  }, [file.path])

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }

  // 可预览的文件（docx/xlsx）
  if (canPreview && htmlContent) {
    const officeTypeLabel: Record<string, string> = {
      docx: 'Word',
      xlsx: 'Excel',
      xls: 'Excel',
    }
    return (
      <div className="h-full flex flex-col">
        <div className="flex items-center justify-between px-2 py-1.5 border-b bg-muted/30 shrink-0 relative z-panel">
          <span className="text-xs text-muted-foreground">{officeTypeLabel[ext || ''] || 'Office'} {t('预览')}</span>
          <Button variant="ghost" size="sm" className="h-6 text-xs titlebar-no-drag" onClick={handleOpen}>
            <ExternalLink className="h-3 w-3 mr-1" />
            {t('使用默认应用打开')}
          </Button>
        </div>
        <ScrollArea className="flex-1">
          <div
            className="p-4 prose prose-sm dark:prose-invert max-w-none [&_table]:border-collapse [&_td]:border [&_td]:border-border [&_td]:px-2 [&_td]:py-1 [&_th]:border [&_th]:border-border [&_th]:px-2 [&_th]:py-1 [&_th]:bg-muted/50"
            dangerouslySetInnerHTML={{ __html: htmlContent }}
          />
        </ScrollArea>
      </div>
    )
  }

  // 不可预览的文件（doc/ppt/pptx）或加载失败
  const officeIcons: Record<string, string> = {
    doc: 'Word',
    docx: 'Word',
    xls: 'Excel',
    xlsx: 'Excel',
    ppt: 'PowerPoint',
    pptx: 'PowerPoint',
  }

  return (
    <div className="h-full flex flex-col items-center justify-center p-4 text-center">
      <FileSpreadsheet className="h-16 w-16 text-muted-foreground mb-4" />
      <p className="text-sm font-medium mb-1">{file.name}</p>
      <p className="text-xs text-muted-foreground mb-4">{officeIcons[ext || ''] || 'Office'} {t('文档')}</p>
      {error && <p className="text-xs text-destructive mb-2">{error}</p>}
      <Button variant="outline" size="sm" onClick={handleOpen}>
        <ExternalLink className="h-3.5 w-3.5 mr-1.5" />
        {t('使用默认应用打开')}
      </Button>
    </div>
  )
}

/**
 * Unknown file type preview
 */
function UnknownPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  
  return (
    <div className="flex items-center justify-center h-full text-center p-4">
      <div>
        <File className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
        <p className="text-sm text-muted-foreground">{t('无法预览此文件类型')}</p>
        <p className="text-xs text-muted-foreground mt-1">
          {file.extension ? `.${file.extension}` : t('未知类型')}
        </p>
      </div>
    </div>
  )
}

/**
 * FilePreviewPanel component
 */
export function FilePreviewPanel({
  file,
  className,
  onConvertRequest,
  onClose,
}: FilePreviewPanelProps) {
  const t = useT()
  const [isEditing, setIsEditing] = useState(false)

  // 文件切换时退出编辑模式
  useEffect(() => {
    setIsEditing(false)
  }, [file?.path])

  const handleOpenExternal = useCallback(() => {
    if (file) {
      window.electronAPI.openFile(file.path)
    }
  }, [file])

  // Empty state
  if (!file) {
    return (
      <div className={cn("flex flex-col h-full scenic-opaque", className)}>
        <div className="flex items-center justify-center h-full text-center p-4">
          <div>
            <FileText className="h-12 w-12 text-muted-foreground/50 mx-auto mb-3" />
            <p className="text-sm text-muted-foreground">{t('选择文件以预览')}</p>
          </div>
        </div>
      </div>
    )
  }

  // Directory selected (shouldn't happen but handle gracefully)
  if (file.isDirectory) {
    return (
      <div className={cn("flex flex-col h-full scenic-opaque", className)}>
        <div className="flex items-center justify-center h-full text-center p-4">
          <div>
            <FileText className="h-12 w-12 text-muted-foreground/50 mx-auto mb-3" />
            <p className="text-sm text-muted-foreground">{t('选择文件以预览')}</p>
          </div>
        </div>
      </div>
    )
  }

  const fileType = getFileType(file)

  return (
    <div className={cn("flex flex-col h-full scenic-opaque", className)}>
      {/* Header */}
      <div className="flex items-center justify-between px-3 py-2 border-b shrink-0 relative z-panel">
        <div className="flex items-center gap-2 min-w-0 flex-1">
          {fileType === 'text' && <FileText className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'markdown' && <FileText className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'code' && <FileCode className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'image' && <ImageIcon className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'video' && <FileVideo className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'audio' && <FileAudio className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'pdf' && <FileType className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'office' && <FileSpreadsheet className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'unknown' && <File className="h-4 w-4 text-muted-foreground shrink-0" />}
          <span className="text-sm font-medium truncate">{file.name}</span>
        </div>
        <div className="flex items-center gap-1 shrink-0 titlebar-no-drag">
          {/* 可编辑文件且非编辑状态时显示编辑按钮 */}
          {isEditable(fileType) && !isEditing && (
            <Button
              variant="ghost"
              size="icon"
              className="h-6 w-6"
              onClick={() => setIsEditing(true)}
              title={t('编辑')}
            >
              <Pencil className="h-3.5 w-3.5" />
            </Button>
          )}
          <Button
            variant="ghost"
            size="icon"
            className="h-6 w-6"
            onClick={handleOpenExternal}
            title={t('使用默认应用打开')}
          >
            <ExternalLink className="h-3.5 w-3.5" />
          </Button>
          {onClose && (
            <Button
              variant="ghost"
              size="icon"
              className="h-6 w-6"
              onClick={onClose}
              title={t('关闭预览')}
            >
              <X className="h-3.5 w-3.5" />
            </Button>
          )}
        </div>
      </div>

      {/* Preview content */}
      <div className="flex-1 min-h-0">
        {fileType === 'markdown' && (
          <MarkdownPreview
            file={file}
            isEditing={isEditing}
            onEditingChange={setIsEditing}
            onConvertRequest={onConvertRequest}
          />
        )}
        {(fileType === 'text' || fileType === 'code') && (
          <TextPreviewEdit
            file={file}
            isEditing={isEditing}
            onEditingChange={setIsEditing}
          />
        )}
        {fileType === 'image' && <ImagePreview file={file} />}
        {fileType === 'video' && <VideoPreview file={file} />}
        {fileType === 'audio' && <AudioPreview file={file} />}
        {fileType === 'pdf' && <PdfPreview file={file} />}
        {fileType === 'office' && <OfficePreview file={file} />}
        {fileType === 'unknown' && <UnknownPreview file={file} />}
      </div>

      {/* File info footer */}
      <div className="px-3 py-2 border-t text-xs text-muted-foreground shrink-0 space-y-0.5">
        <div className="flex justify-between">
          <span>{t('大小')}</span>
          <span>{formatFileSize(file.size)}</span>
        </div>
        <div className="flex justify-between">
          <span>{t('修改时间')}</span>
          <span>{formatDate(file.modifiedTime)}</span>
        </div>
      </div>
    </div>
  )
}
