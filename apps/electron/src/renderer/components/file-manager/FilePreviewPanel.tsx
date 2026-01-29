/**
 * FilePreviewPanel - Preview panel for selected files
 *
 * Features:
 * - Text file preview with syntax highlighting
 * - Image preview
 * - File metadata display
 * - Open in external app button
 */

import * as React from 'react'
import { useState, useEffect, useCallback } from 'react'
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
  ZoomIn,
  ZoomOut,
  RotateCw,
  Maximize2,
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Button } from '@/components/ui/button'
import { useT } from '@/context/LocaleContext'
import type { FMFileEntry } from '../../../shared/types'

export interface FilePreviewPanelProps {
  /** Selected file to preview */
  file?: FMFileEntry | null
  className?: string
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
function getFileType(file: FMFileEntry): 'text' | 'image' | 'code' | 'video' | 'audio' | 'pdf' | 'unknown' {
  const ext = file.name.split('.').pop()?.toLowerCase()
  
  // Text types
  if (['txt', 'md', 'markdown', 'log', 'csv', 'ini', 'conf', 'cfg'].includes(ext || '')) {
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
  
  return 'unknown'
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
 * Text/Code preview component
 */
function TextPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [content, setContent] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [isTruncated, setIsTruncated] = useState(false)
  
  const MAX_PREVIEW_SIZE = 100 * 1024 // 100KB limit for preview
  
  useEffect(() => {
    const loadContent = async () => {
      setIsLoading(true)
      setError(null)
      setIsTruncated(false)
      
      try {
        // Check file size first
        if (file.size && file.size > MAX_PREVIEW_SIZE) {
          setIsTruncated(true)
        }
        
        // Read file content via IPC
        const text = await window.electronAPI.readFile(file.path)
        
        // Truncate if too large
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
  
  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }
  
  if (error) {
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
      <div className="flex items-center justify-center gap-1 px-2 py-1.5 border-b bg-muted/30 shrink-0">
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
  const [videoSrc, setVideoSrc] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  
  useEffect(() => {
    const loadVideo = async () => {
      setIsLoading(true)
      setError(null)
      
      try {
        // Load video via IPC as base64
        const { data, mimeType } = await window.electronAPI.fm.readFileBase64(file.path, 50 * 1024 * 1024) // 50MB limit
        const src = `data:${mimeType};base64,${data}`
        setVideoSrc(src)
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }
    
    loadVideo()
  }, [file.path])
  
  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }
  
  if (error || !videoSrc) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileVideo className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此视频')}</p>
          {error && <p className="text-xs text-destructive mt-1">{error}</p>}
        </div>
      </div>
    )
  }
  
  return (
    <div className="h-full flex items-center justify-center p-4 bg-black/5">
      <video
        src={videoSrc}
        controls
        className="max-w-full max-h-full rounded shadow-sm"
        onError={() => setError('Failed to load video')}
      />
    </div>
  )
}

/**
 * Audio preview component
 */
function AudioPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  const [audioSrc, setAudioSrc] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  
  useEffect(() => {
    const loadAudio = async () => {
      setIsLoading(true)
      setError(null)
      
      try {
        // Load audio via IPC as base64
        const { data, mimeType } = await window.electronAPI.fm.readFileBase64(file.path, 20 * 1024 * 1024) // 20MB limit
        const src = `data:${mimeType};base64,${data}`
        setAudioSrc(src)
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err))
      } finally {
        setIsLoading(false)
      }
    }
    
    loadAudio()
  }, [file.path])
  
  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <RefreshCw className="h-5 w-5 animate-spin text-muted-foreground" />
      </div>
    )
  }
  
  if (error || !audioSrc) {
    return (
      <div className="flex items-center justify-center h-full text-center p-4">
        <div>
          <FileAudio className="h-10 w-10 text-muted-foreground mx-auto mb-2" />
          <p className="text-sm text-muted-foreground">{t('无法预览此音频')}</p>
          {error && <p className="text-xs text-destructive mt-1">{error}</p>}
        </div>
      </div>
    )
  }
  
  return (
    <div className="h-full flex flex-col items-center justify-center p-4">
      <FileAudio className="h-16 w-16 text-muted-foreground mb-4" />
      <audio
        src={audioSrc}
        controls
        className="w-full max-w-md"
        onError={() => setError('Failed to load audio')}
      />
    </div>
  )
}

/**
 * PDF preview component (shows info and open button since PDF rendering is complex)
 */
function PdfPreview({ file }: { file: FMFileEntry }) {
  const t = useT()
  
  const handleOpen = useCallback(() => {
    window.electronAPI.openFile(file.path)
  }, [file.path])
  
  return (
    <div className="h-full flex flex-col items-center justify-center p-4 text-center">
      <FileType className="h-16 w-16 text-muted-foreground mb-4" />
      <p className="text-sm font-medium mb-1">{file.name}</p>
      <p className="text-xs text-muted-foreground mb-4">{t('PDF 文档')}</p>
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
}: FilePreviewPanelProps) {
  const t = useT()
  
  const handleOpenExternal = useCallback(() => {
    if (file) {
      window.electronAPI.openFile(file.path)
    }
  }, [file])
  
  // Empty state
  if (!file) {
    return (
      <div className={cn("flex flex-col h-full", className)}>
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
      <div className={cn("flex flex-col h-full", className)}>
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
    <div className={cn("flex flex-col h-full", className)}>
      {/* Header */}
      <div className="flex items-center justify-between px-3 py-2 border-b shrink-0">
      <div className="flex items-center gap-2 min-w-0 flex-1">
          {fileType === 'text' && <FileText className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'code' && <FileCode className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'image' && <ImageIcon className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'video' && <FileVideo className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'audio' && <FileAudio className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'pdf' && <FileType className="h-4 w-4 text-muted-foreground shrink-0" />}
          {fileType === 'unknown' && <File className="h-4 w-4 text-muted-foreground shrink-0" />}
          <span className="text-sm font-medium truncate">{file.name}</span>
        </div>
        <Button
          variant="ghost"
          size="icon"
          className="h-6 w-6 shrink-0"
          onClick={handleOpenExternal}
          title={t('使用默认应用打开')}
        >
          <ExternalLink className="h-3.5 w-3.5" />
        </Button>
      </div>
      
      {/* Preview content */}
      <div className="flex-1 min-h-0">
        {(fileType === 'text' || fileType === 'code') && <TextPreview file={file} />}
        {fileType === 'image' && <ImagePreview file={file} />}
        {fileType === 'video' && <VideoPreview file={file} />}
        {fileType === 'audio' && <AudioPreview file={file} />}
        {fileType === 'pdf' && <PdfPreview file={file} />}
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
