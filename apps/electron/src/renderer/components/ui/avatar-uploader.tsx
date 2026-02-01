/**
 * AvatarUploader
 *
 * 头像上传组件，支持图片选择、裁剪、旋转和上传
 */

import { useState, useRef, useCallback, useEffect } from 'react'
import ReactCrop, { type Crop, centerCrop, makeAspectCrop } from 'react-image-crop'
import 'react-image-crop/dist/ReactCrop.css'
import { Camera, RotateCw, Check, X, Upload, Loader2 } from 'lucide-react'
import { Button } from './button'
import { Dialog, DialogContent, DialogHeader, DialogTitle } from './dialog'
import { useT } from '@/context/LocaleContext'
import { cn } from '@/lib/utils'

interface AvatarUploaderProps {
  value?: string
  onChange?: (url: string) => void
  onUpload?: (file: File) => Promise<{ url: string }>
  size?: number
  className?: string
  disabled?: boolean
}

function centerAspectCrop(mediaWidth: number, mediaHeight: number, aspect: number) {
  return centerCrop(
    makeAspectCrop(
      {
        unit: '%',
        width: 90,
      },
      aspect,
      mediaWidth,
      mediaHeight
    ),
    mediaWidth,
    mediaHeight
  )
}

export function AvatarUploader({
  value,
  onChange,
  onUpload,
  size = 96,
  className,
  disabled = false,
}: AvatarUploaderProps) {
  const t = useT()
  const [isOpen, setIsOpen] = useState(false)
  const [imgSrc, setImgSrc] = useState('')
  const [crop, setCrop] = useState<Crop>()
  const [rotation, setRotation] = useState(0)
  const [isUploading, setIsUploading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  
  const imgRef = useRef<HTMLImageElement>(null)
  const fileInputRef = useRef<HTMLInputElement>(null)
  const canvasRef = useRef<HTMLCanvasElement>(null)

  // 处理文件选择
  const handleFileSelect = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      const file = e.target.files[0]
      
      // 验证文件类型
      if (!file.type.startsWith('image/')) {
        setError(t('请选择图片文件'))
        return
      }
      
      // 验证文件大小 (5MB)
      if (file.size > 5 * 1024 * 1024) {
        setError(t('文件大小不能超过 5MB'))
        return
      }
      
      setError(null)
      const reader = new FileReader()
      reader.addEventListener('load', () => {
        setImgSrc(reader.result?.toString() || '')
        setRotation(0)
        setCrop(undefined)
        setIsOpen(true)
      })
      reader.readAsDataURL(file)
    }
    // 清除 input value 以便再次选择同一文件
    e.target.value = ''
  }, [t])

  // 图片加载完成后设置默认裁剪区域
  const onImageLoad = useCallback((e: React.SyntheticEvent<HTMLImageElement>) => {
    const { width, height } = e.currentTarget
    setCrop(centerAspectCrop(width, height, 1))
  }, [])

  // 旋转图片
  const handleRotate = useCallback(() => {
    setRotation((prev) => (prev + 90) % 360)
  }, [])

  // 获取裁剪后的图片 Blob
  const getCroppedImage = useCallback(async (): Promise<Blob | null> => {
    if (!imgRef.current || !crop || !canvasRef.current) return null

    const image = imgRef.current
    const canvas = canvasRef.current
    const ctx = canvas.getContext('2d')
    if (!ctx) return null

    // 设置输出尺寸
    const outputSize = 256
    canvas.width = outputSize
    canvas.height = outputSize

    // 计算缩放比例
    const scaleX = image.naturalWidth / image.width
    const scaleY = image.naturalHeight / image.height

    // 裁剪区域
    const cropX = crop.x * scaleX
    const cropY = crop.y * scaleY
    const cropWidth = crop.width * scaleX
    const cropHeight = crop.height * scaleY

    // 清除画布
    ctx.fillStyle = '#ffffff'
    ctx.fillRect(0, 0, outputSize, outputSize)

    // 应用旋转
    ctx.save()
    ctx.translate(outputSize / 2, outputSize / 2)
    ctx.rotate((rotation * Math.PI) / 180)
    ctx.translate(-outputSize / 2, -outputSize / 2)

    // 绘制裁剪后的图片
    ctx.drawImage(
      image,
      cropX,
      cropY,
      cropWidth,
      cropHeight,
      0,
      0,
      outputSize,
      outputSize
    )
    ctx.restore()

    // 转换为 Blob
    return new Promise((resolve) => {
      canvas.toBlob(
        (blob) => resolve(blob),
        'image/jpeg',
        0.9
      )
    })
  }, [crop, rotation])

  // 确认并上传
  const handleConfirm = useCallback(async () => {
    if (!onUpload) return

    try {
      setIsUploading(true)
      setError(null)

      const blob = await getCroppedImage()
      if (!blob) {
        setError(t('裁剪图片失败'))
        return
      }

      // 创建 File 对象
      const file = new File([blob], 'avatar.jpg', { type: 'image/jpeg' })

      // 上传
      const result = await onUpload(file)
      onChange?.(result.url)
      setIsOpen(false)
      setImgSrc('')
    } catch (err: any) {
      console.error('Upload avatar failed:', err)
      setError(err.message || t('上传失败'))
    } finally {
      setIsUploading(false)
    }
  }, [getCroppedImage, onUpload, onChange, t])

  // 取消
  const handleCancel = useCallback(() => {
    setIsOpen(false)
    setImgSrc('')
    setError(null)
  }, [])

  // 触发文件选择
  const triggerFileSelect = useCallback(() => {
    if (!disabled) {
      fileInputRef.current?.click()
    }
  }, [disabled])

  return (
    <>
      {/* 隐藏的文件输入 */}
      <input
        ref={fileInputRef}
        type="file"
        accept="image/jpeg,image/png,image/gif,image/webp"
        className="hidden"
        onChange={handleFileSelect}
      />

      {/* 隐藏的 Canvas 用于裁剪 */}
      <canvas ref={canvasRef} className="hidden" />

      {/* 头像预览和上传按钮 */}
      <div className={cn('flex items-center gap-4', className)}>
        <div
          className={cn(
            'relative rounded-full overflow-hidden bg-muted cursor-pointer group',
            disabled && 'cursor-not-allowed opacity-50'
          )}
          style={{ width: size, height: size }}
          onClick={triggerFileSelect}
        >
          {value ? (
            <img
              src={value}
              alt="Avatar"
              className="w-full h-full object-cover"
            />
          ) : (
            <div className="w-full h-full flex items-center justify-center text-muted-foreground">
              <Camera className="size-6" />
            </div>
          )}
          {!disabled && (
            <div className="absolute inset-0 bg-black/50 opacity-0 group-hover:opacity-100 transition-opacity flex items-center justify-center">
              <Camera className="size-6 text-white" />
            </div>
          )}
        </div>
        <Button
          variant="outline"
          size="sm"
          onClick={triggerFileSelect}
          disabled={disabled}
        >
          <Upload className="size-4 mr-1.5" />
          {t('上传头像')}
        </Button>
      </div>

      {/* 裁剪对话框 */}
      <Dialog open={isOpen} onOpenChange={(open) => !open && handleCancel()}>
        <DialogContent className="max-w-md">
          <DialogHeader>
            <DialogTitle>{t('裁剪头像')}</DialogTitle>
          </DialogHeader>

          <div className="space-y-4">
            {error && (
              <div className="text-sm text-destructive bg-destructive/10 px-3 py-2 rounded">
                {error}
              </div>
            )}

            {imgSrc && (
              <div
                className="relative bg-muted rounded-lg overflow-hidden"
                style={{
                  transform: `rotate(${rotation}deg)`,
                  transformOrigin: 'center center',
                }}
              >
                <ReactCrop
                  crop={crop}
                  onChange={(c) => setCrop(c)}
                  aspect={1}
                  circularCrop
                >
                  <img
                    ref={imgRef}
                    src={imgSrc}
                    alt="Crop preview"
                    onLoad={onImageLoad}
                    className="max-h-[300px] w-auto mx-auto"
                  />
                </ReactCrop>
              </div>
            )}

            <div className="flex items-center justify-between">
              <Button
                variant="outline"
                size="sm"
                onClick={handleRotate}
                disabled={isUploading}
              >
                <RotateCw className="size-4 mr-1.5" />
                {t('旋转')}
              </Button>

              <div className="flex gap-2">
                <Button
                  variant="outline"
                  size="sm"
                  onClick={handleCancel}
                  disabled={isUploading}
                >
                  <X className="size-4 mr-1.5" />
                  {t('取消')}
                </Button>
                <Button
                  size="sm"
                  onClick={handleConfirm}
                  disabled={isUploading || !crop}
                >
                  {isUploading ? (
                    <Loader2 className="size-4 mr-1.5 animate-spin" />
                  ) : (
                    <Check className="size-4 mr-1.5" />
                  )}
                  {isUploading ? t('上传中...') : t('确认')}
                </Button>
              </div>
            </div>
          </div>
        </DialogContent>
      </Dialog>
    </>
  )
}
