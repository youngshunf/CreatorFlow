/**
 * VideoTemplatePicker - 精简版视频模板选择器
 *
 * 用于 CreateContentDialog 中的横向滚动模板卡片列表。
 * 不使用 Remotion Player（太重），用色块+文字代替缩略图。
 */

import { useT } from '@/context/LocaleContext'
import { cn } from '@/lib/utils'
import { Check } from 'lucide-react'
import { ALL_TEMPLATES, type VideoTemplate } from '@sprouty-ai/video'

interface VideoTemplatePickerProps {
  selected?: string
  onSelect: (templateId: string) => void
  className?: string
}

/** 根据模板默认配色生成缩略图背景 */
function TemplateThumbnail({ template }: { template: VideoTemplate }) {
  const colors = template.defaultProps.colors
  const bg = colors?.background ?? '#1a1a2e'
  const primary = colors?.primary ?? '#6366f1'
  const text = colors?.text ?? '#ffffff'

  return (
    <div
      className="w-full h-full flex flex-col items-center justify-center gap-1 p-2"
      style={{ background: `linear-gradient(135deg, ${bg} 0%, ${primary}30 100%)` }}
    >
      <div
        className="text-[10px] font-bold leading-tight text-center line-clamp-2"
        style={{ color: text }}
      >
        {template.name}
      </div>
      <div
        className="text-[8px] opacity-60"
        style={{ color: text }}
      >
        {template.aspectRatio ?? `${template.defaultConfig.width}x${template.defaultConfig.height}`}
      </div>
    </div>
  )
}

export function VideoTemplatePicker({ selected, onSelect, className }: VideoTemplatePickerProps) {
  const t = useT()

  return (
    <div className={cn('space-y-2', className)}>
      {/* 横向滚动容器 */}
      <div className="flex gap-2 overflow-x-auto pb-1 -mx-1 px-1">
        {ALL_TEMPLATES.map((tpl) => {
          const isSelected = selected === tpl.id
          return (
            <button
              key={tpl.id}
              type="button"
              onClick={() => onSelect(tpl.id)}
              className={cn(
                'relative flex-shrink-0 w-28 rounded-md border overflow-hidden transition-all',
                'hover:border-primary/50 hover:shadow-sm',
                isSelected
                  ? 'border-primary ring-1 ring-primary'
                  : 'border-border'
              )}
            >
              {/* 缩略图区域 */}
              <div className="h-16">
                <TemplateThumbnail template={tpl} />
              </div>
              {/* 名称 + 分辨率 */}
              <div className="px-1.5 py-1 bg-background">
                <div className="text-[10px] font-medium truncate">{t(tpl.name)}</div>
                <div className="text-[9px] text-muted-foreground">
                  {tpl.defaultConfig.width}x{tpl.defaultConfig.height}
                </div>
              </div>
              {/* 选中标记 */}
              {isSelected && (
                <div className="absolute top-1 right-1 h-4 w-4 rounded-full bg-primary flex items-center justify-center">
                  <Check className="h-2.5 w-2.5 text-primary-foreground" />
                </div>
              )}
            </button>
          )
        })}
      </div>
    </div>
  )
}
