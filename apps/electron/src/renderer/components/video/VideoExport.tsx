/**
 * VideoExport - Export options component
 *
 * Features:
 * - Output format selection (MP4, WebM, GIF)
 * - Quality preset selection
 * - Composition selection
 * - Export progress display
 *
 * @requirements 9.7
 */

import * as React from 'react';
import { useState, useEffect, useCallback } from 'react';
import { Download, Film, Settings, Loader2, CheckCircle, XCircle } from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Label } from '@/components/ui/label';
import { Progress } from '@/components/ui/progress';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { Separator } from '@/components/ui/separator';
import { useT } from '@/context/LocaleContext';
import type { VideoProject, OutputFormat, QualityPreset, RenderProgress } from '@creator-flow/video';

export interface VideoExportProps {
  /** Current video project */
  project: VideoProject;
  /** Callback when export is triggered */
  onExport: (options: {
    compositionId: string;
    outputFormat: OutputFormat;
    quality: QualityPreset;
  }) => void;
  /** Callback to cancel/close export panel */
  onCancel: () => void;
  /** Optional class name */
  className?: string;
}

/**
 * Format descriptions
 */
const FORMAT_INFO: Record<OutputFormat, { name: string; description: string }> = {
  mp4: { name: 'MP4', description: '通用格式，兼容性最好' },
  webm: { name: 'WebM', description: '网页优化，文件较小' },
  gif: { name: 'GIF', description: '动图格式，无声音' },
};

/**
 * Quality descriptions
 */
const QUALITY_INFO: Record<QualityPreset, { name: string; description: string }> = {
  draft: { name: '草稿', description: '快速预览，低质量' },
  standard: { name: '标准', description: '平衡质量和大小' },
  high: { name: '高质量', description: '最佳质量，文件较大' },
};

/**
 * VideoExport component
 */
export function VideoExport({
  project,
  onExport,
  onCancel,
  className,
}: VideoExportProps) {
  const t = useT();

  // Export options
  const [compositionId, setCompositionId] = useState<string>(
    project.compositions[0]?.id || ''
  );
  const [outputFormat, setOutputFormat] = useState<OutputFormat>('mp4');
  const [quality, setQuality] = useState<QualityPreset>('standard');

  // Render progress
  const [isRendering, setIsRendering] = useState(false);
  const [renderProgress, setRenderProgress] = useState<RenderProgress | null>(null);

  // Subscribe to render progress
  useEffect(() => {
    const unsubscribe = window.electronAPI.video.onRenderProgress((progress) => {
      setRenderProgress(progress);
      if (progress.status === 'completed' || progress.status === 'failed') {
        setIsRendering(false);
      }
    });

    return () => {
      unsubscribe();
    };
  }, []);

  // Handle export
  const handleExport = useCallback(() => {
    if (!compositionId) return;
    setIsRendering(true);
    setRenderProgress({ status: 'bundling', progress: 0 });
    onExport({ compositionId, outputFormat, quality });
  }, [compositionId, outputFormat, quality, onExport]);

  // Handle cancel render
  const handleCancelRender = useCallback(() => {
    window.electronAPI.video.cancelRender();
    setIsRendering(false);
    setRenderProgress(null);
  }, []);

  // Get status text
  const getStatusText = (status: RenderProgress['status']): string => {
    switch (status) {
      case 'bundling':
        return t('打包中...');
      case 'preparing':
        return t('准备中...');
      case 'rendering':
        return t('渲染中...');
      case 'completed':
        return t('完成');
      case 'failed':
        return t('失败');
      default:
        return '';
    }
  };

  return (
    <div className={cn('p-3 space-y-4', className)}>
      {/* Header */}
      <div className="flex items-center gap-2">
        <Download className="h-5 w-5 text-primary" />
        <h3 className="font-medium">{t('导出视频')}</h3>
      </div>

      {/* Composition selection */}
      {project.compositions.length > 1 && (
        <div>
          <Label className="text-xs text-muted-foreground mb-1.5 block">
            {t('选择组件')}
          </Label>
          <Select value={compositionId} onValueChange={setCompositionId}>
            <SelectTrigger className="h-9">
              <SelectValue />
            </SelectTrigger>
            <SelectContent>
              {project.compositions.map((comp) => (
                <SelectItem key={comp.id} value={comp.id}>
                  <div className="flex items-center gap-2">
                    <Film className="h-4 w-4" />
                    {comp.name}
                  </div>
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
      )}

      <Separator />

      {/* Format selection */}
      <div>
        <Label className="text-xs text-muted-foreground mb-1.5 block">
          {t('输出格式')}
        </Label>
        <div className="grid grid-cols-3 gap-2">
          {(Object.keys(FORMAT_INFO) as OutputFormat[]).map((format) => (
            <button
              key={format}
              className={cn(
                'p-2 rounded-lg border text-center transition-colors',
                outputFormat === format
                  ? 'border-primary bg-primary/10'
                  : 'border-border hover:border-primary/50'
              )}
              onClick={() => setOutputFormat(format)}
            >
              <div className="font-medium text-sm">{FORMAT_INFO[format].name}</div>
              <div className="text-[10px] text-muted-foreground mt-0.5">
                {t(FORMAT_INFO[format].description)}
              </div>
            </button>
          ))}
        </div>
      </div>

      {/* Quality selection */}
      <div>
        <Label className="text-xs text-muted-foreground mb-1.5 block">
          {t('质量预设')}
        </Label>
        <div className="grid grid-cols-3 gap-2">
          {(Object.keys(QUALITY_INFO) as QualityPreset[]).map((q) => (
            <button
              key={q}
              className={cn(
                'p-2 rounded-lg border text-center transition-colors',
                quality === q
                  ? 'border-primary bg-primary/10'
                  : 'border-border hover:border-primary/50'
              )}
              onClick={() => setQuality(q)}
            >
              <div className="font-medium text-sm">{t(QUALITY_INFO[q].name)}</div>
              <div className="text-[10px] text-muted-foreground mt-0.5">
                {t(QUALITY_INFO[q].description)}
              </div>
            </button>
          ))}
        </div>
      </div>

      <Separator />

      {/* Export info */}
      <div className="text-xs text-muted-foreground space-y-1">
        <div className="flex justify-between">
          <span>{t('分辨率')}</span>
          <span>
            {project.config.width}×{project.config.height}
          </span>
        </div>
        <div className="flex justify-between">
          <span>{t('帧率')}</span>
          <span>{project.config.fps} fps</span>
        </div>
        <div className="flex justify-between">
          <span>{t('时长')}</span>
          <span>
            {(project.config.durationInFrames / project.config.fps).toFixed(1)}s
          </span>
        </div>
      </div>

      {/* Render progress */}
      {renderProgress && (
        <div className="space-y-2">
          <div className="flex items-center gap-2">
            {renderProgress.status === 'completed' ? (
              <CheckCircle className="h-4 w-4 text-green-500" />
            ) : renderProgress.status === 'failed' ? (
              <XCircle className="h-4 w-4 text-destructive" />
            ) : (
              <Loader2 className="h-4 w-4 animate-spin" />
            )}
            <span className="text-sm">{getStatusText(renderProgress.status)}</span>
            <span className="text-sm text-muted-foreground ml-auto">
              {renderProgress.progress.toFixed(0)}%
            </span>
          </div>
          <Progress value={renderProgress.progress} className="h-2" />
          {renderProgress.error && (
            <p className="text-xs text-destructive">{renderProgress.error}</p>
          )}
        </div>
      )}

      {/* Actions */}
      <div className="flex gap-2">
        <Button variant="outline" className="flex-1" onClick={onCancel}>
          {t('取消')}
        </Button>
        {isRendering ? (
          <Button
            variant="destructive"
            className="flex-1"
            onClick={handleCancelRender}
          >
            {t('停止渲染')}
          </Button>
        ) : (
          <Button
            className="flex-1"
            onClick={handleExport}
            disabled={!compositionId}
          >
            <Download className="h-4 w-4 mr-1" />
            {t('开始导出')}
          </Button>
        )}
      </div>
    </div>
  );
}
