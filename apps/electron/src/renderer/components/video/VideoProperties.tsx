/**
 * VideoProperties - Properties panel for video configuration
 *
 * Features:
 * - Project name and description editing
 * - Video configuration (dimensions, fps, duration)
 * - Composition props editing
 * - Asset management
 *
 * @requirements 9.4
 */

import * as React from 'react';
import { useState, useCallback } from 'react';
import { Settings, Image, Music, Type, Palette, Clock, Maximize } from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { Separator } from '@/components/ui/separator';
import { useT } from '@/context/LocaleContext';
import type { VideoProject } from '@creator-flow/video';

export interface VideoPropertiesProps {
  /** Current video project */
  project: VideoProject | null;
  /** Callback when project is updated */
  onUpdate: (updates: Partial<VideoProject>) => void;
  /** Callback to open render dialog */
  onRender: () => void;
  /** Optional class name */
  className?: string;
}

/**
 * Section header component
 */
function SectionHeader({ icon: Icon, title }: { icon: React.ElementType; title: string }) {
  return (
    <div className="flex items-center gap-2 text-sm font-medium text-foreground mb-2">
      <Icon className="h-4 w-4 text-muted-foreground" />
      {title}
    </div>
  );
}

/**
 * Property row component
 */
function PropertyRow({
  label,
  children,
}: {
  label: string;
  children: React.ReactNode;
}) {
  return (
    <div className="flex items-center justify-between gap-2 py-1">
      <Label className="text-xs text-muted-foreground shrink-0">{label}</Label>
      <div className="flex-1 max-w-[140px]">{children}</div>
    </div>
  );
}

/**
 * VideoProperties component
 */
export function VideoProperties({
  project,
  onUpdate,
  onRender,
  className,
}: VideoPropertiesProps) {
  const t = useT();
  const [editingName, setEditingName] = useState(false);
  const [nameValue, setNameValue] = useState('');

  // Handle name edit start
  const handleNameEditStart = useCallback(() => {
    if (project) {
      setNameValue(project.name);
      setEditingName(true);
    }
  }, [project]);

  // Handle name edit submit
  const handleNameEditSubmit = useCallback(() => {
    if (nameValue.trim() && nameValue !== project?.name) {
      onUpdate({ name: nameValue.trim() });
    }
    setEditingName(false);
  }, [nameValue, project?.name, onUpdate]);

  // Handle config change
  const handleConfigChange = useCallback(
    (key: keyof VideoProject['config'], value: number) => {
      if (!project) return;
      onUpdate({
        config: {
          ...project.config,
          [key]: value,
        },
      });
    },
    [project, onUpdate]
  );

  // Handle description change
  const handleDescriptionChange = useCallback(
    (e: React.ChangeEvent<HTMLTextAreaElement>) => {
      onUpdate({ description: e.target.value });
    },
    [onUpdate]
  );

  // Preset dimensions
  const handlePresetChange = useCallback(
    (preset: string) => {
      if (!project) return;
      const presets: Record<string, { width: number; height: number }> = {
        '1080p': { width: 1920, height: 1080 },
        '720p': { width: 1280, height: 720 },
        '4k': { width: 3840, height: 2160 },
        'vertical': { width: 1080, height: 1920 },
        'square': { width: 1080, height: 1080 },
        'portrait': { width: 1080, height: 1350 },
      };
      const dimensions = presets[preset];
      if (dimensions) {
        onUpdate({
          config: {
            ...project.config,
            ...dimensions,
          },
        });
      }
    },
    [project, onUpdate]
  );

  if (!project) {
    return (
      <div className={cn('p-4 text-center text-muted-foreground text-sm', className)}>
        {t('选择一个项目查看属性')}
      </div>
    );
  }

  return (
    <div className={cn('p-3 space-y-4', className)}>
      {/* Project Info */}
      <div>
        <SectionHeader icon={Settings} title={t('项目信息')} />
        <div className="space-y-2">
          {/* Name */}
          <div>
            <Label className="text-xs text-muted-foreground">{t('名称')}</Label>
            {editingName ? (
              <Input
                value={nameValue}
                onChange={(e) => setNameValue(e.target.value)}
                onBlur={handleNameEditSubmit}
                onKeyDown={(e) => {
                  if (e.key === 'Enter') handleNameEditSubmit();
                  if (e.key === 'Escape') setEditingName(false);
                }}
                className="h-8 mt-1"
                autoFocus
              />
            ) : (
              <div
                className="text-sm font-medium cursor-pointer hover:bg-muted/50 rounded px-2 py-1 mt-1"
                onClick={handleNameEditStart}
              >
                {project.name}
              </div>
            )}
          </div>

          {/* Description */}
          <div>
            <Label className="text-xs text-muted-foreground">{t('描述')}</Label>
            <Textarea
              value={project.description || ''}
              onChange={handleDescriptionChange}
              placeholder={t('添加项目描述...')}
              className="h-16 mt-1 text-sm resize-none"
            />
          </div>
        </div>
      </div>

      <Separator />

      {/* Video Config */}
      <div>
        <SectionHeader icon={Maximize} title={t('视频配置')} />
        <div className="space-y-1">
          {/* Preset */}
          <PropertyRow label={t('预设')}>
            <Select onValueChange={handlePresetChange}>
              <SelectTrigger className="h-7 text-xs">
                <SelectValue placeholder={t('选择预设')} />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="1080p">1080p (16:9)</SelectItem>
                <SelectItem value="720p">720p (16:9)</SelectItem>
                <SelectItem value="4k">4K (16:9)</SelectItem>
                <SelectItem value="vertical">{t('竖屏')} (9:16)</SelectItem>
                <SelectItem value="square">{t('方形')} (1:1)</SelectItem>
                <SelectItem value="portrait">{t('肖像')} (4:5)</SelectItem>
              </SelectContent>
            </Select>
          </PropertyRow>

          {/* Width */}
          <PropertyRow label={t('宽度')}>
            <Input
              type="number"
              value={project.config.width}
              onChange={(e) => handleConfigChange('width', parseInt(e.target.value) || 1920)}
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* Height */}
          <PropertyRow label={t('高度')}>
            <Input
              type="number"
              value={project.config.height}
              onChange={(e) => handleConfigChange('height', parseInt(e.target.value) || 1080)}
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* FPS */}
          <PropertyRow label={t('帧率')}>
            <Select
              value={project.config.fps.toString()}
              onValueChange={(v) => handleConfigChange('fps', parseInt(v))}
            >
              <SelectTrigger className="h-7 text-xs">
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="24">24 fps</SelectItem>
                <SelectItem value="25">25 fps</SelectItem>
                <SelectItem value="30">30 fps</SelectItem>
                <SelectItem value="50">50 fps</SelectItem>
                <SelectItem value="60">60 fps</SelectItem>
              </SelectContent>
            </Select>
          </PropertyRow>

          {/* Duration */}
          <PropertyRow label={t('时长(帧)')}>
            <Input
              type="number"
              value={project.config.durationInFrames}
              onChange={(e) =>
                handleConfigChange('durationInFrames', parseInt(e.target.value) || 300)
              }
              className="h-7 text-xs"
            />
          </PropertyRow>

          {/* Duration in seconds (calculated) */}
          <PropertyRow label={t('时长(秒)')}>
            <div className="text-xs text-muted-foreground text-right">
              {(project.config.durationInFrames / project.config.fps).toFixed(2)}s
            </div>
          </PropertyRow>
        </div>
      </div>

      <Separator />

      {/* Assets */}
      <div>
        <SectionHeader icon={Image} title={t('素材')} />
        <div className="space-y-1">
          <div className="text-xs text-muted-foreground">
            {project.assets.length === 0
              ? t('暂无素材')
              : t('{{count}} 个素材', { count: project.assets.length })}
          </div>
          {project.assets.slice(0, 5).map((asset) => (
            <div
              key={asset.id}
              className="flex items-center gap-2 text-xs py-1 px-2 rounded hover:bg-muted/50"
            >
              {asset.type === 'image' && <Image className="h-3 w-3" />}
              {asset.type === 'audio' && <Music className="h-3 w-3" />}
              {asset.type === 'font' && <Type className="h-3 w-3" />}
              <span className="truncate flex-1">{asset.name}</span>
            </div>
          ))}
          {project.assets.length > 5 && (
            <div className="text-xs text-muted-foreground text-center">
              +{project.assets.length - 5} {t('更多')}
            </div>
          )}
        </div>
      </div>

      <Separator />

      {/* Render Button */}
      <Button onClick={onRender} className="w-full">
        {t('导出视频')}
      </Button>
    </div>
  );
}
