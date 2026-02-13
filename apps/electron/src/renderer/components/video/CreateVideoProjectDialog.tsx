/**
 * CreateVideoProjectDialog - 从模板创建视频项目的对话框
 */

import { useState, useEffect } from 'react';
import { useT } from '@/context/LocaleContext';
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
} from '@/components/ui/dialog';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import { Badge } from '@/components/ui/badge';
import type { VideoTemplate } from '@sprouty-ai/video';

export interface CreateVideoProjectDialogProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  template: VideoTemplate | null;
  onConfirm: (data: { name: string; description: string }) => void;
}

export function CreateVideoProjectDialog({
  open,
  onOpenChange,
  template,
  onConfirm,
}: CreateVideoProjectDialogProps) {
  const t = useT();
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');

  // 当模板变化或对话框打开时，预填模板信息
  useEffect(() => {
    if (open && template) {
      setName(`${template.name} - ${new Date().toLocaleDateString()}`);
      setDescription(template.description || '');
    }
  }, [open, template]);

  const handleClose = (v: boolean) => {
    if (!v) {
      setName('');
      setDescription('');
    }
    onOpenChange(v);
  };

  const handleConfirm = () => {
    if (!name.trim()) return;
    onConfirm({ name: name.trim(), description: description.trim() });
    handleClose(false);
  };

  if (!template) return null;

  return (
    <Dialog open={open} onOpenChange={handleClose}>
      <DialogContent className="sm:max-w-md">
        <DialogHeader>
          <DialogTitle>{t('从模板创建项目')}</DialogTitle>
          <DialogDescription>
            {t('填写项目信息，确认后将基于模板创建新的视频项目')}
          </DialogDescription>
        </DialogHeader>

        {/* 模板信息预览 */}
        <div className="flex items-center gap-3 rounded-lg border bg-muted/30 p-3">
          <div className="flex-1 min-w-0">
            <div className="text-sm font-medium truncate">{t(template.name)}</div>
            <div className="flex items-center gap-2 mt-1 text-xs text-muted-foreground">
              <span>{template.defaultConfig.width}×{template.defaultConfig.height}</span>
              <span>·</span>
              <span>{template.defaultConfig.fps}fps</span>
              <span>·</span>
              <span>{(template.defaultConfig.durationInFrames / template.defaultConfig.fps).toFixed(1)}s</span>
            </div>
          </div>
          {template.aspectRatio && (
            <Badge variant="secondary" className="text-[10px] shrink-0">
              {template.aspectRatio}
            </Badge>
          )}
        </div>

        <div className="space-y-4">
          <div className="space-y-2">
            <Label>{t('项目名称')} *</Label>
            <Input
              value={name}
              onChange={(e) => setName(e.target.value)}
              placeholder={t('输入项目名称')}
              autoFocus
              onKeyDown={(e) => {
                if (e.key === 'Enter' && name.trim()) handleConfirm();
              }}
            />
          </div>

          <div className="space-y-2">
            <Label>{t('项目描述')}</Label>
            <Textarea
              value={description}
              onChange={(e) => setDescription(e.target.value)}
              placeholder={t('简要描述项目内容...')}
              rows={3}
            />
          </div>
        </div>

        <DialogFooter>
          <Button variant="outline" onClick={() => handleClose(false)}>
            {t('取消')}
          </Button>
          <Button onClick={handleConfirm} disabled={!name.trim()}>
            {t('创建项目')}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
}
