/**
 * TransitionEditor - 过渡效果编辑器
 *
 * 功能：
 * - 过渡类型选择（grid 按钮组）
 * - 时长调节（Slider + 数字显示）
 * - 方向选择（仅 slide/wipe/flip 时显示）
 *
 * TransitionIndicator - 场景列表中的过渡效果指示条
 */

import * as React from 'react';
import { useCallback } from 'react';
import {
  Blend,
  ArrowRightLeft,
  Columns,
  RotateCcw,
  Clock,
  Ban,
  ArrowLeft,
  ArrowRight,
  ArrowUp,
  ArrowDown,
  ChevronRight,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Label } from '@/components/ui/label';
import { Slider } from '@/components/ui/slider';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { useT } from '@/context/LocaleContext';
import type {
  Transition,
  TransitionType,
  TransitionDirection,
} from '@sprouty-ai/video';

// ============================================================================
// 常量
// ============================================================================

/** 过渡类型配置 */
const TRANSITION_TYPES: {
  value: TransitionType;
  labelKey: string;
  icon: React.ElementType;
}[] = [
  { value: 'fade', labelKey: '淡入淡出', icon: Blend },
  { value: 'slide', labelKey: '滑动', icon: ArrowRightLeft },
  { value: 'wipe', labelKey: '擦除', icon: Columns },
  { value: 'flip', labelKey: '翻转', icon: RotateCcw },
  { value: 'clock-wipe', labelKey: '时钟擦除', icon: Clock },
  { value: 'none', labelKey: '无', icon: Ban },
];

/** 需要方向选择的过渡类型 */
const DIRECTIONAL_TYPES: TransitionType[] = ['slide', 'wipe', 'flip'];

/** 方向选项配置 */
const DIRECTION_OPTIONS: {
  value: TransitionDirection;
  labelKey: string;
  icon: React.ElementType;
}[] = [
  { value: 'from-left', labelKey: '从左', icon: ArrowRight },
  { value: 'from-right', labelKey: '从右', icon: ArrowLeft },
  { value: 'from-top', labelKey: '从上', icon: ArrowDown },
  { value: 'from-bottom', labelKey: '从下', icon: ArrowUp },
];

/** 过渡类型对应的简短图标（用于 Indicator） */
const TYPE_ICON_MAP: Record<TransitionType, React.ElementType> = {
  fade: Blend,
  slide: ArrowRightLeft,
  wipe: Columns,
  flip: RotateCcw,
  'clock-wipe': Clock,
  none: Ban,
};

// ============================================================================
// TransitionEditor
// ============================================================================

export interface TransitionEditorProps {
  /** 当前过渡效果 */
  transition: Transition;
  /** 更新回调 */
  onUpdate: (updates: Partial<Transition>) => void;
  /** 可选 class name */
  className?: string;
}

/**
 * 过渡效果编辑器组件
 */
export function TransitionEditor({
  transition,
  onUpdate,
  className,
}: TransitionEditorProps) {
  const t = useT();

  // 是否显示方向选择
  const showDirection = DIRECTIONAL_TYPES.includes(transition.type);

  // 切换过渡类型
  const handleTypeChange = useCallback(
    (type: TransitionType) => {
      const updates: Partial<Transition> = { type };
      // 切换到有方向的类型时，设置默认方向
      if (DIRECTIONAL_TYPES.includes(type) && !transition.direction) {
        updates.direction = 'from-right';
      }
      onUpdate(updates);
    },
    [transition.direction, onUpdate],
  );

  // 调整时长
  const handleDurationChange = useCallback(
    ([value]: number[]) => {
      onUpdate({ durationInFrames: value });
    },
    [onUpdate],
  );

  // 切换方向
  const handleDirectionChange = useCallback(
    (direction: string) => {
      onUpdate({ direction: direction as TransitionDirection });
    },
    [onUpdate],
  );

  return (
    <div className={cn('space-y-3', className)}>
      {/* 过渡类型选择 */}
      <div>
        <Label className="text-xs text-muted-foreground mb-1.5 block">
          {t('过渡类型')}
        </Label>
        <div className="grid grid-cols-3 gap-1.5">
          {TRANSITION_TYPES.map(({ value, labelKey, icon: Icon }) => (
            <Button
              key={value}
              variant={transition.type === value ? 'default' : 'outline'}
              size="sm"
              className="h-8 text-xs gap-1 px-2"
              onClick={() => handleTypeChange(value)}
            >
              <Icon className="h-3.5 w-3.5" />
              {t(labelKey)}
            </Button>
          ))}
        </div>
      </div>

      {/* 时长调节 */}
      <div>
        <div className="flex items-center justify-between mb-1.5">
          <Label className="text-xs text-muted-foreground">
            {t('过渡时长')}
          </Label>
          <span className="text-xs text-muted-foreground tabular-nums">
            {transition.durationInFrames} {t('帧')}
          </span>
        </div>
        <Slider
          value={[transition.durationInFrames]}
          onValueChange={handleDurationChange}
          min={5}
          max={60}
          step={1}
        />
      </div>

      {/* 方向选择（仅 slide/wipe/flip） */}
      {showDirection && (
        <div>
          <Label className="text-xs text-muted-foreground mb-1.5 block">
            {t('方向')}
          </Label>
          <Select
            value={transition.direction ?? 'from-right'}
            onValueChange={handleDirectionChange}
          >
            <SelectTrigger className="h-7 text-xs">
              <SelectValue />
            </SelectTrigger>
            <SelectContent>
              {DIRECTION_OPTIONS.map(({ value, labelKey, icon: Icon }) => (
                <SelectItem key={value} value={value}>
                  <span className="flex items-center gap-1.5">
                    <Icon className="h-3.5 w-3.5" />
                    {t(labelKey)}
                  </span>
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
      )}
    </div>
  );
}

// ============================================================================
// TransitionIndicator
// ============================================================================

export interface TransitionIndicatorProps {
  /** 过渡效果 */
  transition: Transition;
  /** 点击回调 */
  onClick: () => void;
  /** 是否选中 */
  isSelected?: boolean;
}

/**
 * 过渡效果指示条 - 用于场景列表中显示过渡效果
 */
export function TransitionIndicator({
  transition,
  onClick,
  isSelected,
}: TransitionIndicatorProps) {
  const t = useT();
  const Icon = TYPE_ICON_MAP[transition.type];

  // 无过渡时显示简化样式
  if (transition.type === 'none') {
    return (
      <button
        type="button"
        onClick={onClick}
        className={cn(
          'w-full flex items-center justify-center gap-1 py-0.5 text-xs',
          'text-muted-foreground/50 hover:text-muted-foreground hover:bg-muted/30',
          'rounded transition-colors cursor-pointer',
          isSelected && 'bg-muted/50 text-muted-foreground',
        )}
      >
        <ChevronRight className="h-3 w-3" />
      </button>
    );
  }

  return (
    <button
      type="button"
      onClick={onClick}
      className={cn(
        'w-full flex items-center justify-center gap-1.5 py-1 text-xs',
        'rounded transition-colors cursor-pointer',
        'text-muted-foreground hover:bg-accent/50',
        isSelected && 'bg-accent text-accent-foreground',
      )}
    >
      <Icon className="h-3 w-3" />
      <span className="tabular-nums">{transition.durationInFrames}{t('帧')}</span>
    </button>
  );
}
