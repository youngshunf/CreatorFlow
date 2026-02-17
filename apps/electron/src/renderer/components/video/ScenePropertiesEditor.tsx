/**
 * ScenePropertiesEditor - 场景属性编辑器
 *
 * 根据场景的 compositionId 动态渲染对应的属性表单，
 * 支持基础属性（名称、组合类型、时长）和动态 props 编辑。
 *
 * @requirements 9.4
 */

import * as React from 'react';
import { useCallback, useMemo } from 'react';
import { Palette, Plus, Trash2, GripVertical, ChevronDown, ChevronRight } from 'lucide-react';
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
import type { Scene } from '@sprouty-ai/video';
import { COMPOSITION_IDS } from '@sprouty-ai/video';
import { DEFAULT_COLOR_SCHEMES } from '@sprouty-ai/video';
import type { ColorSchemeConfig } from '@sprouty-ai/video';

// ============================================================================
// 类型定义
// ============================================================================

interface PropFieldConfig {
  key: string;
  label: string;
  type: 'text' | 'number' | 'select' | 'colorScheme' | 'itemList' | 'textarea';
  options?: { value: string; label: string }[];
  placeholder?: string;
}

export interface ScenePropertiesEditorProps {
  /** 当前选中的场景 */
  scene: Scene;
  /** 视频帧率 */
  fps: number;
  /** 场景更新回调 */
  onUpdate: (updates: Partial<Scene>) => void;
}

// ============================================================================
// 每个 compositionId 对应的可编辑字段配置
// ============================================================================

const ANIMATION_STYLE_OPTIONS = [
  { value: 'fade', label: '淡入' },
  { value: 'slide', label: '滑入' },
  { value: 'scale', label: '缩放' },
  { value: 'spring', label: '弹性' },
  { value: 'zoom', label: '变焦' },
];

const COMPOSITION_PROPS_CONFIG: Record<string, PropFieldConfig[]> = {
  TitleAnimation: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  Slideshow: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'items', label: '幻灯片', type: 'itemList' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  DataVisualization: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'chartType', label: '图表类型', type: 'select', options: [
      { value: 'bar', label: '柱状图' },
      { value: 'line', label: '折线图' },
      { value: 'pie', label: '饼图' },
      { value: 'donut', label: '环形图' },
    ]},
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  ProductShowcase: [
    { key: 'title', label: '产品名称', type: 'text' },
    { key: 'subtitle', label: '产品标语', type: 'text' },
    { key: 'items', label: '功能列表', type: 'itemList' },
    { key: 'layout', label: '布局', type: 'select', options: [
      { value: 'centered', label: '居中' },
      { value: 'split', label: '分栏' },
      { value: 'features-grid', label: '网格' },
    ]},
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  SocialMediaVertical: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  SocialMediaSquare: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  StepByStepTutorial: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'items', label: '步骤列表', type: 'itemList' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  Explainer: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'items', label: '要点列表', type: 'itemList' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  Tips: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'items', label: '技巧列表', type: 'itemList' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  ProductMarketing: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'items', label: '功能列表', type: 'itemList' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
  PromoAd: [
    { key: 'title', label: '标题', type: 'text' },
    { key: 'subtitle', label: '副标题', type: 'text' },
    { key: 'discount', label: '折扣标签', type: 'text', placeholder: '例如: -50%' },
    { key: 'colors', label: '配色方案', type: 'colorScheme' },
    { key: 'animationStyle', label: '动画风格', type: 'select', options: ANIMATION_STYLE_OPTIONS },
  ],
};

// ============================================================================
// compositionId 中文名称映射
// ============================================================================

const COMPOSITION_LABELS: Record<string, string> = {
  TitleAnimation: '标题动画',
  Slideshow: '幻灯片',
  DataVisualization: '数据可视化',
  ProductShowcase: '产品展示',
  SocialMediaVertical: '竖版社交媒体',
  SocialMediaSquare: '方形社交媒体',
  StepByStepTutorial: '分步教程',
  Explainer: '概念讲解',
  Tips: '技巧清单',
  ProductMarketing: '产品营销',
  PromoAd: '促销广告',
};

// ============================================================================
// 预设配色方案名称
// ============================================================================

const COLOR_SCHEME_LABELS: Record<string, string> = {
  modern: '现代',
  minimal: '极简',
  playful: '活泼',
  corporate: '商务',
  cinematic: '电影',
  vibrant: '鲜艳',
  nature: '自然',
};

// ============================================================================
// ColorSchemeEditor 子组件
// ============================================================================

interface ColorSchemeEditorProps {
  value: ColorSchemeConfig | undefined;
  onChange: (colors: ColorSchemeConfig) => void;
}

function ColorSchemeEditor({ value, onChange }: ColorSchemeEditorProps) {
  const t = useT();
  const colors = value ?? DEFAULT_COLOR_SCHEMES.modern;

  const handleColorChange = useCallback(
    (key: keyof ColorSchemeConfig, color: string) => {
      onChange({ ...colors, [key]: color });
    },
    [colors, onChange]
  );

  const handlePresetClick = useCallback(
    (schemeName: keyof typeof DEFAULT_COLOR_SCHEMES) => {
      onChange({ ...DEFAULT_COLOR_SCHEMES[schemeName] });
    },
    [onChange]
  );

  const colorFields: { key: keyof ColorSchemeConfig; label: string }[] = [
    { key: 'primary', label: '主色' },
    { key: 'secondary', label: '辅色' },
    { key: 'background', label: '背景' },
    { key: 'text', label: '文字' },
  ];

  return (
    <div className="space-y-2">
      {/* 预设配色按钮 */}
      <div className="flex flex-wrap gap-1">
        {(Object.keys(DEFAULT_COLOR_SCHEMES) as Array<keyof typeof DEFAULT_COLOR_SCHEMES>).map(
          (name) => (
            <button
              key={name}
              type="button"
              className="flex items-center gap-1 px-2 py-1 text-xs rounded border border-border hover:bg-muted/50 transition-colors"
              onClick={() => handlePresetClick(name)}
              title={t(COLOR_SCHEME_LABELS[name] ?? name)}
            >
              <span
                className="w-3 h-3 rounded-full border border-border"
                style={{ backgroundColor: DEFAULT_COLOR_SCHEMES[name].primary }}
              />
              <span>{t(COLOR_SCHEME_LABELS[name] ?? name)}</span>
            </button>
          )
        )}
      </div>

      {/* 颜色输入 */}
      <div className="grid grid-cols-2 gap-2">
        {colorFields.map(({ key, label }) => (
          <div key={key} className="flex items-center gap-1.5">
            <input
              type="color"
              value={colors[key]}
              onChange={(e) => handleColorChange(key, e.target.value)}
              className="w-6 h-6 rounded border border-border cursor-pointer p-0"
            />
            <div className="flex-1 min-w-0">
              <Label className="text-[10px] text-muted-foreground">{t(label)}</Label>
              <Input
                value={colors[key]}
                onChange={(e) => handleColorChange(key, e.target.value)}
                className="h-6 text-[10px] font-mono px-1"
              />
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

// ============================================================================
// ItemListEditor 子组件
// ============================================================================

interface ItemData {
  title: string;
  description?: string;
  image?: string;
  icon?: string;
}

interface ItemListEditorProps {
  value: ItemData[] | undefined;
  onChange: (items: ItemData[]) => void;
  label: string;
}

function ItemListEditor({ value, onChange, label }: ItemListEditorProps) {
  const t = useT();
  const items = value ?? [];
  const [expandedIndex, setExpandedIndex] = React.useState<number | null>(null);

  const handleAdd = useCallback(() => {
    onChange([...items, { title: '' }]);
    setExpandedIndex(items.length);
  }, [items, onChange]);

  const handleRemove = useCallback(
    (index: number) => {
      const next = items.filter((_, i) => i !== index);
      onChange(next);
      if (expandedIndex === index) setExpandedIndex(null);
    },
    [items, onChange, expandedIndex]
  );

  const handleItemChange = useCallback(
    (index: number, field: keyof ItemData, fieldValue: string) => {
      const next = items.map((item, i) =>
        i === index ? { ...item, [field]: fieldValue } : item
      );
      onChange(next);
    },
    [items, onChange]
  );

  const handleMoveUp = useCallback(
    (index: number) => {
      if (index === 0) return;
      const next = [...items];
      [next[index - 1], next[index]] = [next[index], next[index - 1]];
      onChange(next);
      setExpandedIndex(index - 1);
    },
    [items, onChange]
  );

  const handleMoveDown = useCallback(
    (index: number) => {
      if (index >= items.length - 1) return;
      const next = [...items];
      [next[index], next[index + 1]] = [next[index + 1], next[index]];
      onChange(next);
      setExpandedIndex(index + 1);
    },
    [items, onChange]
  );

  return (
    <div className="space-y-1.5">
      <div className="flex items-center justify-between">
        <Label className="text-xs text-muted-foreground">{t(label)}</Label>
        <Button variant="ghost" size="sm" className="h-6 px-1.5" onClick={handleAdd}>
          <Plus className="h-3 w-3 mr-1" />
          <span className="text-xs">{t('添加')}</span>
        </Button>
      </div>

      {items.length === 0 && (
        <div className="text-xs text-muted-foreground py-2 text-center">
          {t('暂无项目，点击添加')}
        </div>
      )}

      {items.map((item, index) => (
        <div key={index} className="border border-border rounded p-1.5 space-y-1">
          {/* 项目头部 */}
          <div className="flex items-center gap-1">
            <button
              type="button"
              className="p-0.5 hover:bg-muted/50 rounded"
              onClick={() => setExpandedIndex(expandedIndex === index ? null : index)}
            >
              {expandedIndex === index ? (
                <ChevronDown className="h-3 w-3" />
              ) : (
                <ChevronRight className="h-3 w-3" />
              )}
            </button>
            <span className="text-xs text-muted-foreground w-5">{index + 1}.</span>
            <Input
              value={item.title}
              onChange={(e) => handleItemChange(index, 'title', e.target.value)}
              placeholder={t('标题')}
              className="h-6 text-xs flex-1"
            />
            <div className="flex items-center gap-0.5">
              <button
                type="button"
                className="p-0.5 hover:bg-muted/50 rounded text-muted-foreground disabled:opacity-30"
                onClick={() => handleMoveUp(index)}
                disabled={index === 0}
              >
                <GripVertical className="h-3 w-3 rotate-90" />
              </button>
              <button
                type="button"
                className="p-0.5 hover:bg-destructive/10 rounded text-destructive"
                onClick={() => handleRemove(index)}
              >
                <Trash2 className="h-3 w-3" />
              </button>
            </div>
          </div>

          {/* 展开编辑区 */}
          {expandedIndex === index && (
            <div className="pl-6 space-y-1.5">
              <div>
                <Label className="text-[10px] text-muted-foreground">{t('描述')}</Label>
                <Textarea
                  value={item.description ?? ''}
                  onChange={(e) => handleItemChange(index, 'description', e.target.value)}
                  placeholder={t('可选描述...')}
                  className="h-14 text-xs resize-none mt-0.5"
                />
              </div>
              <div>
                <Label className="text-[10px] text-muted-foreground">{t('图标')}</Label>
                <Input
                  value={item.icon ?? ''}
                  onChange={(e) => handleItemChange(index, 'icon', e.target.value)}
                  placeholder={t('图标名称')}
                  className="h-6 text-xs mt-0.5"
                />
              </div>
              <div>
                <Label className="text-[10px] text-muted-foreground">{t('图片路径')}</Label>
                <Input
                  value={item.image ?? ''}
                  onChange={(e) => handleItemChange(index, 'image', e.target.value)}
                  placeholder={t('images/example.jpg')}
                  className="h-6 text-xs mt-0.5"
                />
              </div>
            </div>
          )}
        </div>
      ))}
    </div>
  );
}

// ============================================================================
// ScenePropertiesEditor 主组件
// ============================================================================

export function ScenePropertiesEditor({ scene, fps, onUpdate }: ScenePropertiesEditorProps) {
  const t = useT();

  // 当前 compositionId 对应的字段配置
  const fieldConfigs = useMemo(
    () => COMPOSITION_PROPS_CONFIG[scene.compositionId] ?? [],
    [scene.compositionId]
  );

  // 更新场景 props 中的某个字段
  const handlePropChange = useCallback(
    (key: string, value: unknown) => {
      onUpdate({
        props: { ...scene.props, [key]: value },
      });
    },
    [scene.props, onUpdate]
  );

  // 切换 compositionId 时重置 props
  const handleCompositionChange = useCallback(
    (compositionId: string) => {
      onUpdate({
        compositionId,
        props: {},
      });
    },
    [onUpdate]
  );

  // 更新场景名称
  const handleNameChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      onUpdate({ name: e.target.value });
    },
    [onUpdate]
  );

  // 更新时长（帧数）
  const handleDurationChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const frames = parseInt(e.target.value) || 1;
      onUpdate({ durationInFrames: Math.max(1, frames) });
    },
    [onUpdate]
  );

  // 渲染单个 prop 字段
  const renderField = useCallback(
    (field: PropFieldConfig) => {
      const propValue = scene.props[field.key];

      switch (field.type) {
        case 'text':
          return (
            <div key={field.key}>
              <Label className="text-xs text-muted-foreground">{t(field.label)}</Label>
              <Input
                value={(propValue as string) ?? ''}
                onChange={(e) => handlePropChange(field.key, e.target.value)}
                placeholder={field.placeholder ? t(field.placeholder) : undefined}
                className="h-7 text-xs mt-0.5"
              />
            </div>
          );

        case 'textarea':
          return (
            <div key={field.key}>
              <Label className="text-xs text-muted-foreground">{t(field.label)}</Label>
              <Textarea
                value={(propValue as string) ?? ''}
                onChange={(e) => handlePropChange(field.key, e.target.value)}
                placeholder={field.placeholder ? t(field.placeholder) : undefined}
                className="h-16 text-xs resize-none mt-0.5"
              />
            </div>
          );

        case 'number':
          return (
            <div key={field.key}>
              <Label className="text-xs text-muted-foreground">{t(field.label)}</Label>
              <Input
                type="number"
                value={(propValue as number) ?? 0}
                onChange={(e) => handlePropChange(field.key, parseFloat(e.target.value) || 0)}
                className="h-7 text-xs mt-0.5"
              />
            </div>
          );

        case 'select':
          return (
            <div key={field.key}>
              <Label className="text-xs text-muted-foreground">{t(field.label)}</Label>
              <Select
                value={(propValue as string) ?? ''}
                onValueChange={(v) => handlePropChange(field.key, v)}
              >
                <SelectTrigger className="h-7 text-xs mt-0.5">
                  <SelectValue placeholder={t('选择...')} />
                </SelectTrigger>
                <SelectContent>
                  {field.options?.map((opt) => (
                    <SelectItem key={opt.value} value={opt.value}>
                      {t(opt.label)}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>
          );

        case 'colorScheme':
          return (
            <div key={field.key}>
              <Label className="text-xs text-muted-foreground flex items-center gap-1 mb-1">
                <Palette className="h-3 w-3" />
                {t(field.label)}
              </Label>
              <ColorSchemeEditor
                value={propValue as ColorSchemeConfig | undefined}
                onChange={(colors) => handlePropChange(field.key, colors)}
              />
            </div>
          );

        case 'itemList':
          return (
            <ItemListEditor
              key={field.key}
              value={propValue as ItemData[] | undefined}
              onChange={(items) => handlePropChange(field.key, items)}
              label={field.label}
            />
          );

        default:
          return null;
      }
    },
    [scene.props, handlePropChange, t]
  );

  return (
    <div className="space-y-3">
      {/* 基础属性区 */}
      <div className="space-y-2">
        <div className="text-xs font-medium text-foreground">{t('场景属性')}</div>

        {/* 场景名称 */}
        <div>
          <Label className="text-xs text-muted-foreground">{t('场景名称')}</Label>
          <Input
            value={scene.name}
            onChange={handleNameChange}
            className="h-7 text-xs mt-0.5"
          />
        </div>

        {/* compositionId 选择器 */}
        <div>
          <Label className="text-xs text-muted-foreground">{t('组合类型')}</Label>
          <Select value={scene.compositionId} onValueChange={handleCompositionChange}>
            <SelectTrigger className="h-7 text-xs mt-0.5">
              <SelectValue />
            </SelectTrigger>
            <SelectContent>
              {COMPOSITION_IDS.map((id) => (
                <SelectItem key={id} value={id}>
                  {t(COMPOSITION_LABELS[id] ?? id)}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>

        {/* 时长编辑 */}
        <div className="flex items-end gap-2">
          <div className="flex-1">
            <Label className="text-xs text-muted-foreground">{t('时长(帧)')}</Label>
            <Input
              type="number"
              value={scene.durationInFrames}
              onChange={handleDurationChange}
              min={1}
              className="h-7 text-xs mt-0.5"
            />
          </div>
          <div className="text-xs text-muted-foreground pb-1.5">
            {(scene.durationInFrames / fps).toFixed(2)}s
          </div>
        </div>
      </div>

      {/* 动态 props 表单 */}
      {fieldConfigs.length > 0 && (
        <>
          <Separator />
          <div className="space-y-2">
            <div className="text-xs font-medium text-foreground">
              {t('组合参数')} - {t(COMPOSITION_LABELS[scene.compositionId] ?? scene.compositionId)}
            </div>
            {fieldConfigs.map(renderField)}
          </div>
        </>
      )}
    </div>
  );
}
