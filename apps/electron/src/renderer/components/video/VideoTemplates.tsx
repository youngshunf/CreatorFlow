/**
 * VideoTemplates - Template selector component
 *
 * Features:
 * - Display available video templates
 * - Template preview and selection
 * - Category filtering
 *
 * @requirements 9.6
 */

import * as React from "react";
import { useState, useMemo, Component } from "react";
import { Player } from "@remotion/player";
import {
  Film,
  Smartphone,
  Megaphone,
  GraduationCap,
  Check,
  AlertTriangle,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@/components/ui/tabs";
import { useT } from "@/context/LocaleContext";
import {
  ALL_TEMPLATES,
  TEMPLATES_BY_CATEGORY,
  TitleAnimation,
  Slideshow,
  DataVisualization,
  ProductShowcase,
  SocialMediaVertical,
  SocialMediaSquare,
  StepByStepTutorial,
  Explainer,
  Tips,
  ProductMarketing,
  PromoAd,
  type VideoTemplate,
  type TemplateCategory,
} from "@sprouty-ai/video";

/**
 * compositionId → React 组件的映射（全部 11 个内置组合）
 */
const COMPOSITION_MAP: Record<string, React.FC<any>> = {
  TitleAnimation,
  Slideshow,
  DataVisualization,
  ProductShowcase,
  SocialMediaVertical,
  SocialMediaSquare,
  StepByStepTutorial,
  Explainer,
  Tips,
  ProductMarketing,
  PromoAd,
};

export interface VideoTemplatesProps {
  /** Callback when template is selected for preview */
  onSelect: (template: VideoTemplate) => void;
  /** Callback when user wants to create a project from template */
  onCreate?: (template: VideoTemplate) => void;
  /** Optional class name */
  className?: string;
}

/**
 * Category icon mapping
 */
const CATEGORY_ICONS: Record<TemplateCategory, React.ElementType> = {
  "social-media": Smartphone,
  marketing: Megaphone,
  tutorial: GraduationCap,
};

/**
 * Category label mapping
 */
const CATEGORY_LABELS: Record<TemplateCategory, string> = {
  "social-media": "社交媒体",
  marketing: "营销推广",
  tutorial: "教程",
};

/**
 * Error boundary for template thumbnail Player
 */
class TemplateThumbnailErrorBoundary extends Component<
  { children: React.ReactNode; compositionId: string },
  { hasError: boolean; error: Error | null }
> {
  state = { hasError: false, error: null as Error | null };

  static getDerivedStateFromError(error: Error) {
    return { hasError: true, error };
  }

  render() {
    if (this.state.hasError) {
      return (
        <div
          style={{
            height: "100%",
            display: "flex",
            flexDirection: "column",
            alignItems: "center",
            justifyContent: "center",
            backgroundColor: "#1a1a2e",
            color: "#94a3b8",
            fontSize: 12,
            padding: 16,
            textAlign: "center",
          }}
        >
          <AlertTriangle
            style={{ width: 24, height: 24, marginBottom: 8, opacity: 0.5 }}
          />
          <div>{this.props.compositionId}</div>
          <div style={{ fontSize: 10, opacity: 0.6, marginTop: 4 }}>
            {this.state.error?.message?.slice(0, 80) || "预览加载失败"}
          </div>
        </div>
      );
    }
    return this.props.children;
  }
}

/**
 * Template card component
 */
function TemplateCard({
  template,
  onSelect,
  onCreate,
}: {
  template: VideoTemplate;
  onSelect: () => void;
  onCreate?: () => void;
}) {
  const t = useT();
  const CategoryIcon = CATEGORY_ICONS[template.category];

  // Resolve the composition component for thumbnail
  const CompositionComponent = useMemo(() => {
    return COMPOSITION_MAP[template.compositionId] || TitleAnimation;
  }, [template.compositionId]);

  return (
    <div
      className={cn(
        "group relative rounded-lg border bg-card overflow-hidden cursor-pointer",
        "hover:border-primary/50 hover:shadow-md transition-all",
      )}
      onClick={onSelect}
    >
      {/* Thumbnail preview */}
      <div
        className="relative bg-black overflow-hidden flex items-center justify-center"
        style={{
          height: "160px",
        }}
      >
        <TemplateThumbnailErrorBoundary compositionId={template.compositionId}>
          <Player
            component={CompositionComponent}
            inputProps={template.defaultProps || {}}
            durationInFrames={template.defaultConfig.durationInFrames}
            fps={template.defaultConfig.fps}
            compositionWidth={template.defaultConfig.width}
            compositionHeight={template.defaultConfig.height}
            style={{
              height: "100%",
              aspectRatio: `${template.defaultConfig.width} / ${template.defaultConfig.height}`,
            }}
            acknowledgeRemotionLicense
            autoPlay={false}
            controls={false}
          />
        </TemplateThumbnailErrorBoundary>
        {/* Aspect ratio badge */}
        {template.aspectRatio && (
          <Badge
            variant="secondary"
            className="absolute top-2 right-2 text-[10px] px-1.5 py-0"
          >
            {template.aspectRatio}
          </Badge>
        )}
      </div>

      {/* Info */}
      <div className="p-3">
        <div className="flex items-start justify-between gap-2">
          <div className="flex-1 min-w-0">
            <h4 className="font-medium text-sm truncate">{t(template.name)}</h4>
            <p className="text-xs text-muted-foreground mt-0.5 line-clamp-2">
              {t(template.description)}
            </p>
          </div>
          <CategoryIcon className="h-4 w-4 text-muted-foreground shrink-0" />
        </div>

        {/* Tags */}
        {template.tags && template.tags.length > 0 && (
          <div className="flex flex-wrap gap-1 mt-2">
            {template.tags.slice(0, 3).map((tag) => (
              <Badge
                key={tag}
                variant="outline"
                className="text-[10px] px-1.5 py-0"
              >
                {tag}
              </Badge>
            ))}
          </div>
        )}

        {/* Config info */}
        <div className="flex items-center gap-2 mt-2 text-[10px] text-muted-foreground">
          <span>
            {template.defaultConfig.width}×{template.defaultConfig.height}
          </span>
          <span>•</span>
          <span>{template.defaultConfig.fps}fps</span>
          <span>•</span>
          <span>
            {(
              template.defaultConfig.durationInFrames /
              template.defaultConfig.fps
            ).toFixed(1)}
            s
          </span>
        </div>
      </div>

      {/* Hover overlay */}
      <div
        className={cn(
          "absolute inset-0 bg-primary/10 opacity-0 group-hover:opacity-100",
          "flex items-center justify-center transition-opacity pointer-events-none",
        )}
      >
        <Button
          size="sm"
          className="shadow-lg pointer-events-auto"
          onClick={(e) => {
            e.stopPropagation();
            onCreate?.();
          }}
        >
          <Check className="h-4 w-4 mr-1" />
          {t("使用模板")}
        </Button>
      </div>
    </div>
  );
}

/**
 * VideoTemplates component
 */
export function VideoTemplates({
  onSelect,
  onCreate,
  className,
}: VideoTemplatesProps) {
  const t = useT();
  const [category, setCategory] = useState<"all" | TemplateCategory>("all");

  // Filter templates by category
  const filteredTemplates = useMemo(() => {
    if (category === "all") {
      return ALL_TEMPLATES;
    }
    return TEMPLATES_BY_CATEGORY[category] || [];
  }, [category]);

  return (
    <div className={cn("p-2", className)}>
      {/* Category tabs */}
      <Tabs
        value={category}
        onValueChange={(v) => setCategory(v as "all" | TemplateCategory)}
        className="mb-3"
      >
        <TabsList className="h-8 w-full grid grid-cols-4">
          <TabsTrigger value="all" className="text-xs px-1">
            {t("全部")}
          </TabsTrigger>
          <TabsTrigger value="social-media" className="text-xs px-1">
            <Smartphone className="h-3 w-3 mr-1" />
            {t("社交")}
          </TabsTrigger>
          <TabsTrigger value="marketing" className="text-xs px-1">
            <Megaphone className="h-3 w-3 mr-1" />
            {t("营销")}
          </TabsTrigger>
          <TabsTrigger value="tutorial" className="text-xs px-1">
            <GraduationCap className="h-3 w-3 mr-1" />
            {t("教程")}
          </TabsTrigger>
        </TabsList>
      </Tabs>

      {/* Template grid */}
      <div className="grid grid-cols-1 gap-3">
        {filteredTemplates.map((template) => (
          <TemplateCard
            key={template.id}
            template={template}
            onSelect={() => onSelect(template)}
            onCreate={() => onCreate?.(template)}
          />
        ))}
      </div>

      {filteredTemplates.length === 0 && (
        <div className="text-center text-muted-foreground py-8">
          <Film className="h-12 w-12 mx-auto mb-2 opacity-30" />
          <p className="text-sm">{t("暂无模板")}</p>
        </div>
      )}
    </div>
  );
}
