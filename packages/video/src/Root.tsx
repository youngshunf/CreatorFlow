/**
 * Remotion Root Component
 *
 * This is the entry point for Remotion. It registers all available
 * compositions that can be rendered or previewed.
 *
 * @requirements 1.3
 */
import { Composition, Folder } from "remotion";

// 场景编排器（核心组合）
import {
  SceneComposer,
  SceneComposerSchema,
  calculateSceneComposerMetadata,
} from "./compositions/SceneComposer";

// Import compositions
import { TitleAnimation } from "./compositions/TitleAnimation";
import { Slideshow } from "./compositions/Slideshow";
import { DataVisualization } from "./compositions/DataVisualization";
import { ProductShowcase } from "./compositions/ProductShowcase";
import { SocialMediaVertical } from "./compositions/SocialMediaVertical";
import { SocialMediaSquare } from "./compositions/SocialMediaSquare";
import { StepByStepTutorial } from "./compositions/StepByStepTutorial";
import { Explainer } from "./compositions/Explainer";
import { Tips } from "./compositions/Tips";
import { ProductMarketing } from "./compositions/ProductMarketing";
import { PromoAd } from "./compositions/PromoAd";

/**
 * Default video configuration
 */
const DEFAULT_CONFIG = {
  width: 1920,
  height: 1080,
  fps: 30,
  durationInFrames: 300, // 10 seconds at 30fps
};

/**
 * Default colors for compositions
 */
const DEFAULT_COLORS = {
  primary: "#6366f1",
  secondary: "#8b5cf6",
  background: "#1a1a2e",
  text: "#ffffff",
};

/**
 * RemotionRoot - Registers all available compositions
 *
 * This component is used by Remotion Studio and the renderer
 * to discover and render video compositions.
 */
export const RemotionRoot: React.FC = () => {
  return (
    <>
      {/* 核心：场景编排器组合 */}
      <Composition
        id="SceneComposer"
        component={SceneComposer}
        schema={SceneComposerSchema}
        calculateMetadata={calculateSceneComposerMetadata}
        durationInFrames={300}
        fps={DEFAULT_CONFIG.fps}
        width={DEFAULT_CONFIG.width}
        height={DEFAULT_CONFIG.height}
        defaultProps={{
          scenes: [
            {
              id: "demo-title",
              name: "标题",
              compositionId: "TitleAnimation",
              durationInFrames: 90,
              props: {
                title: "SceneComposer 演示",
                subtitle: "多场景编排",
                colors: DEFAULT_COLORS,
                animationStyle: "spring",
              },
            },
          ],
          transitions: [],
        }}
      />

      {/* 基础动画 */}
      <Folder name="基础动画">
        {/* TitleAnimation composition */}
        <Composition
          id="TitleAnimation"
          component={TitleAnimation}
          durationInFrames={150} // 5 seconds at 30fps
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Welcome",
            subtitle: "Your subtitle here",
            colors: DEFAULT_COLORS,
            animationStyle: "spring",
            titleFontSize: 72,
            subtitleFontSize: 36,
          }}
        />

        {/* Slideshow composition */}
        <Composition
          id="Slideshow"
          component={Slideshow}
          durationInFrames={DEFAULT_CONFIG.durationInFrames}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Photo Gallery",
            items: [
              { title: "Slide 1", description: "First slide description" },
              { title: "Slide 2", description: "Second slide description" },
              { title: "Slide 3", description: "Third slide description" },
            ],
            colors: DEFAULT_COLORS,
            slideDuration: 90,
            transitionDuration: 30,
            transitionStyle: "crossfade",
            showTitles: true,
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const slideDuration = props.slideDuration || 90;
            const transitionDuration = props.transitionDuration || 30;
            const itemCount = props.items?.length || 3;
            const totalDuration =
              itemCount * (slideDuration + transitionDuration);

            return {
              durationInFrames: totalDuration,
            };
          }}
        />
      </Folder>

      {/* 数据可视化 */}
      <Folder name="数据可视化">
        {/* DataVisualization composition */}
        <Composition
          id="DataVisualization"
          component={DataVisualization}
          durationInFrames={DEFAULT_CONFIG.durationInFrames}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Data Overview",
            subtitle: "Monthly Statistics",
            data: [
              { label: "Jan", value: 65 },
              { label: "Feb", value: 85 },
              { label: "Mar", value: 45 },
              { label: "Apr", value: 95 },
              { label: "May", value: 75 },
              { label: "Jun", value: 55 },
            ],
            colors: DEFAULT_COLORS,
            chartType: "bar",
            showLabels: true,
            showValues: true,
            staggerDelay: 5,
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2); // 2秒介绍
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 6;
            const chartDuration =
              Math.floor(fps * 3) + dataCount * staggerDelay;
            const outroDuration = Math.floor(fps * 1); // 1秒结尾

            return {
              durationInFrames: introDuration + chartDuration + outroDuration,
            };
          }}
        />

        {/* DataVisualization - Line Chart variant */}
        <Composition
          id="DataVisualization-Line"
          component={DataVisualization}
          durationInFrames={DEFAULT_CONFIG.durationInFrames}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Trend Analysis",
            subtitle: "Growth Over Time",
            data: [
              { label: "Q1", value: 30 },
              { label: "Q2", value: 55 },
              { label: "Q3", value: 70 },
              { label: "Q4", value: 95 },
            ],
            colors: DEFAULT_COLORS,
            chartType: "line",
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration =
              Math.floor(fps * 3) + dataCount * staggerDelay;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames: introDuration + chartDuration + outroDuration,
            };
          }}
        />

        {/* DataVisualization - Pie Chart variant */}
        <Composition
          id="DataVisualization-Pie"
          component={DataVisualization}
          durationInFrames={DEFAULT_CONFIG.durationInFrames}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Market Share",
            subtitle: "Distribution by Category",
            data: [
              { label: "Product A", value: 35 },
              { label: "Product B", value: 25 },
              { label: "Product C", value: 20 },
              { label: "Product D", value: 20 },
            ],
            colors: DEFAULT_COLORS,
            chartType: "pie",
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration =
              Math.floor(fps * 3) + dataCount * staggerDelay;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames: introDuration + chartDuration + outroDuration,
            };
          }}
        />

        {/* DataVisualization - Donut Chart variant */}
        <Composition
          id="DataVisualization-Donut"
          component={DataVisualization}
          durationInFrames={DEFAULT_CONFIG.durationInFrames}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Budget Allocation",
            subtitle: "Annual Spending",
            data: [
              { label: "Marketing", value: 40 },
              { label: "Development", value: 30 },
              { label: "Operations", value: 20 },
              { label: "Other", value: 10 },
            ],
            colors: DEFAULT_COLORS,
            chartType: "donut",
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration =
              Math.floor(fps * 3) + dataCount * staggerDelay;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames: introDuration + chartDuration + outroDuration,
            };
          }}
        />
      </Folder>

      {/* 产品展示 */}
      <Folder name="产品展示">
        {/* ProductShowcase composition - Centered layout */}
        <Composition
          id="ProductShowcase"
          component={ProductShowcase}
          durationInFrames={450} // 15 seconds at 30fps
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Amazing Product",
            subtitle: "The future of productivity",
            items: [
              {
                title: "Feature One",
                description: "Amazing feature that helps you achieve more",
                icon: "🚀",
              },
              {
                title: "Feature Two",
                description: "Powerful capability for better results",
                icon: "⚡",
              },
              {
                title: "Feature Three",
                description: "Smart solution for modern challenges",
                icon: "🎯",
              },
              {
                title: "Feature Four",
                description: "Innovative approach to common problems",
                icon: "💡",
              },
            ],
            colors: DEFAULT_COLORS,
            layout: "centered",
            animationStyle: "spring",
            cta: {
              text: "Get Started",
            },
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2); // 2秒介绍
            const featureCount = props.items?.length || 4;
            const featureDuration = featureCount * Math.floor(fps * 2); // 每个功能2秒
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0; // 3秒CTA
            const outroDuration = Math.floor(fps * 1); // 1秒结尾

            return {
              durationInFrames:
                introDuration + featureDuration + ctaDuration + outroDuration,
            };
          }}
        />

        {/* ProductShowcase - Split layout variant */}
        <Composition
          id="ProductShowcase-Split"
          component={ProductShowcase}
          durationInFrames={450}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Product Name",
            subtitle: "Your tagline here",
            items: [
              {
                title: "Feature 1",
                description: "Description here",
                icon: "✨",
              },
              {
                title: "Feature 2",
                description: "Description here",
                icon: "🔥",
              },
              {
                title: "Feature 3",
                description: "Description here",
                icon: "💎",
              },
            ],
            colors: DEFAULT_COLORS,
            layout: "split",
            animationStyle: "slide",
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const featureCount = props.items?.length || 3;
            const featureDuration = featureCount * Math.floor(fps * 2);
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames:
                introDuration + featureDuration + ctaDuration + outroDuration,
            };
          }}
        />

        {/* ProductShowcase - Features Grid layout variant */}
        <Composition
          id="ProductShowcase-Grid"
          component={ProductShowcase}
          durationInFrames={450}
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "Key Features",
            subtitle: "Everything you need",
            items: [
              {
                title: "Speed",
                description: "Lightning fast performance",
                icon: "⚡",
              },
              {
                title: "Security",
                description: "Enterprise-grade protection",
                icon: "🔒",
              },
              {
                title: "Scale",
                description: "Grows with your needs",
                icon: "📈",
              },
              {
                title: "Support",
                description: "24/7 expert assistance",
                icon: "🤝",
              },
            ],
            colors: DEFAULT_COLORS,
            layout: "features-grid",
            animationStyle: "spring",
          }}
          calculateMetadata={({ props }: { props: any }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const featureCount = props.items?.length || 4;
            const featureDuration = featureCount * Math.floor(fps * 2);
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames:
                introDuration + featureDuration + ctaDuration + outroDuration,
            };
          }}
        />
      </Folder>

      {/* 社交媒体 */}
      <Folder name="社交媒体">
        {/* SocialMediaVertical - 9:16 竖版 */}
        <Composition
          id="SocialMediaVertical"
          component={SocialMediaVertical}
          durationInFrames={150} // 5 seconds
          fps={DEFAULT_CONFIG.fps}
          width={1080}
          height={1920}
          defaultProps={{
            title: "你的精彩标题",
            subtitle: "上滑了解更多",
            colors: {
              primary: "#ec4899",
              secondary: "#8b5cf6",
              background: "#0f0f23",
              text: "#ffffff",
            },
            cta: {
              text: "了解更多",
            },
          }}
        />

        {/* SocialMediaSquare - 1:1 方形 */}
        <Composition
          id="SocialMediaSquare"
          component={SocialMediaSquare}
          durationInFrames={150} // 5 seconds
          fps={DEFAULT_CONFIG.fps}
          width={1080}
          height={1080}
          defaultProps={{
            title: "你的精彩标题",
            subtitle: "双击点赞",
            colors: {
              primary: "#ec4899",
              secondary: "#8b5cf6",
              background: "#0f0f23",
              text: "#ffffff",
            },
            cta: {
              text: "立即购买",
            },
          }}
        />
      </Folder>

      {/* 教程 */}
      <Folder name="教程">
        {/* StepByStepTutorial - 分步教程 */}
        <Composition
          id="StepByStepTutorial"
          component={StepByStepTutorial}
          durationInFrames={450} // 15 seconds
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "操作指南",
            subtitle: "跟着这些简单步骤操作",
            colors: {
              primary: "#10b981",
              secondary: "#059669",
              background: "#0f172a",
              text: "#f8fafc",
            },
            items: [
              {
                title: "准备工作环境",
                description: "收集所有必要材料并搭建环境",
                icon: "1️⃣",
              },
              {
                title: "按步骤执行",
                description: "仔细且有条理地执行每个步骤",
                icon: "2️⃣",
              },
              {
                title: "检查与优化",
                description: "检查成果并进行必要的调整",
                icon: "3️⃣",
              },
              {
                title: "完成与分享",
                description: "完善项目并分享你的成果",
                icon: "4️⃣",
              },
            ],
          }}
        />

        {/* Explainer - 概念讲解 */}
        <Composition
          id="Explainer"
          component={Explainer}
          durationInFrames={300} // 10 seconds
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "概念解析",
            subtitle: "简单易懂的讲解",
            colors: {
              primary: "#6366f1",
              secondary: "#8b5cf6",
              background: "#ffffff",
              text: "#1f2937",
            },
            items: [
              {
                title: "核心要点一",
                description: "你需要理解的基础概念",
                icon: "💡",
              },
              {
                title: "核心要点二",
                description: "在基础上深入了解更多细节",
                icon: "📊",
              },
              {
                title: "核心要点三",
                description: "实际应用与关键收获",
                icon: "🎯",
              },
            ],
          }}
        />

        {/* Tips - 技巧清单 */}
        <Composition
          id="Tips"
          component={Tips}
          durationInFrames={240} // 8 seconds
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "5 个实用技巧",
            subtitle: "今天就提升你的效率",
            colors: {
              primary: "#f59e0b",
              secondary: "#d97706",
              background: "#1f2937",
              text: "#ffffff",
            },
            items: [
              { title: "早起行动", description: "带着目标开始新的一天" },
              { title: "保持专注", description: "排除干扰，集中精力" },
              { title: "适当休息", description: "休息是为了保持充沛精力" },
              { title: "回顾进展", description: "追踪你的成就和进步" },
              { title: "持续学习", description: "永远不要停止成长" },
            ],
          }}
        />
      </Folder>

      {/* 营销 */}
      <Folder name="营销">
        {/* ProductMarketing - 产品营销 */}
        <Composition
          id="ProductMarketing"
          component={ProductMarketing}
          durationInFrames={300} // 10 seconds
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "全新产品发布",
            subtitle: "创新引领未来",
            colors: {
              primary: "#0ea5e9",
              secondary: "#0284c7",
              background: "#0f172a",
              text: "#f8fafc",
            },
            items: [
              {
                title: "极速性能",
                description: "满足你所有需求的极致性能",
                icon: "⚡",
              },
              {
                title: "智能设计",
                description: "直觉式交互，开箱即用",
                icon: "🎯",
              },
              {
                title: "安全可靠",
                description: "企业级安全防护，内置保障",
                icon: "🔒",
              },
            ],
            cta: {
              text: "立即体验",
            },
          }}
        />

        {/* PromoAd - 促销广告 */}
        <Composition
          id="PromoAd"
          component={PromoAd}
          durationInFrames={150} // 5 seconds
          fps={DEFAULT_CONFIG.fps}
          width={DEFAULT_CONFIG.width}
          height={DEFAULT_CONFIG.height}
          defaultProps={{
            title: "超级大促",
            subtitle: "限时优惠，不容错过！",
            colors: {
              primary: "#f43f5e",
              secondary: "#fb923c",
              background: "#1a1a2e",
              text: "#ffffff",
            },
            cta: {
              text: "立即抢购",
            },
            discount: "50% OFF",
          }}
        />
      </Folder>
    </>
  );
};

export default RemotionRoot;
