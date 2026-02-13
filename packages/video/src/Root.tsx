/**
 * Remotion Root Component
 *
 * This is the entry point for Remotion. It registers all available
 * compositions that can be rendered or previewed.
 *
 * @requirements 1.3
 */
import { Composition, Folder } from 'remotion';

// Import compositions
import { TitleAnimation } from './compositions/TitleAnimation';
import { Slideshow } from './compositions/Slideshow';
import { DataVisualization } from './compositions/DataVisualization';
import { ProductShowcase } from './compositions/ProductShowcase';

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
  primary: '#6366f1',
  secondary: '#8b5cf6',
  background: '#1a1a2e',
  text: '#ffffff',
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
            title: 'Welcome',
            subtitle: 'Your subtitle here',
            colors: DEFAULT_COLORS,
            animationStyle: 'spring',
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
            title: 'Photo Gallery',
            items: [
              { title: 'Slide 1', description: 'First slide description' },
              { title: 'Slide 2', description: 'Second slide description' },
              { title: 'Slide 3', description: 'Third slide description' },
            ],
            colors: DEFAULT_COLORS,
            slideDuration: 90,
            transitionDuration: 30,
            transitionStyle: 'crossfade',
            showTitles: true,
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const slideDuration = props.slideDuration || 90;
            const transitionDuration = props.transitionDuration || 30;
            const itemCount = props.items?.length || 3;
            const totalDuration = itemCount * (slideDuration + transitionDuration);

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
            title: 'Data Overview',
            subtitle: 'Monthly Statistics',
            data: [
              { label: 'Jan', value: 65 },
              { label: 'Feb', value: 85 },
              { label: 'Mar', value: 45 },
              { label: 'Apr', value: 95 },
              { label: 'May', value: 75 },
              { label: 'Jun', value: 55 },
            ],
            colors: DEFAULT_COLORS,
            chartType: 'bar',
            showLabels: true,
            showValues: true,
            staggerDelay: 5,
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2); // 2秒介绍
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 6;
            const chartDuration = Math.floor(fps * 3) + dataCount * staggerDelay;
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
            title: 'Trend Analysis',
            subtitle: 'Growth Over Time',
            data: [
              { label: 'Q1', value: 30 },
              { label: 'Q2', value: 55 },
              { label: 'Q3', value: 70 },
              { label: 'Q4', value: 95 },
            ],
            colors: DEFAULT_COLORS,
            chartType: 'line',
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration = Math.floor(fps * 3) + dataCount * staggerDelay;
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
            title: 'Market Share',
            subtitle: 'Distribution by Category',
            data: [
              { label: 'Product A', value: 35 },
              { label: 'Product B', value: 25 },
              { label: 'Product C', value: 20 },
              { label: 'Product D', value: 20 },
            ],
            colors: DEFAULT_COLORS,
            chartType: 'pie',
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration = Math.floor(fps * 3) + dataCount * staggerDelay;
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
            title: 'Budget Allocation',
            subtitle: 'Annual Spending',
            data: [
              { label: 'Marketing', value: 40 },
              { label: 'Development', value: 30 },
              { label: 'Operations', value: 20 },
              { label: 'Other', value: 10 },
            ],
            colors: DEFAULT_COLORS,
            chartType: 'donut',
            showLabels: true,
            showValues: true,
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const staggerDelay = props.staggerDelay || 5;
            const dataCount = props.data?.length || 4;
            const chartDuration = Math.floor(fps * 3) + dataCount * staggerDelay;
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
            title: 'Amazing Product',
            subtitle: 'The future of productivity',
            items: [
              {
                title: 'Feature One',
                description: 'Amazing feature that helps you achieve more',
                icon: '🚀',
              },
              {
                title: 'Feature Two',
                description: 'Powerful capability for better results',
                icon: '⚡',
              },
              {
                title: 'Feature Three',
                description: 'Smart solution for modern challenges',
                icon: '🎯',
              },
              {
                title: 'Feature Four',
                description: 'Innovative approach to common problems',
                icon: '💡',
              },
            ],
            colors: DEFAULT_COLORS,
            layout: 'centered',
            animationStyle: 'spring',
            cta: {
              text: 'Get Started',
            },
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2); // 2秒介绍
            const featureCount = props.items?.length || 4;
            const featureDuration = featureCount * Math.floor(fps * 2); // 每个功能2秒
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0; // 3秒CTA
            const outroDuration = Math.floor(fps * 1); // 1秒结尾

            return {
              durationInFrames: introDuration + featureDuration + ctaDuration + outroDuration,
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
            title: 'Product Name',
            subtitle: 'Your tagline here',
            items: [
              { title: 'Feature 1', description: 'Description here', icon: '✨' },
              { title: 'Feature 2', description: 'Description here', icon: '🔥' },
              { title: 'Feature 3', description: 'Description here', icon: '💎' },
            ],
            colors: DEFAULT_COLORS,
            layout: 'split',
            animationStyle: 'slide',
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const featureCount = props.items?.length || 3;
            const featureDuration = featureCount * Math.floor(fps * 2);
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames: introDuration + featureDuration + ctaDuration + outroDuration,
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
            title: 'Key Features',
            subtitle: 'Everything you need',
            items: [
              { title: 'Speed', description: 'Lightning fast performance', icon: '⚡' },
              { title: 'Security', description: 'Enterprise-grade protection', icon: '🔒' },
              { title: 'Scale', description: 'Grows with your needs', icon: '📈' },
              { title: 'Support', description: '24/7 expert assistance', icon: '🤝' },
            ],
            colors: DEFAULT_COLORS,
            layout: 'features-grid',
            animationStyle: 'spring',
          }}
          calculateMetadata={({ props }) => {
            const fps = 30;
            const introDuration = Math.floor(fps * 2);
            const featureCount = props.items?.length || 4;
            const featureDuration = featureCount * Math.floor(fps * 2);
            const ctaDuration = props.cta ? Math.floor(fps * 3) : 0;
            const outroDuration = Math.floor(fps * 1);

            return {
              durationInFrames: introDuration + featureDuration + ctaDuration + outroDuration,
            };
          }}
        />
      </Folder>
    </>
  );
};

export default RemotionRoot;
