/**
 * Slideshow Composition
 *
 * A video composition that displays images in a slideshow format
 * with smooth transition effects between slides.
 *
 * @requirements 3.2, 3.5, 3.6
 */
import React from 'react';
import {
  AbsoluteFill,
  useCurrentFrame,
  useVideoConfig,
  interpolate,
  spring,
  Img,
} from 'remotion';
import type { TemplateProps } from '../templates/types';

/**
 * Props for Slideshow composition
 */
export interface SlideshowProps extends TemplateProps {
  /** Array of image URLs or items with images */
  items?: Array<{
    title: string;
    description?: string;
    image?: string;
  }>;
  /** Duration of each slide in frames */
  slideDuration?: number;
  /** Duration of transition between slides in frames */
  transitionDuration?: number;
  /** Transition style */
  transitionStyle?: 'fade' | 'slide' | 'zoom' | 'crossfade';
  /** Show slide titles */
  showTitles?: boolean;
}

/**
 * Default colors for the composition
 */
const DEFAULT_COLORS = {
  primary: '#6366f1',
  secondary: '#8b5cf6',
  background: '#1a1a2e',
  text: '#ffffff',
};

/**
 * Default items for demo purposes
 */
const DEFAULT_ITEMS = [
  {
    title: 'Slide 1',
    description: 'First slide description',
    image: '',
  },
  {
    title: 'Slide 2',
    description: 'Second slide description',
    image: '',
  },
  {
    title: 'Slide 3',
    description: 'Third slide description',
    image: '',
  },
];

/**
 * Slideshow - Displays images in a slideshow format
 *
 * This composition creates a professional slideshow with customizable
 * transitions and optional title overlays.
 */
export const Slideshow: React.FC<SlideshowProps> = ({
  title,
  items = DEFAULT_ITEMS,
  colors = DEFAULT_COLORS,
  slideDuration = 90, // 3 seconds at 30fps
  transitionDuration = 30, // 1 second transition
  transitionStyle = 'crossfade',
  showTitles = true,
  logo,
}) => {
  const frame = useCurrentFrame();
  const { fps, durationInFrames } = useVideoConfig();

  // Ensure we have at least one item
  const safeItems = items.length > 0 ? items : DEFAULT_ITEMS;

  // Calculate which slide is currently active
  const totalSlideTime = slideDuration + transitionDuration;
  const currentSlideIndex = Math.floor(frame / totalSlideTime);
  const frameInSlide = frame % totalSlideTime;

  // Get current and next slide (guaranteed to exist due to safeItems)
  const currentSlide = safeItems[currentSlideIndex % safeItems.length]!;
  const nextSlide = safeItems[(currentSlideIndex + 1) % safeItems.length]!;

  // Calculate transition progress
  const isTransitioning = frameInSlide >= slideDuration;
  const transitionProgress = isTransitioning
    ? (frameInSlide - slideDuration) / transitionDuration
    : 0;

  // Calculate opacity for current slide
  const getCurrentSlideOpacity = () => {
    if (!isTransitioning) return 1;

    switch (transitionStyle) {
      case 'fade':
      case 'crossfade':
        return interpolate(transitionProgress, [0, 1], [1, 0]);
      case 'slide':
      case 'zoom':
        return 1;
      default:
        return 1;
    }
  };

  // Calculate opacity for next slide
  const getNextSlideOpacity = () => {
    if (!isTransitioning) return 0;

    switch (transitionStyle) {
      case 'fade':
        return interpolate(transitionProgress, [0, 1], [0, 1]);
      case 'crossfade':
        return interpolate(transitionProgress, [0, 1], [0, 1]);
      case 'slide':
      case 'zoom':
        return 1;
      default:
        return 0;
    }
  };

  // Calculate transform for current slide
  const getCurrentSlideTransform = () => {
    if (!isTransitioning) return 'none';

    switch (transitionStyle) {
      case 'slide':
        const translateX = interpolate(transitionProgress, [0, 1], [0, -100]);
        return `translateX(${translateX}%)`;
      case 'zoom':
        const scale = interpolate(transitionProgress, [0, 1], [1, 1.2]);
        return `scale(${scale})`;
      default:
        return 'none';
    }
  };

  // Calculate transform for next slide
  const getNextSlideTransform = () => {
    if (!isTransitioning) return 'translateX(100%)';

    switch (transitionStyle) {
      case 'slide':
        const translateX = interpolate(transitionProgress, [0, 1], [100, 0]);
        return `translateX(${translateX}%)`;
      case 'zoom':
        const scale = interpolate(transitionProgress, [0, 1], [0.8, 1]);
        return `scale(${scale})`;
      default:
        return 'none';
    }
  };

  // Title animation
  const titleOpacity = spring({
    frame: frameInSlide,
    fps,
    config: {
      damping: 200,
      stiffness: 100,
      mass: 0.5,
    },
  });

  const titleTranslateY = spring({
    frame: frameInSlide,
    fps,
    config: {
      damping: 200,
      stiffness: 100,
      mass: 0.5,
    },
    from: 30,
    to: 0,
  });

  // Render a slide
  const renderSlide = (
    slide: (typeof items)[0],
    opacity: number,
    transform: string,
    zIndex: number
  ) => (
    <AbsoluteFill
      style={{
        opacity,
        transform,
        zIndex,
        backgroundColor: colors.background,
        justifyContent: 'center',
        alignItems: 'center',
        overflow: 'hidden',
      }}
    >
      {/* Background Image */}
      {slide.image ? (
        <Img
          src={slide.image}
          style={{
            width: '100%',
            height: '100%',
            objectFit: 'cover',
          }}
        />
      ) : (
        <div
          style={{
            width: '100%',
            height: '100%',
            background: `linear-gradient(135deg, ${colors.primary} 0%, ${colors.secondary} 100%)`,
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
          }}
        >
          <span
            style={{
              fontSize: 120,
              color: 'rgba(255, 255, 255, 0.2)',
              fontWeight: 'bold',
            }}
          >
            {slide.title}
          </span>
        </div>
      )}

      {/* Title Overlay */}
      {showTitles && (
        <div
          style={{
            position: 'absolute',
            bottom: 0,
            left: 0,
            right: 0,
            padding: '60px',
            background:
              'linear-gradient(transparent, rgba(0, 0, 0, 0.7))',
          }}
        >
          <h2
            style={{
              fontSize: 48,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              marginBottom: 10,
              opacity: titleOpacity,
              transform: `translateY(${titleTranslateY}px)`,
            }}
          >
            {slide.title}
          </h2>
          {slide.description && (
            <p
              style={{
                fontSize: 24,
                color: 'rgba(255, 255, 255, 0.8)',
                margin: 0,
                opacity: titleOpacity,
                transform: `translateY(${titleTranslateY}px)`,
              }}
            >
              {slide.description}
            </p>
          )}
        </div>
      )}
    </AbsoluteFill>
  );

  return (
    <AbsoluteFill
      style={{
        backgroundColor: colors.background,
        fontFamily: 'system-ui, -apple-system, sans-serif',
        overflow: 'hidden',
      }}
    >
      {/* Current Slide */}
      {renderSlide(
        currentSlide,
        getCurrentSlideOpacity(),
        getCurrentSlideTransform(),
        1
      )}

      {/* Next Slide (during transition) */}
      {isTransitioning &&
        renderSlide(
          nextSlide,
          getNextSlideOpacity(),
          getNextSlideTransform(),
          2
        )}

      {/* Optional Logo */}
      {logo && (
        <img
          src={logo}
          alt="Logo"
          style={{
            position: 'absolute',
            top: 40,
            right: 40,
            height: 50,
            zIndex: 10,
            opacity: 0.9,
          }}
        />
      )}

      {/* Main Title (if provided) */}
      {title && (
        <div
          style={{
            position: 'absolute',
            top: 40,
            left: 40,
            zIndex: 10,
          }}
        >
          <h1
            style={{
              fontSize: 32,
              fontWeight: 'bold',
              color: colors.text,
              margin: 0,
              textShadow: '0 2px 4px rgba(0, 0, 0, 0.5)',
            }}
          >
            {title}
          </h1>
        </div>
      )}

      {/* Slide Indicator */}
      <div
        style={{
          position: 'absolute',
          bottom: 20,
          left: '50%',
          transform: 'translateX(-50%)',
          display: 'flex',
          gap: 10,
          zIndex: 10,
        }}
      >
        {safeItems.map((_, index) => (
          <div
            key={index}
            style={{
              width: 10,
              height: 10,
              borderRadius: '50%',
              backgroundColor:
                index === currentSlideIndex % safeItems.length
                  ? colors.primary
                  : 'rgba(255, 255, 255, 0.5)',
              transition: 'background-color 0.3s',
            }}
          />
        ))}
      </div>
    </AbsoluteFill>
  );
};

export default Slideshow;
