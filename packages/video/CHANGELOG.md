# Changelog

All notable changes to the @creator-flow/video package will be documented in this file.

## [Unreleased]

### Added
- **Dynamic Video Duration**: Implemented `calculateMetadata` for 7 compositions (Slideshow, DataVisualization variants, ProductShowcase variants) to automatically calculate video duration based on content length
- **Improved Organization**: Added `<Folder>` structure in Root.tsx to group compositions by category (Basic Animations, Data Visualization, Product Showcase)
- **Reusable Components**: Created 4 new sub-components for better code organization:
  - `TitleText` - Animated title text component
  - `SubtitleText` - Animated subtitle text component
  - `SlideItem` - Individual slide item for Slideshow
  - `BarItem` - Animated bar chart item for DataVisualization

### Changed
- **Resource Management**: All image paths now use `staticFile()` to ensure correct resolution in both development and production environments
- **Type Safety**: Converted all Props definitions from `interface` to `type` for better type composition and consistency
- **Type Exports**: Fixed `DataPoint` type duplication by importing from DataVisualization composition

### Fixed
- **Deployment Issues**: Fixed asset path resolution in bundled/deployed environments using Remotion's `staticFile()` API
- **Type Conflicts**: Resolved duplicate `DataPoint` export between compositions and components

### Performance
- All 192 tests passing with no performance regression
- Test execution time: ~295ms

## [0.1.0] - 2025-01-XX

### Added
- Initial release with 9 video compositions
- Template system for video generation
- MCP server integration for AI-driven video creation
- Comprehensive test suite with 192 tests
