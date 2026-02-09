import React from 'react'

const ProjectDashboard = React.lazy(() => import('./ProjectDashboard'))
const CreationWorkspace = React.lazy(() => import('./CreationWorkspace'))
const HotTopicsBoard = React.lazy(() => import('./HotTopicsBoard'))
const ContentCalendar = React.lazy(() => import('./ContentCalendar'))
const AnalyticsDashboard = React.lazy(() => import('./AnalyticsDashboard'))

/**
 * APP 视图注册表
 * key: appId → Record<viewId, LazyComponent>
 */
export const APP_VIEW_REGISTRY: Record<string, Record<string, React.LazyExoticComponent<React.ComponentType>>> = {
  'app.creator-media': {
    dashboard: ProjectDashboard,
    workspace: CreationWorkspace,
    'hot-topics': HotTopicsBoard,
    calendar: ContentCalendar,
    analytics: AnalyticsDashboard,
  },
}
