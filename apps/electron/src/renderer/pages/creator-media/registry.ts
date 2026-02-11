import React from 'react'

const ProjectDashboard = React.lazy(() => import('./ProjectDashboard'))
const CreationWorkspace = React.lazy(() => import('./CreationWorkspace'))
const HotTopicsBoard = React.lazy(() => import('./HotTopicsBoard'))
const VideoStudio = React.lazy(() => import('./VideoStudio'))
const ContentCalendar = React.lazy(() => import('./ContentCalendar'))
const AnalyticsDashboard = React.lazy(() => import('./AnalyticsDashboard'))
const PlatformAccountsView = React.lazy(() => import('./PlatformAccountsView'))
const PublishWorkbench = React.lazy(() => import('./PublishWorkbench'))
const ScheduledTasksView = React.lazy(() => import('./ScheduledTasksView'))

/**
 * APP 视图注册表
 * key: appId → Record<viewId, LazyComponent>
 */
export const APP_VIEW_REGISTRY: Record<string, Record<string, React.LazyExoticComponent<React.ComponentType>>> = {
  'app.creator-media': {
    dashboard: ProjectDashboard,
    workspace: CreationWorkspace,
    'hot-topics': HotTopicsBoard,
    'video-studio': VideoStudio,
    calendar: ContentCalendar,
    analytics: AnalyticsDashboard,
    accounts: PlatformAccountsView,
    publisher: PublishWorkbench,
    'scheduled-tasks': ScheduledTasksView,
  },
}
