/**
 * Lucide 图标名称映射工具
 *
 * 将 kebab-case 图标名称映射到 Lucide React 组件，
 * 用于从 manifest.json 的 views.sidebar 配置动态加载图标。
 */
import {
  LayoutDashboard,
  PenTool,
  TrendingUp,
  Send,
  Calendar,
  BarChart3,
  Users,
  Clapperboard,
  Home,
  Inbox,
  Settings,
  Files,
  Store,
  Film,
  Globe,
  Search,
  Plus,
  Zap,
  Flag,
  Tag,
  FolderOpen,
  ExternalLink,
  HelpCircle,
  CheckCircle2,
  type LucideIcon,
} from 'lucide-react'

const ICON_MAP: Record<string, LucideIcon> = {
  'layout-dashboard': LayoutDashboard,
  'pen-tool': PenTool,
  'trending-up': TrendingUp,
  'send': Send,
  'calendar': Calendar,
  'bar-chart-3': BarChart3,
  'users': Users,
  'clapperboard': Clapperboard,
  'home': Home,
  'inbox': Inbox,
  'settings': Settings,
  'files': Files,
  'store': Store,
  'film': Film,
  'globe': Globe,
  'search': Search,
  'plus': Plus,
  'zap': Zap,
  'flag': Flag,
  'tag': Tag,
  'folder-open': FolderOpen,
  'external-link': ExternalLink,
  'help-circle': HelpCircle,
  'check-circle-2': CheckCircle2,
}

/**
 * 根据 kebab-case 图标名称获取对应的 Lucide 图标组件
 */
export function getLucideIcon(name: string): LucideIcon | undefined {
  return ICON_MAP[name]
}
