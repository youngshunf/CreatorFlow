/**
 * Workspace Types
 *
 * Workspaces are the top-level organizational unit. Everything (sources, sessions)
 * is scoped to a workspace.
 *
 * Directory structure (all workspace data is under .sprouty-ai/ to avoid mixing with user files):
 * {workspaceRoot}/
 *   └── .sprouty-ai/       - All workspace data (hidden from user's project files)
 *       ├── config.json      - Workspace settings
 *       ├── sources/         - Data sources (MCP, API, local)
 *       ├── sessions/        - Conversation sessions
 *       ├── skills/          - Agent skills
 *       ├── statuses/        - Status configuration and icons
 *       ├── labels/          - Label configuration
 *       └── permissions.json - Workspace-level permissions
 */

import type { PermissionMode } from '../agent/mode-manager.ts';
import type { ThinkingLevel } from '../agent/thinking-levels.ts';

/**
 * Local MCP server configuration
 * Controls whether stdio-based (local subprocess) MCP servers can be spawned.
 */
export interface LocalMcpConfig {
  /**
   * Whether local (stdio) MCP servers are enabled for this workspace.
   * When false, only HTTP-based MCP servers will be used.
   * Default: true (can be overridden by CRAFT_LOCAL_MCP_ENABLED env var)
   */
  enabled: boolean;
}

/**
 * Workspace configuration (stored in config.json)
 */
export interface WorkspaceConfig {
  id: string;
  name: string;
  slug: string; // Folder name (URL-safe)

  /**
   * Default settings for new sessions in this workspace
   */
  defaults?: {
    model?: string;
    enabledSourceSlugs?: string[]; // Sources to enable by default
    permissionMode?: PermissionMode; // Default permission mode ('safe', 'ask', 'allow-all')
    cyclablePermissionModes?: PermissionMode[]; // Which modes can be cycled with SHIFT+TAB (min 2, default: all 3)
    workingDirectory?: string;
    thinkingLevel?: ThinkingLevel; // Default thinking level ('off', 'think', 'max') - default: 'think'
    colorTheme?: string; // Color theme override for this workspace (preset ID). Undefined = inherit from app default.
  };

  /**
   * Local MCP server configuration.
   * Controls whether stdio-based MCP servers can be spawned in this workspace.
   * Resolution order: ENV (CRAFT_LOCAL_MCP_ENABLED) > workspace config > default (true)
   */
  localMcpServers?: LocalMcpConfig;

  /**
   * Application framework fields
   */
  
  /**
   * ID of the main application bound to this workspace.
   * Example: 'app.creator-media' for the content creation app.
   * When set, workspace initialization uses the app's configuration.
   */
  appId?: string;

  /**
   * List of installed plugin application IDs.
   * Plugin apps extend the main app with additional capabilities.
   * Example: ['plugin.data-analytics', 'plugin.seo-optimizer']
   */
  installedPluginApps?: string[];

  /**
   * Application-specific settings.
   * Each app can store its own configuration here.
   * Structure depends on the app's requirements.
   */
  appSettings?: Record<string, any>;

  createdAt: number;
  updatedAt: number;
}

/**
 * Workspace creation input
 */
export interface CreateWorkspaceInput {
  name: string;
  defaults?: WorkspaceConfig['defaults'];
}

/**
 * Loaded workspace with resolved sources
 */
export interface LoadedWorkspace {
  config: WorkspaceConfig;
  sourceSlugs: string[]; // Available source slugs (not fully loaded to save memory)
  sessionCount: number; // Number of sessions
}

/**
 * Workspace summary for listing (lightweight)
 */
export interface WorkspaceSummary {
  slug: string;
  name: string;
  sourceCount: number;
  sessionCount: number;
  createdAt: number;
  updatedAt: number;
}
