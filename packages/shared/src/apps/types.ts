/**
 * Application Types
 *
 * Applications are pre-configured bundles that provide complete solutions for specific use cases.
 * An application combines:
 * - Skills: AI capabilities (using Claude Agent Skill format)
 * - Platform dependencies: UI components, data models, integrations
 * - Workspace configuration: directory structure, default settings, preset data
 *
 * Applications follow the "capability combiner" pattern - they don't implement features,
 * they declare dependencies and configuration.
 */

import type { PermissionMode } from '../agent/mode-manager.ts';
import type { ThinkingLevel } from '../agent/thinking-levels.ts';
import type { LabelConfig } from '../labels/types.ts';
import type { StatusConfig, StatusCategory } from '../statuses/types.ts';

/**
 * Application pricing model
 */
export interface AppPricing {
  type: 'free' | 'paid' | 'freemium';
  /** One-time purchase price (if type is 'paid') */
  price?: number;
  /** Subscription pricing (if applicable) */
  subscription?: {
    monthly: number;
    yearly: number;
  };
  /** Trial configuration */
  trial?: {
    enabled: boolean;
    days: number;
  };
}

/**
 * Application author information
 */
export interface AppAuthor {
  name: string;
  email?: string;
  url?: string;
}

/**
 * Directory structure template for workspace initialization
 * Keys are directory names, values are nested structures or empty objects for leaf directories
 */
export interface DirectoryStructure {
  [key: string]: DirectoryStructure;
}

/**
 * Preset label configuration for app manifest.
 * Simplified version of LabelConfig - id and children are optional.
 * If id is omitted, it will be generated from name.
 */
export interface PresetLabelConfig {
  /** Unique ID (if omitted, generated from name as slug) */
  id?: string;
  /** Display name */
  name: string;
  /** Color in EntityColor format (hex string or { light, dark } object) */
  color?: string | { light: string; dark: string };
  /** Child labels */
  children?: PresetLabelConfig[];
  /** Value type for valued labels */
  valueType?: 'string' | 'number' | 'date';
}

/**
 * Preset status configuration for app manifest.
 * Simplified version of StatusConfig with sensible defaults.
 */
export interface PresetStatusConfig {
  /** Unique ID (if omitted, generated from label as slug) */
  id?: string;
  /** Display name */
  label: string;
  /** Color in EntityColor format */
  color?: string | { light: string; dark: string };
  /** Icon (emoji or URL) */
  icon?: string;
  /** Category: 'open' (inbox) or 'closed' (archive). Default: 'open' */
  category?: StatusCategory;
  /** If true, cannot be deleted. Default: false */
  isFixed?: boolean;
  /** Display order (lower = first). If omitted, uses array index */
  order?: number;
}

/**
 * Preset data for workspace initialization
 */
export interface PresetData {
  /**
   * Labels to create in the workspace.
   * These will be merged with/replace default labels.
   * If 'replaceDefaults' is true, these replace the default labels entirely.
   */
  labels?: PresetLabelConfig[];
  /** If true, preset labels replace defaults instead of merging. Default: false */
  replaceDefaultLabels?: boolean;
  
  /**
   * Statuses to create in the workspace.
   * These will be merged with default statuses.
   * Fixed statuses (todo, done, cancelled) cannot be replaced.
   */
  statuses?: PresetStatusConfig[];
  /** If true, preset statuses replace non-fixed defaults. Default: false */
  replaceDefaultStatuses?: boolean;
  
  /** Custom preset data (app-specific) */
  [key: string]: unknown;
}

/**
 * Local MCP server configuration for workspace
 */
export interface WorkspaceLocalMcpConfig {
  /**
   * Whether local (stdio) MCP servers are enabled for this workspace.
   * When false, only HTTP-based MCP servers will be used.
   * Default: true (can be overridden by CRAFT_LOCAL_MCP_ENABLED env var)
   */
  enabled: boolean;
}

/**
 * Workspace initialization configuration
 */
export interface WorkspaceInitConfig {
  /** Default directory for new workspaces (can use ~ for home) */
  defaultDirectory?: string;
  /** Directory structure template to create */
  directoryStructure?: DirectoryStructure;
  /** Default settings for new sessions */
  defaultSettings?: {
    // AI model configuration
    defaultModel?: string;
    
    // Permission configuration
    /** Default permission mode ('safe', 'ask', 'allow-all'). Default: from global config */
    permissionMode?: PermissionMode;
    /** Which modes can be cycled with SHIFT+TAB (min 2). Default: ['safe', 'ask', 'allow-all'] */
    cyclablePermissionModes?: PermissionMode[];
    
    // Thinking configuration
    /** Default thinking level ('off', 'think', 'max'). Default: 'think' */
    thinkingLevel?: ThinkingLevel;
    
    // File management
    autoSave?: boolean;
    autoBackup?: boolean;
    backupInterval?: number;
    
    // UI preferences
    theme?: string;
    language?: string;
    
    [key: string]: unknown;
  };
  
  /**
   * Local MCP server configuration.
   * Controls whether stdio-based MCP servers can be spawned in workspaces using this app.
   */
  localMcpServers?: WorkspaceLocalMcpConfig;
  
  /** Preset data to initialize */
  presetData?: PresetData;
}

/**
 * Marketplace metadata for app discovery
 */
export interface MarketplaceInfo {
  /** Category for marketplace browsing */
  category: string;
  /** Tags for search */
  tags: string[];
  /** Screenshot paths (relative to app directory) */
  screenshots?: string[];
  /** Demo video path */
  video?: string;
  /** Changelog in markdown format */
  changelog?: string;
}

/**
 * Plugin compatibility information
 */
export interface PluginCompatibility {
  /** List of compatible plugin app IDs */
  compatible?: string[];
  /** Recommended plugins */
  recommended?: string[];
}

/**
 * OS compatibility
 */
export type OSType = 'macos' | 'windows' | 'linux';

/**
 * Application compatibility requirements
 */
export interface AppCompatibility {
  /** Minimum platform version */
  platform: string;
  /** Supported operating systems */
  os?: OSType[];
  /** Node.js version requirement */
  node?: string;
}

/**
 * Application dependencies
 */
export interface AppDependencies {
  /** Platform version requirement */
  platform?: string;
  /** UI component dependencies (format: "ui.component-name@version") */
  ui?: string[];
  /** Data model dependencies (format: "model.model-name@version") */
  models?: string[];
  /** Integration dependencies (format: "integration.name@version") */
  integrations?: string[];
}

/**
 * Application capabilities (Skills binding)
 */
export interface AppCapabilities {
  /**
   * Skills bound to this application.
   * Format: "skill-slug@version" or just "skill-slug" for bundled skills.
   * These skills will be copied/installed to workspaces using this app.
   */
  skills?: string[];
}

/**
 * Application Manifest
 *
 * The complete definition of an application, stored as manifest.json in the app directory.
 * Applications are installed globally at ~/.sprouty-ai/apps/{app-id}/
 */
export interface AppManifest {
  /** Unique application ID (format: "app.name" for main apps, "plugin.name" for plugins) */
  id: string;
  /** Display name */
  name: string;
  /** Semantic version */
  version: string;
  /** Application type: 'app' for main applications, 'plugin' for extensions */
  type?: 'app' | 'plugin';

  // Metadata
  author?: AppAuthor;
  description?: string;
  license?: string;
  homepage?: string;
  repository?: string;

  // Pricing
  pricing?: AppPricing;

  // Dependencies and Capabilities
  dependencies?: AppDependencies;
  capabilities?: AppCapabilities;

  // Workspace Configuration
  workspace?: WorkspaceInitConfig;

  // Permissions
  permissions?: string[];

  // Compatibility
  compatibility?: AppCompatibility;

  // Marketplace
  marketplace?: MarketplaceInfo;

  // Plugin support
  plugins?: PluginCompatibility;

  // APP 视图配置
  /** APP 视图配置（侧边栏导航项等） */
  views?: {
    /** 默认视图 ID */
    defaultView?: string;
    /** 侧边栏视图列表 */
    sidebar?: Array<{
      viewId: string;
      title: string;
      icon: string;    // lucide icon name, e.g. "layout-dashboard", "send"
      order: number;
    }>;
  };

  // Plugin-specific fields (only for type: 'plugin')
  /** For plugins: list of compatible main app IDs with versions */
  compatibleApps?: string[];
  /** For plugins: UI extensions */
  extensions?: {
    menu?: Array<{
      id: string;
      label: string;
      icon?: string;
      path?: string;
    }>;
    toolbar?: Array<{
      id: string;
      label: string;
      icon?: string;
    }>;
  };
}

/**
 * Loaded application with resolved paths
 */
export interface LoadedApp {
  /** The manifest data */
  manifest: AppManifest;
  /** Absolute path to the app directory */
  path: string;
  /** Absolute path to icon file (if exists) */
  iconPath?: string;
  /** Whether this is a bundled (built-in) app */
  bundled?: boolean;
}

/**
 * Application summary for listing (lightweight)
 */
export interface AppSummary {
  id: string;
  name: string;
  description?: string;
  version: string;
  type: 'app' | 'plugin';
  iconPath?: string;
  bundled?: boolean;
}

/**
 * Application installation options
 */
export interface InstallAppOptions {
  /** Whether to install bundled skills to the app directory */
  installSkills?: boolean;
  /** Force reinstall even if already installed */
  force?: boolean;
}

/**
 * Skill reference in app manifest
 */
export interface SkillReference {
  /** Skill slug or ID */
  slug: string;
  /** Optional version constraint */
  version?: string;
  /** Whether this skill is bundled with the app */
  bundled?: boolean;
  /** Path to bundled skill (relative to app directory) */
  bundledPath?: string;
}
