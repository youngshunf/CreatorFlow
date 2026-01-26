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
 * Preset data for workspace initialization
 */
export interface PresetData {
  /** Labels to create in the workspace */
  labels?: Array<{
    name: string;
    color: string;
  }>;
  /** Statuses to create in the workspace */
  statuses?: Array<{
    name: string;
    type: string;
  }>;
  /** Custom preset data (app-specific) */
  [key: string]: unknown;
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
    defaultModel?: string;
    autoSave?: boolean;
    autoBackup?: boolean;
    backupInterval?: number;
    theme?: string;
    language?: string;
    [key: string]: unknown;
  };
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
 * Applications are installed globally at ~/.creator-flow/apps/{app-id}/
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
