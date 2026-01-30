/**
 * Workspace Storage
 *
 * CRUD operations for workspaces.
 * Workspaces can be stored anywhere on disk via rootPath.
 * Default location: ~/.creator-flow/workspaces/
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  readdirSync,
  rmSync,
  statSync,
} from 'fs';
import { join } from 'path';
import { homedir } from 'os';
import { randomUUID } from 'crypto';
import { expandPath, toPortablePath } from '../utils/paths.ts';
import { getDefaultStatusConfig, saveStatusConfig, ensureDefaultIconFiles } from '../statuses/storage.ts';
import { getDefaultLabelConfig, saveLabelConfig } from '../labels/storage.ts';
import { getDefaultViews } from '../views/defaults.ts';
import { saveViewsConfig } from '../views/storage.ts';
import type { ViewsConfig } from '../views/storage.ts';
import { loadConfigDefaults } from '../config/storage.ts';
import { DEFAULT_MODEL } from '../config/models.ts';
import type {
  WorkspaceConfig,
  CreateWorkspaceInput,
  LoadedWorkspace,
  WorkspaceSummary,
} from './types.ts';

const CONFIG_DIR = join(homedir(), '.creator-flow');
const DEFAULT_WORKSPACES_DIR = join(CONFIG_DIR, 'workspaces');

/** Subdirectory name for workspace internal data (sessions, sources, skills, config) */
const WORKSPACE_DATA_DIR = '.creator-flow';

// ============================================================
// Path Utilities
// ============================================================

/**
 * Get the default workspaces directory (~/.creator-flow/workspaces/)
 */
export function getDefaultWorkspacesDir(): string {
  return DEFAULT_WORKSPACES_DIR;
}

/**
 * Ensure default workspaces directory exists
 */
export function ensureDefaultWorkspacesDir(): void {
  if (!existsSync(DEFAULT_WORKSPACES_DIR)) {
    mkdirSync(DEFAULT_WORKSPACES_DIR, { recursive: true });
  }
}

/**
 * Get workspace root path from ID
 * @param workspaceId - Workspace ID
 * @returns Absolute path to workspace root in default location
 */
export function getWorkspacePath(workspaceId: string): string {
  return join(DEFAULT_WORKSPACES_DIR, workspaceId);
}

/**
 * Get path to workspace data directory (.creator-flow)
 * @param rootPath - Absolute path to workspace root folder
 */
export function getWorkspaceDataPath(rootPath: string): string {
  return join(rootPath, WORKSPACE_DATA_DIR);
}

/**
 * Get path to workspace sources directory
 * @param rootPath - Absolute path to workspace root folder
 */
export function getWorkspaceSourcesPath(rootPath: string): string {
  return join(rootPath, WORKSPACE_DATA_DIR, 'sources');
}

/**
 * Get path to workspace sessions directory
 * @param rootPath - Absolute path to workspace root folder
 */
export function getWorkspaceSessionsPath(rootPath: string): string {
  return join(rootPath, WORKSPACE_DATA_DIR, 'sessions');
}

/**
 * Get path to workspace skills directory
 * @param rootPath - Absolute path to workspace root folder
 */
export function getWorkspaceSkillsPath(rootPath: string): string {
  return join(rootPath, WORKSPACE_DATA_DIR, 'skills');
}

// ============================================================
// Config Operations
// ============================================================

/**
 * Get path to workspace config.json
 * @param rootPath - Absolute path to workspace root folder
 */
function getWorkspaceConfigPath(rootPath: string): string {
  return join(rootPath, WORKSPACE_DATA_DIR, 'config.json');
}

/**
 * Load workspace config.json from a workspace folder
 * @param rootPath - Absolute path to workspace root folder
 */
export function loadWorkspaceConfig(rootPath: string): WorkspaceConfig | null {
  const configPath = getWorkspaceConfigPath(rootPath);
  if (!existsSync(configPath)) return null;

  try {
    const config = JSON.parse(readFileSync(configPath, 'utf-8')) as WorkspaceConfig;

    // Expand path variables in defaults for portability
    if (config.defaults?.workingDirectory) {
      config.defaults.workingDirectory = expandPath(config.defaults.workingDirectory);
    }

    return config;
  } catch {
    return null;
  }
}

/**
 * Save workspace config.json to a workspace folder
 * @param rootPath - Absolute path to workspace root folder
 */
export function saveWorkspaceConfig(rootPath: string, config: WorkspaceConfig): void {
  const dataPath = getWorkspaceDataPath(rootPath);
  if (!existsSync(dataPath)) {
    mkdirSync(dataPath, { recursive: true });
  }

  // Convert paths to portable form for cross-machine compatibility
  const storageConfig: WorkspaceConfig = {
    ...config,
    updatedAt: Date.now(),
  };

  if (storageConfig.defaults?.workingDirectory) {
    storageConfig.defaults = {
      ...storageConfig.defaults,
      workingDirectory: toPortablePath(storageConfig.defaults.workingDirectory),
    };
  }

  writeFileSync(getWorkspaceConfigPath(rootPath), JSON.stringify(storageConfig, null, 2));
}

// ============================================================
// Load Operations
// ============================================================

/**
 * Count subdirectories in a path
 */
function countSubdirs(dirPath: string): number {
  if (!existsSync(dirPath)) return 0;
  try {
    return readdirSync(dirPath, { withFileTypes: true }).filter((d) => d.isDirectory()).length;
  } catch {
    return 0;
  }
}

/**
 * List subdirectory names in a path
 */
function listSubdirNames(dirPath: string): string[] {
  if (!existsSync(dirPath)) return [];
  try {
    return readdirSync(dirPath, { withFileTypes: true })
      .filter((d) => d.isDirectory())
      .map((d) => d.name);
  } catch {
    return [];
  }
}

/**
 * Load workspace with summary info from a rootPath
 * @param rootPath - Absolute path to workspace root folder
 */
export function loadWorkspace(rootPath: string): LoadedWorkspace | null {
  const config = loadWorkspaceConfig(rootPath);
  if (!config) return null;

  // Ensure plugin manifest exists (migration for existing workspaces)
  ensurePluginManifest(rootPath, config.name);

  // Ensure skills directory exists (migration for existing workspaces)
  const skillsPath = getWorkspaceSkillsPath(rootPath);
  if (!existsSync(skillsPath)) {
    mkdirSync(skillsPath, { recursive: true });
  }

  return {
    config,
    sourceSlugs: listSubdirNames(getWorkspaceSourcesPath(rootPath)),
    sessionCount: countSubdirs(getWorkspaceSessionsPath(rootPath)),
  };
}

/**
 * Get workspace summary from a rootPath
 * @param rootPath - Absolute path to workspace root folder
 */
export function getWorkspaceSummary(rootPath: string): WorkspaceSummary | null {
  const config = loadWorkspaceConfig(rootPath);
  if (!config) return null;

  return {
    slug: config.slug,
    name: config.name,
    sourceCount: countSubdirs(getWorkspaceSourcesPath(rootPath)),
    sessionCount: countSubdirs(getWorkspaceSessionsPath(rootPath)),
    createdAt: config.createdAt,
    updatedAt: config.updatedAt,
  };
}

// ============================================================
// Create/Delete Operations
// ============================================================

/**
 * Generate URL-safe slug from name
 */
export function generateSlug(name: string): string {
  let slug = name
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '')
    .substring(0, 50);

  if (!slug) {
    slug = 'workspace';
  }

  return slug;
}

/**
 * Generate a unique folder path for a workspace by appending a numeric suffix
 * if the slug-based folder already exists.
 * E.g., "my-workspace", "my-workspace-2", "my-workspace-3", ...
 *
 * @param name - Display name to derive the slug from
 * @param baseDir - Parent directory where workspace folders live (e.g., ~/.creator-flow/workspaces/)
 * @returns Full path to a unique, non-existing folder
 */
export function generateUniqueWorkspacePath(name: string, baseDir: string): string {
  const slug = generateSlug(name);
  let candidate = join(baseDir, slug);

  if (!existsSync(candidate)) {
    return candidate;
  }

  // Append numeric suffix until we find a non-existing path
  let counter = 2;
  while (existsSync(join(baseDir, `${slug}-${counter}`))) {
    counter++;
  }

  return join(baseDir, `${slug}-${counter}`);
}

/**
 * Options for creating a workspace
 */
export interface CreateWorkspaceAtPathOptions {
  /** Optional default settings for new sessions */
  defaults?: WorkspaceConfig['defaults'];
  /** Skip initializing default labels (when app will provide custom labels) */
  skipDefaultLabels?: boolean;
  /** Skip initializing default statuses (when app will provide custom statuses) */
  skipDefaultStatuses?: boolean;
}

/**
 * Create workspace folder structure at a given path
 * @param rootPath - Absolute path where workspace folder will be created
 * @param name - Display name for the workspace
 * @param options - Optional configuration for workspace creation
 * @returns The created WorkspaceConfig
 */
export function createWorkspaceAtPath(
  rootPath: string,
  name: string,
  options?: CreateWorkspaceAtPathOptions | WorkspaceConfig['defaults']
): WorkspaceConfig {
  // Support legacy signature: createWorkspaceAtPath(path, name, defaults)
  const opts: CreateWorkspaceAtPathOptions = options && 'skipDefaultLabels' in options
    ? options
    : { defaults: options as WorkspaceConfig['defaults'] };
  const now = Date.now();
  const slug = generateSlug(name);

  // Load global defaults from config-defaults.json
  const globalDefaults = loadConfigDefaults();

  // Merge global defaults with provided defaults
  // Set workingDirectory to rootPath by default (workspace folder is the default working directory)
  const workspaceDefaults: WorkspaceConfig['defaults'] = {
    model: DEFAULT_MODEL,
    permissionMode: globalDefaults.workspaceDefaults.permissionMode,
    cyclablePermissionModes: globalDefaults.workspaceDefaults.cyclablePermissionModes,
    thinkingLevel: globalDefaults.workspaceDefaults.thinkingLevel,
    enabledSourceSlugs: [],
    workingDirectory: rootPath, // Default to workspace root path
    ...opts.defaults, // User-provided defaults override global defaults
  };

  const config: WorkspaceConfig = {
    id: `ws_${randomUUID().slice(0, 8)}`,
    name,
    slug,
    defaults: workspaceDefaults,
    localMcpServers: globalDefaults.workspaceDefaults.localMcpServers,
    createdAt: now,
    updatedAt: now,
  };

  // Create workspace directory structure (all under .creator-flow)
  mkdirSync(rootPath, { recursive: true });
  mkdirSync(getWorkspaceDataPath(rootPath), { recursive: true });
  mkdirSync(getWorkspaceSourcesPath(rootPath), { recursive: true });
  mkdirSync(getWorkspaceSessionsPath(rootPath), { recursive: true });
  mkdirSync(getWorkspaceSkillsPath(rootPath), { recursive: true });

  // Save config
  saveWorkspaceConfig(rootPath, config);

  // Initialize status configuration with defaults (skip if app will provide custom statuses)
  if (!opts.skipDefaultStatuses) {
    saveStatusConfig(rootPath, getDefaultStatusConfig());
    ensureDefaultIconFiles(rootPath);
  }

  // Initialize label configuration with defaults (skip if app will provide custom labels)
  if (!opts.skipDefaultLabels) {
    saveLabelConfig(rootPath, getDefaultLabelConfig());
  }

  // Initialize views configuration with defaults
  const defaultViewsConfig: ViewsConfig = { version: 1, views: getDefaultViews() };
  saveViewsConfig(rootPath, defaultViewsConfig);

  // Initialize plugin manifest for SDK integration (enables skills, commands, agents)
  ensurePluginManifest(rootPath, name);

  return config;
}

/**
 * Delete a workspace data folder (.creator-flow) and its contents.
 * Only removes the .creator-flow subdirectory, preserving user's project files.
 * @param rootPath - Absolute path to workspace root folder
 */
export function deleteWorkspaceFolder(rootPath: string): boolean {
  const dataPath = getWorkspaceDataPath(rootPath);
  if (!existsSync(dataPath)) return false;

  try {
    // Only delete the .creator-flow data directory, not the user's project folder
    rmSync(dataPath, { recursive: true });
    return true;
  } catch {
    return false;
  }
}

/**
 * Check if a valid workspace exists at a path
 * @param rootPath - Absolute path to check
 */
export function isValidWorkspace(rootPath: string): boolean {
  return existsSync(getWorkspaceConfigPath(rootPath));
}

/**
 * Rename a workspace (updates config.json in the workspace folder)
 * @param rootPath - Absolute path to workspace root folder
 * @param newName - New display name
 */
export function renameWorkspaceFolder(rootPath: string, newName: string): boolean {
  const config = loadWorkspaceConfig(rootPath);
  if (!config) return false;

  config.name = newName.trim();
  saveWorkspaceConfig(rootPath, config);
  return true;
}

// ============================================================
// Auto-Discovery (for default workspace location)
// ============================================================

/**
 * Discover workspace folders in the default location that have valid config.json
 * Returns paths to valid workspaces found in ~/.creator-flow/workspaces/
 */
export function discoverWorkspacesInDefaultLocation(): string[] {
  const discovered: string[] = [];

  if (!existsSync(DEFAULT_WORKSPACES_DIR)) {
    return discovered;
  }

  try {
    const entries = readdirSync(DEFAULT_WORKSPACES_DIR, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const rootPath = join(DEFAULT_WORKSPACES_DIR, entry.name);
      if (isValidWorkspace(rootPath)) {
        discovered.push(rootPath);
      }
    }
  } catch {
    // Ignore errors scanning directory
  }

  return discovered;
}

// ============================================================
// Local MCP Configuration
// ============================================================

/**
 * Check if local (stdio) MCP servers are enabled for a workspace.
 * Resolution order: ENV (CRAFT_LOCAL_MCP_ENABLED) > workspace config > default (true)
 *
 * @param rootPath - Absolute path to workspace root folder
 * @returns true if local MCP servers should be enabled
 */
export function isLocalMcpEnabled(rootPath: string): boolean {
  // 1. Environment variable override (highest priority)
  const envValue = process.env.CRAFT_LOCAL_MCP_ENABLED;
  if (envValue !== undefined) {
    return envValue.toLowerCase() === 'true';
  }

  // 2. Workspace config
  const config = loadWorkspaceConfig(rootPath);
  if (config?.localMcpServers?.enabled !== undefined) {
    return config.localMcpServers.enabled;
  }

  // 3. Default: enabled
  return true;
}

// ============================================================
// Exports
// ============================================================

// ============================================================
// Plugin Manifest (for SDK plugin integration)
// ============================================================

/**
 * Ensure workspace has a .claude-plugin/plugin.json manifest.
 * This allows the workspace to be loaded as an SDK plugin,
 * enabling skills, commands, and agents from the workspace.
 *
 * @param rootPath - Absolute path to workspace root folder
 * @param workspaceName - Display name for the workspace (used in plugin name)
 */
export function ensurePluginManifest(rootPath: string, workspaceName: string): void {
  const dataPath = getWorkspaceDataPath(rootPath);
  const pluginDir = join(dataPath, '.claude-plugin');
  const manifestPath = join(pluginDir, 'plugin.json');

  if (existsSync(manifestPath)) return;

  // Create .claude-plugin directory
  if (!existsSync(pluginDir)) {
    mkdirSync(pluginDir, { recursive: true });
  }

  // Create minimal plugin manifest
  const manifest = {
    name: `craft-workspace-${workspaceName.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`,
    version: '1.0.0',
  };

  writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
}

export { CONFIG_DIR, DEFAULT_WORKSPACES_DIR };
