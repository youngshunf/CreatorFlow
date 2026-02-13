/**
 * Workspace Initializer
 *
 * Initialize workspaces based on application configuration.
 * Handles directory structure creation, skill installation, and preset data.
 */

import {
  existsSync,
  mkdirSync,
  copyFileSync,
  readdirSync,
  readFileSync,
  writeFileSync,
} from 'fs';
import { join } from 'path';
import type { DirectoryStructure } from './types.ts';
import type { WorkspaceConfig } from '../workspaces/types.ts';
import { installSkill } from '../marketplace/skill-installer.ts';
import {
  createWorkspaceAtPath,
  getWorkspaceSkillsPath,
  getWorkspaceDataPath,
  saveWorkspaceConfig,
  loadWorkspaceConfig,
} from '../workspaces/storage.ts';
import { loadAppById, getAppPath } from './storage.ts';
import { debug } from '../utils/debug.ts';
import { resolveSourceConfigPaths } from './video-mcp-paths.ts';

// ============================================================
// Directory Structure Creation
// ============================================================

/**
 * Recursively create directory structure from template
 */
export function createDirectoryStructure(basePath: string, structure: DirectoryStructure): void {
  for (const [name, children] of Object.entries(structure)) {
    const dirPath = join(basePath, name);
    
    if (!existsSync(dirPath)) {
      mkdirSync(dirPath, { recursive: true });
    }
    
    // Recursively create children
    if (children && Object.keys(children).length > 0) {
      createDirectoryStructure(dirPath, children);
    }
  }
}

// ============================================================
// Skill Installation
// ============================================================

/**
 * Copy a skill directory from app to workspace
 */
function copySkillToWorkspace(
  skillSourcePath: string,
  workspaceSkillsPath: string,
  skillSlug: string
): boolean {
  const targetPath = join(workspaceSkillsPath, skillSlug);
  
  if (!existsSync(skillSourcePath)) {
    return false;
  }
  
  // Create target directory
  mkdirSync(targetPath, { recursive: true });
  
  // Copy all files from source to target
  try {
    const entries = readdirSync(skillSourcePath, { withFileTypes: true });
    for (const entry of entries) {
      const sourcePath = join(skillSourcePath, entry.name);
      const destPath = join(targetPath, entry.name);
      
      if (entry.isFile()) {
        copyFileSync(sourcePath, destPath);
      } else if (entry.isDirectory()) {
        // Recursively copy directories
        copyDirectory(sourcePath, destPath);
      }
    }
    return true;
  } catch {
    return false;
  }
}

/**
 * Recursively copy a directory
 */
function copyDirectory(source: string, target: string): void {
  mkdirSync(target, { recursive: true });
  
  const entries = readdirSync(source, { withFileTypes: true });
  for (const entry of entries) {
    const sourcePath = join(source, entry.name);
    const destPath = join(target, entry.name);
    
    if (entry.isFile()) {
      copyFileSync(sourcePath, destPath);
    } else if (entry.isDirectory()) {
      copyDirectory(sourcePath, destPath);
    }
  }
}

/**
 * Install bundled skills from app to workspace
 * First tries to download from cloud, falls back to local bundled skills
 */
function installBundledSkills(
  appPath: string,
  workspaceRoot: string,
  skillRefs: string[]
): { installed: string[]; failed: string[] } {
  const workspaceSkillsPath = getWorkspaceSkillsPath(workspaceRoot);
  const appSkillsPath = join(appPath, 'skills');
  
  const installed: string[] = [];
  const failed: string[] = [];
  
  // Ensure workspace skills directory exists
  if (!existsSync(workspaceSkillsPath)) {
    mkdirSync(workspaceSkillsPath, { recursive: true });
  }
  
  for (const skillRef of skillRefs) {
    // Parse skill reference (format: "skill-slug" or "skill-slug@version")
    const [skillSlug, version] = skillRef.split('@');
    if (!skillSlug) {
      failed.push(skillRef);
      continue;
    }
    
    const skillSourcePath = join(appSkillsPath, skillSlug);
    
    if (existsSync(skillSourcePath)) {
      const success = copySkillToWorkspace(skillSourcePath, workspaceSkillsPath, skillSlug);
      if (success) {
        installed.push(skillSlug);
      } else {
        failed.push(skillSlug);
      }
    } else {
      // Skill not bundled, might be from marketplace (skip for now)
      failed.push(skillSlug);
    }
  }
  
  return { installed, failed };
}

/**
 * Install skills from cloud (marketplace)
 * Used for both bundled and marketplace apps to ensure latest versions
 */
export async function installSkillsFromCloud(
  workspaceRoot: string,
  skillRefs: string[],
  fallbackAppPath?: string
): Promise<{ installed: string[]; failed: string[] }> {
  const workspaceSkillsPath = getWorkspaceSkillsPath(workspaceRoot);
  const appSkillsPath = fallbackAppPath ? join(fallbackAppPath, 'skills') : null;
  
  const installed: string[] = [];
  const failed: string[] = [];
  
  // Ensure workspace skills directory exists
  if (!existsSync(workspaceSkillsPath)) {
    mkdirSync(workspaceSkillsPath, { recursive: true });
  }
  
  for (const skillRef of skillRefs) {
    // Parse skill reference (format: "skill-slug" or "skill-slug@version")
    const [skillSlug, version] = skillRef.split('@');
    if (!skillSlug) {
      failed.push(skillRef);
      continue;
    }
    
    // Try to download from cloud first
    try {
      debug(`[installSkillsFromCloud] Downloading skill ${skillSlug} from cloud...`);
      const result = await installSkill(
        workspaceRoot,
        skillSlug,
        version || 'latest'
      );
      
      if (result.success) {
        installed.push(skillSlug);
        debug(`[installSkillsFromCloud] Successfully installed ${skillSlug} from cloud`);
        continue;
      }
    } catch (error) {
      debug(`[installSkillsFromCloud] Failed to download ${skillSlug} from cloud: ${error}`);
    }
    
    // Fall back to local bundled skill if cloud download fails
    if (appSkillsPath) {
      const skillSourcePath = join(appSkillsPath, skillSlug);
      if (existsSync(skillSourcePath)) {
        debug(`[installSkillsFromCloud] Falling back to bundled skill for ${skillSlug}`);
        const success = copySkillToWorkspace(skillSourcePath, workspaceSkillsPath, skillSlug);
        if (success) {
          installed.push(skillSlug);
          continue;
        }
      }
    }
    
    failed.push(skillSlug);
  }
  
  return { installed, failed };
}

// ============================================================
// App Data Directory Copy (Labels & Statuses)
// ============================================================

/**
 * Check if app has custom labels configuration.
 * Returns the app path if found, null otherwise.
 */
function getAppPathWithCustomLabels(appId: string): string | null {
  const appPath = getAppPath(appId);
  if (existsSync(appPath) && existsSync(join(appPath, 'labels', 'config.json'))) {
    return appPath;
  }

  return null;
}

/**
 * Check if app has custom statuses configuration.
 * Returns the app path if found, null otherwise.
 */
function getAppPathWithCustomStatuses(appId: string): string | null {
  const appPath = getAppPath(appId);
  if (existsSync(appPath) && existsSync(join(appPath, 'statuses', 'config.json'))) {
    return appPath;
  }

  return null;
}

/**
 * Copy labels directory from app to workspace.
 * If the app has labels/config.json, copy it to workspace.
 * Otherwise, workspace will use default labels.
 */
function copyAppLabelsToWorkspace(appPath: string, workspaceRoot: string): boolean {
  const appLabelsConfig = join(appPath, 'labels', 'config.json');
  
  if (!existsSync(appLabelsConfig)) {
    debug(`[copyAppLabelsToWorkspace] No labels/config.json in app, using defaults`);
    return false;
  }
  
  const workspaceLabelsDir = join(workspaceRoot, '.sprouty-ai', 'labels');
  const workspaceLabelsConfig = join(workspaceLabelsDir, 'config.json');
  
  // Ensure workspace labels directory exists
  if (!existsSync(workspaceLabelsDir)) {
    mkdirSync(workspaceLabelsDir, { recursive: true });
  }
  
  // Copy config file
  copyFileSync(appLabelsConfig, workspaceLabelsConfig);
  debug(`[copyAppLabelsToWorkspace] Copied labels config to ${workspaceLabelsConfig}`);
  
  return true;
}

/**
 * Copy statuses directory from app to workspace.
 * If the app has statuses/config.json, copy it along with any icons.
 * Otherwise, workspace will use default statuses.
 */
function copyAppStatusesToWorkspace(appPath: string, workspaceRoot: string): boolean {
  const appStatusesConfig = join(appPath, 'statuses', 'config.json');
  
  if (!existsSync(appStatusesConfig)) {
    debug(`[copyAppStatusesToWorkspace] No statuses/config.json in app, using defaults`);
    return false;
  }
  
  const workspaceStatusesDir = join(workspaceRoot, '.sprouty-ai', 'statuses');
  const workspaceStatusesConfig = join(workspaceStatusesDir, 'config.json');
  const workspaceIconsDir = join(workspaceStatusesDir, 'icons');
  
  // Ensure workspace statuses directory exists
  if (!existsSync(workspaceStatusesDir)) {
    mkdirSync(workspaceStatusesDir, { recursive: true });
  }
  
  // Copy config file
  copyFileSync(appStatusesConfig, workspaceStatusesConfig);
  debug(`[copyAppStatusesToWorkspace] Copied statuses config to ${workspaceStatusesConfig}`);
  
  // Copy icons directory if exists
  const appIconsDir = join(appPath, 'statuses', 'icons');
  if (existsSync(appIconsDir)) {
    copyDirectory(appIconsDir, workspaceIconsDir);
    debug(`[copyAppStatusesToWorkspace] Copied status icons to ${workspaceIconsDir}`);
  }
  
  return true;
}

/**
 * Copy app-specific labels, statuses, and sources to workspace.
 * This is the simplified approach: directly copy directory contents.
 *
 * @param appId - Application ID
 * @param workspaceRoot - Workspace root path
 * @param isPackaged - Whether the app is packaged (default: false)
 * @param resourcesPath - Resources path for packaged app
 */
function copyAppDataToWorkspace(appId: string, workspaceRoot: string): void {
  const appPath = getAppPath(appId);

  if (!existsSync(appPath)) {
    debug(`[copyAppDataToWorkspace] App path not found: ${appPath}`);
    return;
  }

  // Copy labels (if app has custom labels)
  copyAppLabelsToWorkspace(appPath, workspaceRoot);

  // Copy statuses (if app has custom statuses)
  copyAppStatusesToWorkspace(appPath, workspaceRoot);

  // Copy sources (if app has preset sources)
  copyAppSourcesToWorkspace(appPath, workspaceRoot, isPackaged, resourcesPath);
}

// ============================================================
// Sources Directory Copy
// ============================================================

/**
 * Check if app has sources configuration.
 * Returns the app path if found, null otherwise.
 */
function getAppPathWithSources(appId: string): string | null {
  const appPath = getAppPath(appId);
  if (existsSync(appPath) && existsSync(join(appPath, 'sources'))) {
    return appPath;
  }

  return null;
}

/**
 * Copy sources directory from app to workspace.
 * Each source is a subdirectory with config.json and optionally guide.md.
 * Existing sources in the workspace are NOT overwritten.
 *
 * Path placeholders in config.json are resolved:
 * - {{BUN_PATH}} -> Bun executable path
 * - {{VIDEO_MCP_SERVER_PATH}} -> Video MCP server entry file path
 */
function copyAppSourcesToWorkspace(
  appPath: string,
  workspaceRoot: string,
  isPackaged: boolean = false,
  resourcesPath?: string
): boolean {
  const appSourcesDir = join(appPath, 'sources');

  if (!existsSync(appSourcesDir)) {
    debug(`[copyAppSourcesToWorkspace] No sources directory in app`);
    return false;
  }

  const workspaceSourcesDir = join(workspaceRoot, '.sprouty-ai', 'sources');

  // Ensure workspace sources directory exists
  if (!existsSync(workspaceSourcesDir)) {
    mkdirSync(workspaceSourcesDir, { recursive: true });
  }

  // Iterate over each source subdirectory
  const sourceEntries = readdirSync(appSourcesDir, { withFileTypes: true });
  let copiedCount = 0;

  for (const entry of sourceEntries) {
    if (!entry.isDirectory()) continue;

    const sourceSlug = entry.name;
    const sourceDir = join(appSourcesDir, sourceSlug);
    const targetDir = join(workspaceSourcesDir, sourceSlug);

    // Skip if source already exists in workspace
    if (existsSync(targetDir)) {
      debug(`[copyAppSourcesToWorkspace] Source ${sourceSlug} already exists, skipping`);
      continue;
    }

    // Create target directory
    mkdirSync(targetDir, { recursive: true });

    // Copy files and resolve path placeholders in config.json
    const files = readdirSync(sourceDir, { withFileTypes: true });
    for (const file of files) {
      if (!file.isFile()) continue;

      const sourcePath = join(sourceDir, file.name);
      const targetPath = join(targetDir, file.name);

      if (file.name === 'config.json') {
        // Read, resolve paths, and write config.json
        try {
          const configContent = readFileSync(sourcePath, 'utf-8');
          const config = JSON.parse(configContent);

          // Resolve path placeholders
          const resolvedConfig = resolveSourceConfigPaths(config, isPackaged, resourcesPath);

          // Write resolved config
          writeFileSync(targetPath, JSON.stringify(resolvedConfig, null, 2), 'utf-8');
          debug(`[copyAppSourcesToWorkspace] Resolved paths in ${sourceSlug}/config.json`);
        } catch (error) {
          debug(`[copyAppSourcesToWorkspace] Failed to resolve paths in ${sourceSlug}/config.json:`, error);
          // Fall back to direct copy
          copyFileSync(sourcePath, targetPath);
        }
      } else {
        // Copy other files directly
        copyFileSync(sourcePath, targetPath);
      }
    }

    copiedCount++;
    debug(`[copyAppSourcesToWorkspace] Copied source ${sourceSlug} to workspace`);
  }

  debug(`[copyAppSourcesToWorkspace] Copied ${copiedCount} sources to workspace`);
  return copiedCount > 0;
}

// ============================================================
// App Manifest Copy
// ============================================================

/**
 * Copy app manifest.json to workspace root for app info display.
 * The manifest is copied to .sprouty-ai/app-manifest.json.
 */
function copyAppManifestToWorkspace(appId: string, workspaceRoot: string): void {
  const appPath = getAppPath(appId);

  if (!existsSync(appPath)) {
    debug(`[copyAppManifestToWorkspace] App path not found: ${appPath}`);
    return;
  }

  const manifestSource = join(appPath, 'manifest.json');
  if (!existsSync(manifestSource)) {
    debug(`[copyAppManifestToWorkspace] No manifest.json found at ${manifestSource}`);
    return;
  }

  const workspaceDataPath = getWorkspaceDataPath(workspaceRoot);
  const manifestTarget = join(workspaceDataPath, 'app-manifest.json');

  // Ensure workspace data directory exists
  if (!existsSync(workspaceDataPath)) {
    mkdirSync(workspaceDataPath, { recursive: true });
  }

  copyFileSync(manifestSource, manifestTarget);
  debug(`[copyAppManifestToWorkspace] Copied manifest to ${manifestTarget}`);
}

// ============================================================
// AGENTS.md Copy
// ============================================================

/**
 * Copy AGENTS.md from app source directory to workspace data directory.
 * This provides workspace-specific AI agent guidance based on the app type.
 *
 * The file is copied to .sprouty-ai/AGENTS.md in the workspace,
 * where it will be discovered by the agent's project context file mechanism.
 */
function copyAppAgentsToWorkspace(appId: string, workspaceRoot: string): void {
  const appPath = getAppPath(appId);
  const agentsPath = join(appPath, 'AGENTS.md');

  if (existsSync(agentsPath)) {
    copyAgentsFileToWorkspace(agentsPath, workspaceRoot);
  } else {
    debug(`[copyAppAgentsToWorkspace] No AGENTS.md found for app ${appId}`);
  }
}

/**
 * Copy AGENTS.md file to workspace data directory.
 */
function copyAgentsFileToWorkspace(sourcePath: string, workspaceRoot: string): void {
  const workspaceDataPath = getWorkspaceDataPath(workspaceRoot);
  const targetPath = join(workspaceDataPath, 'AGENTS.md');
  
  // Ensure workspace data directory exists
  if (!existsSync(workspaceDataPath)) {
    mkdirSync(workspaceDataPath, { recursive: true });
  }
  
  copyFileSync(sourcePath, targetPath);
  debug(`[copyAgentsFileToWorkspace] Copied AGENTS.md to ${targetPath}`);
}

// ============================================================
// Main Initialization Functions
// ============================================================

/**
 * Result of workspace initialization
 */
export interface InitializeWorkspaceResult {
  success: boolean;
  workspaceRoot: string;
  config: WorkspaceConfig;
  skillsInstalled: string[];
  skillsFailed: string[];
  errors: string[];
}

/**
 * Options for workspace initialization
 */
export interface InitializeWorkspaceOptions {
  /** Workspace display name */
  name: string;
  /** Root path for the workspace (if not provided, uses app default or generates one) */
  rootPath?: string;
  /** Application ID to use */
  appId: string;
  /** Override default settings */
  customSettings?: Record<string, any>;
  /** Skip skill installation */
  skipSkills?: boolean;
  /** Skip directory structure creation */
  skipDirectoryStructure?: boolean;
  /** Skip preset data */
  skipPresetData?: boolean;
  /** Whether the app is packaged (for path resolution) */
  isPackaged?: boolean;
  /** Resources path for packaged app (app.getPath('resources')) */
  resourcesPath?: string;
}

/**
 * Initialize a new workspace based on application configuration
 */
export function initializeWorkspaceFromApp(
  options: InitializeWorkspaceOptions
): InitializeWorkspaceResult {
  const errors: string[] = [];
  let skillsInstalled: string[] = [];
  let skillsFailed: string[] = [];
  
  // Load the application
  const app = loadAppById(options.appId);
  if (!app) {
    return {
      success: false,
      workspaceRoot: '',
      config: {} as WorkspaceConfig,
      skillsInstalled: [],
      skillsFailed: [],
      errors: [`Application not found: ${options.appId}`],
    };
  }
  
  const manifest = app.manifest;
  
  // Determine workspace root path
  let workspaceRoot = options.rootPath;
  if (!workspaceRoot) {
    // Use app default directory or generate one
    const defaultDir = manifest.workspace?.defaultDirectory || '~/SproutyAI';
    workspaceRoot = defaultDir.replace('~', process.env.HOME || '');
    workspaceRoot = join(workspaceRoot, options.name.replace(/[^a-zA-Z0-9-_]/g, '-'));
  }
  
  // Create workspace with app binding
  // Map app settings to workspace defaults format
  const appSettings = manifest.workspace?.defaultSettings;
  const workspaceDefaults: WorkspaceConfig['defaults'] = {
    // Model configuration
    model: appSettings?.defaultModel as string | undefined,
    workingDirectory: workspaceRoot,
    enabledSourceSlugs: appSettings?.enabledSourceSlugs ?? [],
    
    // Permission configuration from app manifest
    permissionMode: appSettings?.permissionMode,
    cyclablePermissionModes: appSettings?.cyclablePermissionModes,
    
    // Thinking level from app manifest
    thinkingLevel: appSettings?.thinkingLevel,
    
    // Custom settings override everything
    ...options.customSettings,
  };
  
  // Check if app has custom labels/statuses before creating workspace
  const hasCustomLabels = !options.skipPresetData && !!getAppPathWithCustomLabels(options.appId);
  const hasCustomStatuses = !options.skipPresetData && !!getAppPathWithCustomStatuses(options.appId);
  
  let config: WorkspaceConfig;
  try {
    config = createWorkspaceAtPath(workspaceRoot, options.name, {
      defaults: workspaceDefaults,
      skipDefaultLabels: hasCustomLabels,
      skipDefaultStatuses: hasCustomStatuses,
    });
    
    // Add app binding
    config.appId = options.appId;
    config.installedPluginApps = [];
    config.appSettings = {};
    
    // Apply localMcpServers configuration from app manifest
    if (manifest.workspace?.localMcpServers) {
      config.localMcpServers = manifest.workspace.localMcpServers;
    }
    
    saveWorkspaceConfig(workspaceRoot, config);
  } catch (error) {
    return {
      success: false,
      workspaceRoot,
      config: {} as WorkspaceConfig,
      skillsInstalled: [],
      skillsFailed: [],
      errors: [`Failed to create workspace: ${error}`],
    };
  }
  
  // Create directory structure
  if (!options.skipDirectoryStructure && manifest.workspace?.directoryStructure) {
    try {
      createDirectoryStructure(workspaceRoot, manifest.workspace.directoryStructure);
    } catch (error) {
      errors.push(`Failed to create directory structure: ${error}`);
    }
  }
  
  // Install bundled skills
  if (!options.skipSkills && manifest.capabilities?.skills) {
    const result = installBundledSkills(
      app.path,
      workspaceRoot,
      manifest.capabilities.skills
    );
    skillsInstalled = result.installed;
    skillsFailed = result.failed;
    
    if (result.failed.length > 0) {
      errors.push(`Failed to install skills: ${result.failed.join(', ')}`);
    }
  }
  
  // Copy app-specific labels and statuses (direct directory copy)
  if (!options.skipPresetData) {
    try {
      copyAppDataToWorkspace(
        options.appId,
        workspaceRoot,
        options.isPackaged || false,
        options.resourcesPath
      );
    } catch (error) {
      errors.push(`Failed to copy app data: ${error}`);
    }
  }
  
  // Copy AGENTS.md from app to workspace
  try {
    copyAppAgentsToWorkspace(options.appId, workspaceRoot);
  } catch (error) {
    errors.push(`Failed to copy AGENTS.md: ${error}`);
  }
  
  // Copy app manifest.json to workspace for app info display
  try {
    copyAppManifestToWorkspace(options.appId, workspaceRoot);
  } catch (error) {
    errors.push(`Failed to copy app manifest: ${error}`);
  }
  
  return {
    success: errors.length === 0,
    workspaceRoot,
    config,
    skillsInstalled,
    skillsFailed,
    errors,
  };
}

/**
 * Install a plugin to an existing workspace
 */
export function installPluginToWorkspace(
  workspaceRoot: string,
  pluginId: string
): { success: boolean; error?: string } {
  // Load workspace config
  const config = loadWorkspaceConfig(workspaceRoot);
  if (!config) {
    return { success: false, error: 'Workspace not found' };
  }
  
  // Load the plugin
  const plugin = loadAppById(pluginId);
  if (!plugin) {
    return { success: false, error: `Plugin not found: ${pluginId}` };
  }
  
  if (plugin.manifest.type !== 'plugin') {
    return { success: false, error: 'Not a plugin application' };
  }
  
  // Check compatibility with main app
  if (config.appId && plugin.manifest.compatibleApps) {
    const isCompatible = plugin.manifest.compatibleApps.some(compatible => {
      const [appId] = compatible.split('@');
      return appId === config.appId;
    });
    
    if (!isCompatible) {
      return { 
        success: false, 
        error: `Plugin ${pluginId} is not compatible with ${config.appId}` 
      };
    }
  }
  
  // Add plugin to installed list
  if (!config.installedPluginApps) {
    config.installedPluginApps = [];
  }
  
  if (!config.installedPluginApps.includes(pluginId)) {
    config.installedPluginApps.push(pluginId);
  }
  
  // Install plugin skills
  if (plugin.manifest.capabilities?.skills) {
    installBundledSkills(
      plugin.path,
      workspaceRoot,
      plugin.manifest.capabilities.skills
    );
  }
  
  // Save updated config
  saveWorkspaceConfig(workspaceRoot, config);
  
  return { success: true };
}

/**
 * Uninstall a plugin from a workspace
 */
export function uninstallPluginFromWorkspace(
  workspaceRoot: string,
  pluginId: string
): { success: boolean; error?: string } {
  const config = loadWorkspaceConfig(workspaceRoot);
  if (!config) {
    return { success: false, error: 'Workspace not found' };
  }
  
  if (!config.installedPluginApps) {
    return { success: false, error: 'No plugins installed' };
  }
  
  const index = config.installedPluginApps.indexOf(pluginId);
  if (index === -1) {
    return { success: false, error: 'Plugin not installed' };
  }
  
  // Remove from list
  config.installedPluginApps.splice(index, 1);
  
  // Note: We don't remove skills as they might have been modified by user
  // This is a design decision - skills are considered workspace-owned once installed
  
  saveWorkspaceConfig(workspaceRoot, config);
  
  return { success: true };
}

/**
 * Migrate an existing workspace to use an application
 * Used for backward compatibility with workspaces created before app framework
 */
export function migrateWorkspaceToApp(
  workspaceRoot: string,
  appId: string
): { success: boolean; error?: string } {
  const config = loadWorkspaceConfig(workspaceRoot);
  if (!config) {
    return { success: false, error: 'Workspace not found' };
  }
  
  if (config.appId) {
    return { success: false, error: 'Workspace already has an app binding' };
  }
  
  const app = loadAppById(appId);
  if (!app) {
    return { success: false, error: `Application not found: ${appId}` };
  }
  
  // Simply add the app binding
  config.appId = appId;
  config.installedPluginApps = [];
  config.appSettings = {};
  
  saveWorkspaceConfig(workspaceRoot, config);
  
  return { success: true };
}
