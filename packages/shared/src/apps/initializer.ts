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
import { join, basename } from 'path';
import type { AppManifest, DirectoryStructure, PresetData } from './types.ts';
import type { WorkspaceConfig } from '../workspaces/types.ts';
import {
  createWorkspaceAtPath,
  getWorkspaceSkillsPath,
  getWorkspaceDataPath,
  saveWorkspaceConfig,
  loadWorkspaceConfig,
} from '../workspaces/storage.ts';
import { loadAppById, getAppPath } from './storage.ts';
import { getBundledAppSourcePath } from './bundled-apps.ts';
import { saveStatusConfig, getDefaultStatusConfig } from '../statuses/storage.ts';
import { saveLabelConfig, getDefaultLabelConfig } from '../labels/storage.ts';
import { debug } from '../utils/debug.ts';

// ============================================================
// Directory Structure Creation
// ============================================================

/**
 * Recursively create directory structure from template
 */
function createDirectoryStructure(basePath: string, structure: DirectoryStructure): void {
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
    const skillSlug = skillRef.split('@')[0];
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

// ============================================================
// Preset Data Application
// ============================================================

/**
 * Apply preset labels to workspace
 */
function applyPresetLabels(
  workspaceRoot: string,
  labels: PresetData['labels']
): void {
  if (!labels || labels.length === 0) {
    return;
  }
  
  // Load current label config or create default
  const labelConfig = getDefaultLabelConfig();
  
  // Add preset labels to workspace
  // Note: Label config structure may vary, this is a simplified implementation
  // In practice, you'd want to properly integrate with the label system
  
  saveLabelConfig(workspaceRoot, labelConfig);
}

/**
 * Apply preset statuses to workspace
 */
function applyPresetStatuses(
  workspaceRoot: string,
  statuses: PresetData['statuses']
): void {
  if (!statuses || statuses.length === 0) {
    return;
  }
  
  // Load current status config or create default
  const statusConfig = getDefaultStatusConfig();
  
  // Preset statuses can be used to customize the default statuses
  // For now, we use defaults - can be extended later
  
  saveStatusConfig(workspaceRoot, statusConfig);
}

// ============================================================
// AGENTS.md Copy
// ============================================================

/**
 * Copy AGENTS.md from app source directory to workspace data directory.
 * This provides workspace-specific AI agent guidance based on the app type.
 * 
 * The file is copied to .creator-flow/AGENTS.md in the workspace,
 * where it will be discovered by the agent's project context file mechanism.
 */
function copyAppAgentsToWorkspace(appId: string, workspaceRoot: string): void {
  // Try to find the source path for the bundled app
  const sourcePath = getBundledAppSourcePath(appId);
  
  if (!sourcePath) {
    // For non-bundled apps, check the installed app directory
    const appPath = getAppPath(appId, false);
    const agentsPath = join(appPath, 'AGENTS.md');
    
    if (existsSync(agentsPath)) {
      copyAgentsFileToWorkspace(agentsPath, workspaceRoot);
    } else {
      debug(`[copyAppAgentsToWorkspace] No AGENTS.md found for app ${appId}`);
    }
    return;
  }
  
  // Copy from bundled app source directory
  const agentsSourcePath = join(sourcePath, 'AGENTS.md');
  if (existsSync(agentsSourcePath)) {
    copyAgentsFileToWorkspace(agentsSourcePath, workspaceRoot);
  } else {
    // Try the installed bundled app directory as fallback
    const bundledAppPath = getAppPath(appId, true);
    const bundledAgentsPath = join(bundledAppPath, 'AGENTS.md');
    
    if (existsSync(bundledAgentsPath)) {
      copyAgentsFileToWorkspace(bundledAgentsPath, workspaceRoot);
    } else {
      debug(`[copyAppAgentsToWorkspace] No AGENTS.md found for bundled app ${appId}`);
    }
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
    const defaultDir = manifest.workspace?.defaultDirectory || '~/CreatorFlow';
    workspaceRoot = defaultDir.replace('~', process.env.HOME || '');
    workspaceRoot = join(workspaceRoot, options.name.replace(/[^a-zA-Z0-9-_]/g, '-'));
  }
  
  // Create workspace with app binding
  // Map app settings to workspace defaults format
  const workspaceDefaults: WorkspaceConfig['defaults'] = {
    model: manifest.workspace?.defaultSettings?.defaultModel as string | undefined,
    workingDirectory: workspaceRoot,
    enabledSourceSlugs: [],
    ...options.customSettings,
  };
  
  let config: WorkspaceConfig;
  try {
    config = createWorkspaceAtPath(workspaceRoot, options.name, workspaceDefaults);
    
    // Add app binding
    config.appId = options.appId;
    config.installedPluginApps = [];
    config.appSettings = {};
    
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
  
  // Apply preset data
  if (!options.skipPresetData && manifest.workspace?.presetData) {
    const presetData = manifest.workspace.presetData;
    
    try {
      if (presetData.labels) {
        applyPresetLabels(workspaceRoot, presetData.labels);
      }
      if (presetData.statuses) {
        applyPresetStatuses(workspaceRoot, presetData.statuses);
      }
    } catch (error) {
      errors.push(`Failed to apply preset data: ${error}`);
    }
  }
  
  // Copy AGENTS.md from app to workspace
  try {
    copyAppAgentsToWorkspace(options.appId, workspaceRoot);
  } catch (error) {
    errors.push(`Failed to copy AGENTS.md: ${error}`);
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
