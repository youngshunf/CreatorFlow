/**
 * Application Storage
 *
 * CRUD operations for applications.
 * Applications are stored globally at ~/.creator-flow/apps/
 * Built-in (bundled) apps are in ~/.creator-flow/apps/bundled/
 */

import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  readdirSync,
  rmSync,
  copyFileSync,
  statSync,
} from 'fs';
import { join } from 'path';
import { homedir } from 'os';
import type {
  AppManifest,
  LoadedApp,
  AppSummary,
} from './types.ts';
import { findIconFile } from '../utils/icon.ts';

// ============================================================
// Constants and Paths
// ============================================================

const CONFIG_DIR = join(homedir(), '.creator-flow');
const APPS_DIR = join(CONFIG_DIR, 'apps');
const BUNDLED_APPS_DIR = join(APPS_DIR, 'bundled');

/**
 * Get the global apps directory
 */
export function getAppsDir(): string {
  return APPS_DIR;
}

/**
 * Get the bundled apps directory
 */
export function getBundledAppsDir(): string {
  return BUNDLED_APPS_DIR;
}

/**
 * Ensure apps directories exist
 */
export function ensureAppsDir(): void {
  if (!existsSync(APPS_DIR)) {
    mkdirSync(APPS_DIR, { recursive: true });
  }
  if (!existsSync(BUNDLED_APPS_DIR)) {
    mkdirSync(BUNDLED_APPS_DIR, { recursive: true });
  }
}

/**
 * Get path to an app directory
 */
export function getAppPath(appId: string, bundled = false): string {
  const baseDir = bundled ? BUNDLED_APPS_DIR : APPS_DIR;
  // Convert app.name format to directory name
  const dirName = appId.replace(/\./g, '-');
  return join(baseDir, dirName);
}

/**
 * Get path to app manifest.json
 */
function getManifestPath(appPath: string): string {
  return join(appPath, 'manifest.json');
}

// ============================================================
// Load Operations
// ============================================================

/**
 * Parse and validate app manifest
 */
function parseManifest(content: string): AppManifest | null {
  try {
    const manifest = JSON.parse(content) as AppManifest;
    
    // Validate required fields
    if (!manifest.id || !manifest.name || !manifest.version) {
      return null;
    }
    
    // Set default type
    if (!manifest.type) {
      manifest.type = manifest.id.startsWith('plugin.') ? 'plugin' : 'app';
    }
    
    return manifest;
  } catch {
    return null;
  }
}

/**
 * Load a single app from a directory
 */
export function loadApp(appPath: string, bundled = false): LoadedApp | null {
  const manifestPath = getManifestPath(appPath);
  
  if (!existsSync(manifestPath)) {
    return null;
  }
  
  let content: string;
  try {
    content = readFileSync(manifestPath, 'utf-8');
  } catch {
    return null;
  }
  
  const manifest = parseManifest(content);
  if (!manifest) {
    return null;
  }
  
  return {
    manifest,
    path: appPath,
    iconPath: findIconFile(appPath) || undefined,
    bundled,
  };
}

/**
 * Load an app by its ID
 */
export function loadAppById(appId: string): LoadedApp | null {
  // Try bundled first
  const bundledPath = getAppPath(appId, true);
  if (existsSync(bundledPath)) {
    return loadApp(bundledPath, true);
  }
  
  // Then try user-installed
  const userPath = getAppPath(appId, false);
  if (existsSync(userPath)) {
    return loadApp(userPath, false);
  }
  
  return null;
}

/**
 * List all apps from a directory
 */
function listAppsFromDir(dir: string, bundled: boolean): LoadedApp[] {
  if (!existsSync(dir)) {
    return [];
  }
  
  const apps: LoadedApp[] = [];
  
  try {
    const entries = readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;
      // Skip bundled subdirectory when listing from APPS_DIR
      if (entry.name === 'bundled') continue;
      
      const appPath = join(dir, entry.name);
      const app = loadApp(appPath, bundled);
      if (app) {
        apps.push(app);
      }
    }
  } catch {
    // Ignore errors
  }
  
  return apps;
}

/**
 * List all installed apps (bundled + user-installed)
 */
export function listAllApps(): LoadedApp[] {
  ensureAppsDir();
  
  const bundledApps = listAppsFromDir(BUNDLED_APPS_DIR, true);
  const userApps = listAppsFromDir(APPS_DIR, false);
  
  // Merge, user apps can override bundled apps with same ID
  const appMap = new Map<string, LoadedApp>();
  
  for (const app of bundledApps) {
    appMap.set(app.manifest.id, app);
  }
  
  for (const app of userApps) {
    appMap.set(app.manifest.id, app);
  }
  
  return Array.from(appMap.values());
}

/**
 * List main apps only (exclude plugins)
 */
export function listMainApps(): LoadedApp[] {
  return listAllApps().filter(app => app.manifest.type !== 'plugin');
}

/**
 * List plugin apps only
 */
export function listPluginApps(): LoadedApp[] {
  return listAllApps().filter(app => app.manifest.type === 'plugin');
}

/**
 * Get app summaries for UI listing (lightweight)
 */
export function getAppSummaries(): AppSummary[] {
  return listAllApps().map(app => ({
    id: app.manifest.id,
    name: app.manifest.name,
    description: app.manifest.description,
    version: app.manifest.version,
    type: app.manifest.type || 'app',
    iconPath: app.iconPath,
    bundled: app.bundled,
  }));
}

// ============================================================
// Save Operations
// ============================================================

/**
 * Save an app manifest to disk
 */
export function saveAppManifest(appId: string, manifest: AppManifest, bundled = false): void {
  const appPath = getAppPath(appId, bundled);
  
  if (!existsSync(appPath)) {
    mkdirSync(appPath, { recursive: true });
  }
  
  const manifestPath = getManifestPath(appPath);
  writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
}

/**
 * Create app directory structure
 */
export function createAppDirectory(appId: string, bundled = false): string {
  const appPath = getAppPath(appId, bundled);
  
  // Create main directory
  mkdirSync(appPath, { recursive: true });
  
  // Create subdirectories
  mkdirSync(join(appPath, 'assets'), { recursive: true });
  mkdirSync(join(appPath, 'skills'), { recursive: true });
  mkdirSync(join(appPath, 'config'), { recursive: true });
  
  return appPath;
}

// ============================================================
// Delete Operations
// ============================================================

/**
 * Delete an installed app
 */
export function deleteApp(appId: string): boolean {
  // Don't allow deleting bundled apps
  const userPath = getAppPath(appId, false);
  
  if (!existsSync(userPath)) {
    return false;
  }
  
  try {
    rmSync(userPath, { recursive: true });
    return true;
  } catch {
    return false;
  }
}

// ============================================================
// Utility Functions
// ============================================================

/**
 * Check if an app exists
 */
export function appExists(appId: string): boolean {
  return existsSync(getAppPath(appId, true)) || existsSync(getAppPath(appId, false));
}

/**
 * Check if an app is bundled
 */
export function isAppBundled(appId: string): boolean {
  return existsSync(getAppPath(appId, true));
}

/**
 * Get compatible plugins for a main app
 */
export function getCompatiblePlugins(mainAppId: string): LoadedApp[] {
  const plugins = listPluginApps();
  
  return plugins.filter(plugin => {
    const compatibleApps = plugin.manifest.compatibleApps || [];
    return compatibleApps.some(compatible => {
      const [appId] = compatible.split('@');
      return appId === mainAppId;
    });
  });
}

// ============================================================
// Built-in App Registration
// ============================================================

/**
 * Register a bundled app (called at startup)
 * This ensures bundled apps are available in the apps directory
 */
export function registerBundledApp(manifest: AppManifest): void {
  ensureAppsDir();
  saveAppManifest(manifest.id, manifest, true);
}
