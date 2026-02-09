/**
 * Application Storage
 *
 * CRUD operations for applications.
 * Applications are stored globally at ~/.sprouty-ai/apps/
 * Built-in (bundled) apps are in ~/.sprouty-ai/apps/bundled/
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

const CONFIG_DIR = join(homedir(), '.sprouty-ai');
const APPS_DIR = join(CONFIG_DIR, 'apps');

/**
 * Get the global apps directory
 */
export function getAppsDir(): string {
  return APPS_DIR;
}

/**
 * Ensure apps directories exist
 */
export function ensureAppsDir(): void {
  if (!existsSync(APPS_DIR)) {
    mkdirSync(APPS_DIR, { recursive: true });
  }
}

/**
 * Get path to an app directory
 */
export function getAppPath(appId: string): string {
  // Convert app.name format to directory name
  const dirName = appId.replace(/\./g, '-');
  return join(APPS_DIR, dirName);
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
export function loadApp(appPath: string): LoadedApp | null {
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
    bundled: false,
  };
}

/**
 * Load an app by its ID
 */
export function loadAppById(appId: string): LoadedApp | null {
  const appPath = getAppPath(appId);
  if (existsSync(appPath)) {
    return loadApp(appPath);
  }

  return null;
}

/**
 * List all apps from a directory
 */
function listAppsFromDir(dir: string): LoadedApp[] {
  if (!existsSync(dir)) {
    return [];
  }

  const apps: LoadedApp[] = [];

  try {
    const entries = readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const appPath = join(dir, entry.name);
      const app = loadApp(appPath);
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
 * List all installed apps
 */
export function listAllApps(): LoadedApp[] {
  ensureAppsDir();
  return listAppsFromDir(APPS_DIR);
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
export function saveAppManifest(appId: string, manifest: AppManifest): void {
  const appPath = getAppPath(appId);

  if (!existsSync(appPath)) {
    mkdirSync(appPath, { recursive: true });
  }

  const manifestPath = getManifestPath(appPath);
  writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
}

/**
 * Create app directory structure
 */
export function createAppDirectory(appId: string): string {
  const appPath = getAppPath(appId);

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
  const appPath = getAppPath(appId);

  if (!existsSync(appPath)) {
    return false;
  }

  try {
    rmSync(appPath, { recursive: true });
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
  return existsSync(getAppPath(appId));
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
