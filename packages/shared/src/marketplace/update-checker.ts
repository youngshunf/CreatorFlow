/**
 * Update Checker
 *
 * Handles version synchronization and update detection for installed skills.
 * Integrates with the cloud sync API to check for available updates.
 */

import type {
  InstalledItem,
  UpdateItem,
  SyncResponse,
  InstalledSkillInfo,
} from './types.ts';
import { syncInstalled } from './api.ts';
import {
  getInstalledSkills,
  readSkillMeta,
} from './storage.ts';
import { join } from 'path';

// ============================================================
// Update Detection
// ============================================================

/**
 * Check for updates for all installed skills in a workspace
 */
export async function checkForUpdates(
  workspaceRoot: string
): Promise<UpdateItem[]> {
  // Get installed skills
  const installed = getInstalledSkills(workspaceRoot);
  
  if (installed.length === 0) {
    return [];
  }

  // Build sync request
  const syncRequest: InstalledItem[] = installed
    .filter((skill) => skill.version !== 'unknown')
    .map((skill) => ({
      id: skill.skillId,
      version: skill.version.split('-local.')[0]!, // Use base version for comparison
      type: 'skill' as const,
    }));

  if (syncRequest.length === 0) {
    return [];
  }

  // Call sync API
  try {
    const response = await syncInstalled({ installed: syncRequest });
    return response.updates;
  } catch (error) {
    console.error('[UpdateChecker] Failed to check for updates:', error);
    return [];
  }
}

/**
 * Get installed skills with update status
 */
export async function getSkillsWithUpdateStatus(
  workspaceRoot: string
): Promise<InstalledSkillInfo[]> {
  const installed = getInstalledSkills(workspaceRoot);
  
  if (installed.length === 0) {
    return [];
  }

  // Check for updates
  const updates = await checkForUpdates(workspaceRoot);
  const updateMap = new Map(updates.map((u) => [u.id, u]));

  // Merge update info
  return installed.map((skill) => {
    const update = updateMap.get(skill.skillId);
    return {
      ...skill,
      hasUpdate: !!update,
      latestVersion: update?.latest_version,
    };
  });
}

/**
 * Check if a specific skill has an update available
 */
export async function hasSkillUpdate(
  workspaceRoot: string,
  skillId: string
): Promise<{ hasUpdate: boolean; latestVersion?: string; changelog?: string }> {
  const skillDir = join(workspaceRoot, '.sprouty-ai', 'skills', skillId);
  const meta = readSkillMeta(skillDir);
  
  if (!meta) {
    return { hasUpdate: false };
  }

  // Check with sync API
  try {
    const response = await syncInstalled({
      installed: [{
        id: skillId,
        version: meta.baseVersion,
        type: 'skill',
      }],
    });

    const update = response.updates.find((u) => u.id === skillId);
    if (update) {
      return {
        hasUpdate: true,
        latestVersion: update.latest_version,
        changelog: update.changelog || undefined,
      };
    }
  } catch (error) {
    console.error(`[UpdateChecker] Failed to check update for ${skillId}:`, error);
  }

  return { hasUpdate: false };
}

// ============================================================
// Background Sync
// ============================================================

let syncInterval: ReturnType<typeof setInterval> | null = null;
let lastSyncTime = 0;
const MIN_SYNC_INTERVAL = 5 * 60 * 1000; // 5 minutes minimum between syncs

/**
 * Callback type for update notifications
 */
export type UpdateCallback = (updates: UpdateItem[]) => void;

/**
 * Start background update checking
 */
export function startBackgroundSync(
  workspaceRoot: string,
  intervalMs: number = 30 * 60 * 1000, // 30 minutes default
  onUpdates?: UpdateCallback
): void {
  // Stop any existing sync
  stopBackgroundSync();

  // Initial sync
  performSync(workspaceRoot, onUpdates);

  // Schedule periodic sync
  syncInterval = setInterval(() => {
    performSync(workspaceRoot, onUpdates);
  }, intervalMs);
}

/**
 * Stop background update checking
 */
export function stopBackgroundSync(): void {
  if (syncInterval) {
    clearInterval(syncInterval);
    syncInterval = null;
  }
}

/**
 * Perform a sync operation
 */
async function performSync(
  workspaceRoot: string,
  onUpdates?: UpdateCallback
): Promise<void> {
  // Rate limit
  const now = Date.now();
  if (now - lastSyncTime < MIN_SYNC_INTERVAL) {
    return;
  }
  lastSyncTime = now;

  try {
    const updates = await checkForUpdates(workspaceRoot);
    if (updates.length > 0 && onUpdates) {
      onUpdates(updates);
    }
  } catch (error) {
    console.error('[UpdateChecker] Background sync failed:', error);
  }
}

/**
 * Force an immediate sync (respects rate limit)
 */
export async function forceSync(
  workspaceRoot: string,
  onUpdates?: UpdateCallback
): Promise<UpdateItem[]> {
  lastSyncTime = 0; // Reset rate limit
  await performSync(workspaceRoot, onUpdates);
  return checkForUpdates(workspaceRoot);
}

// ============================================================
// Version Comparison Utilities
// ============================================================

/**
 * Compare two semantic versions
 * Returns: -1 if a < b, 0 if a === b, 1 if a > b
 */
export function compareVersions(a: string, b: string): number {
  // Remove any local suffix
  const cleanA = a.split('-local.')[0]!;
  const cleanB = b.split('-local.')[0]!;

  const partsA = cleanA.split('.').map(Number);
  const partsB = cleanB.split('.').map(Number);

  for (let i = 0; i < Math.max(partsA.length, partsB.length); i++) {
    const numA = partsA[i] || 0;
    const numB = partsB[i] || 0;

    if (numA < numB) return -1;
    if (numA > numB) return 1;
  }

  return 0;
}

/**
 * Check if version a is older than version b
 */
export function isVersionOlder(a: string, b: string): boolean {
  return compareVersions(a, b) < 0;
}

/**
 * Check if a version is a local modification
 */
export function isLocalVersion(version: string): boolean {
  return version.includes('-local.');
}
