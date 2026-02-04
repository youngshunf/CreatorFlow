/**
 * Video Project Manager
 *
 * Manages video projects including CRUD operations, asset management,
 * and project persistence.
 *
 * Project Storage Structure:
 * ~/.creator-flow/workspaces/{workspaceId}/
 * └── .creator-flow/
 *     └── video-projects/
 *         └── {projectId}/
 *             ├── project.json          # Project metadata
 *             ├── compositions/         # Composition code files
 *             │   └── {compositionId}.tsx
 *             └── assets/               # Asset files
 *                 ├── images/
 *                 ├── videos/
 *                 ├── audio/
 *                 └── fonts/
 *
 * @requirements 12.1, 12.2, 12.3, 12.4
 */

import { join, basename, extname } from 'path';
import { homedir } from 'os';
import { randomUUID } from 'crypto';
import {
  existsSync,
  mkdirSync,
  readFileSync,
  writeFileSync,
  rmSync,
  copyFileSync,
  readdirSync,
  statSync,
} from 'fs';
// Import types from video package
import type {
  VideoProject,
  VideoConfig,
  Composition,
  Asset,
  AssetType,
} from '@creator-flow/video';
import {
  VideoProjectSchema,
  SUPPORTED_ASSET_EXTENSIONS,
} from '@creator-flow/video';
import log from '../logger';

const videoLog = log.scope('video');

/**
 * Options for creating a new video project
 */
export interface CreateProjectOptions {
  /** Project name */
  name: string;
  /** Workspace ID where the project will be stored */
  workspaceId: string;
  /** Optional template ID to initialize from */
  template?: string;
  /** Optional video configuration overrides */
  config?: Partial<VideoConfig>;
  /** Optional project description */
  description?: string;
}

/**
 * Interface for Project Manager
 */
export interface IProjectManager {
  createProject(options: CreateProjectOptions): Promise<VideoProject>;
  listProjects(workspaceId: string): Promise<VideoProject[]>;
  getProject(projectId: string): Promise<VideoProject | null>;
  updateProject(projectId: string, updates: Partial<VideoProject>): Promise<VideoProject>;
  deleteProject(projectId: string): Promise<boolean>;
  addAsset(projectId: string, assetPath: string, assetType: AssetType): Promise<Asset>;
  removeAsset(projectId: string, assetId: string): Promise<boolean>;
  getProjectPath(projectId: string): string;
}

/**
 * Default video configuration
 */
const DEFAULT_CONFIG: VideoConfig = {
  width: 1920,
  height: 1080,
  fps: 30,
  durationInFrames: 300, // 10 seconds at 30fps
};

/**
 * Get the base path for video projects in a workspace
 */
function getVideoProjectsBasePath(workspaceRootPath: string): string {
  return join(workspaceRootPath, '.creator-flow', 'video-projects');
}

/**
 * Get workspace root path from workspace ID
 * Workspaces are stored in ~/.creator-flow/workspaces/{workspaceId}/
 */
function getWorkspaceRootPath(workspaceId: string): string {
  return join(homedir(), '.creator-flow', 'workspaces', workspaceId);
}

/**
 * Project Manager Implementation
 */
export class ProjectManager implements IProjectManager {
  /** Map of project ID to workspace ID for quick lookup */
  private projectWorkspaceMap: Map<string, string> = new Map();

  constructor() {
    videoLog.info('ProjectManager initialized');
  }

  /**
   * Create a new video project
   * @requirements 12.1, 12.2
   */
  async createProject(options: CreateProjectOptions): Promise<VideoProject> {
    const { name, workspaceId, template, config, description } = options;

    // Generate unique project ID
    const projectId = `proj_${randomUUID().replace(/-/g, '').slice(0, 12)}`;

    // Get workspace root path
    const workspaceRootPath = getWorkspaceRootPath(workspaceId);
    if (!existsSync(workspaceRootPath)) {
      throw new Error(`Workspace not found: ${workspaceId}`);
    }

    // Create project directory structure
    const projectPath = join(getVideoProjectsBasePath(workspaceRootPath), projectId);
    const compositionsPath = join(projectPath, 'compositions');
    const assetsPath = join(projectPath, 'assets');

    mkdirSync(projectPath, { recursive: true });
    mkdirSync(compositionsPath, { recursive: true });
    mkdirSync(join(assetsPath, 'images'), { recursive: true });
    mkdirSync(join(assetsPath, 'videos'), { recursive: true });
    mkdirSync(join(assetsPath, 'audio'), { recursive: true });
    mkdirSync(join(assetsPath, 'fonts'), { recursive: true });

    // Initialize project configuration
    let projectConfig: VideoConfig = { ...DEFAULT_CONFIG };
    let compositions: Composition[] = [];

    // Apply template if specified
    if (template) {
      try {
        const { getTemplateById } = await import('@creator-flow/video');
        const templateData = getTemplateById(template);
        if (templateData) {
          projectConfig = { ...templateData.defaultConfig };
          // Create initial composition from template
          const compositionId = `comp_${randomUUID().replace(/-/g, '').slice(0, 8)}`;
          const compositionFileName = `${compositionId}.tsx`;

          // Write composition code to file
          writeFileSync(
            join(compositionsPath, compositionFileName),
            templateData.compositionCode,
            'utf-8'
          );

          compositions.push({
            id: compositionId,
            name: templateData.name,
            code: `compositions/${compositionFileName}`,
            props: templateData.defaultProps,
          });

          videoLog.info(`Applied template "${template}" to project ${projectId}`);
        }
      } catch (error) {
        videoLog.warn(`Failed to apply template "${template}":`, error);
        // Continue without template
      }
    }

    // Apply config overrides
    if (config) {
      projectConfig = { ...projectConfig, ...config };
    }

    // Create project object
    const now = new Date().toISOString();
    const project: VideoProject = {
      id: projectId,
      name,
      description,
      createdAt: now,
      updatedAt: now,
      config: projectConfig,
      compositions,
      assets: [],
      renders: [],
    };

    // Validate project data
    const validated = VideoProjectSchema.parse(project);

    // Save project to disk
    this.saveProjectToDisk(projectPath, validated);

    // Update project-workspace mapping
    this.projectWorkspaceMap.set(projectId, workspaceId);

    videoLog.info(`Created video project "${name}" (${projectId}) in workspace ${workspaceId}`);

    return validated;
  }

  /**
   * List all video projects in a workspace
   * @requirements 12.4
   */
  async listProjects(workspaceId: string): Promise<VideoProject[]> {
    const workspaceRootPath = getWorkspaceRootPath(workspaceId);
    const projectsBasePath = getVideoProjectsBasePath(workspaceRootPath);

    if (!existsSync(projectsBasePath)) {
      return [];
    }

    const projects: VideoProject[] = [];
    const entries = readdirSync(projectsBasePath, { withFileTypes: true });

    for (const entry of entries) {
      if (!entry.isDirectory()) continue;

      const projectPath = join(projectsBasePath, entry.name);
      const projectJsonPath = join(projectPath, 'project.json');

      if (!existsSync(projectJsonPath)) {
        videoLog.warn(`Project directory ${entry.name} missing project.json, skipping`);
        continue;
      }

      try {
        const data = readFileSync(projectJsonPath, 'utf-8');
        const parsed = JSON.parse(data);
        const validated = VideoProjectSchema.parse(parsed);
        projects.push(validated);

        // Update project-workspace mapping
        this.projectWorkspaceMap.set(validated.id, workspaceId);
      } catch (error) {
        // Log error and skip corrupted project (Requirement 12.5)
        videoLog.error(`Failed to load project ${entry.name}:`, error);
        continue;
      }
    }

    videoLog.info(`Loaded ${projects.length} video projects from workspace ${workspaceId}`);
    return projects;
  }

  /**
   * Get a single video project by ID
   */
  async getProject(projectId: string): Promise<VideoProject | null> {
    const projectPath = this.findProjectPath(projectId);
    if (!projectPath) {
      return null;
    }

    const projectJsonPath = join(projectPath, 'project.json');
    if (!existsSync(projectJsonPath)) {
      return null;
    }

    try {
      const data = readFileSync(projectJsonPath, 'utf-8');
      const parsed = JSON.parse(data);
      const validated = VideoProjectSchema.parse(parsed);
      return validated;
    } catch (error) {
      videoLog.error(`Failed to load project ${projectId}:`, error);
      return null;
    }
  }

  /**
   * Update a video project
   */
  async updateProject(
    projectId: string,
    updates: Partial<VideoProject>
  ): Promise<VideoProject> {
    const project = await this.getProject(projectId);
    if (!project) {
      throw new Error(`Project not found: ${projectId}`);
    }

    // Merge updates (excluding id, createdAt which are immutable)
    const updatedProject: VideoProject = {
      ...project,
      ...updates,
      id: project.id, // Preserve original ID
      createdAt: project.createdAt, // Preserve creation time
      updatedAt: new Date().toISOString(),
    };

    // Validate updated project
    const validated = VideoProjectSchema.parse(updatedProject);

    // Save to disk
    const projectPath = this.findProjectPath(projectId);
    if (projectPath) {
      this.saveProjectToDisk(projectPath, validated);
    }

    videoLog.info(`Updated video project ${projectId}`);
    return validated;
  }

  /**
   * Delete a video project
   */
  async deleteProject(projectId: string): Promise<boolean> {
    const projectPath = this.findProjectPath(projectId);
    if (!projectPath || !existsSync(projectPath)) {
      return false;
    }

    try {
      rmSync(projectPath, { recursive: true, force: true });
      this.projectWorkspaceMap.delete(projectId);
      videoLog.info(`Deleted video project ${projectId}`);
      return true;
    } catch (error) {
      videoLog.error(`Failed to delete project ${projectId}:`, error);
      return false;
    }
  }

  /**
   * Add an asset to a project
   * @requirements 12.3
   */
  async addAsset(
    projectId: string,
    assetPath: string,
    assetType: AssetType
  ): Promise<Asset> {
    // Validate source file exists
    if (!existsSync(assetPath)) {
      throw new Error(`Asset file not found: ${assetPath}`);
    }

    // Validate file extension
    const ext = extname(assetPath).toLowerCase();
    const supportedExts = SUPPORTED_ASSET_EXTENSIONS[assetType];
    if (!supportedExts.includes(ext as any)) {
      throw new Error(
        `Unsupported ${assetType} format: ${ext}. Supported formats: ${supportedExts.join(', ')}`
      );
    }

    // Get project
    const project = await this.getProject(projectId);
    if (!project) {
      throw new Error(`Project not found: ${projectId}`);
    }

    // Get project path
    const projectPath = this.findProjectPath(projectId);
    if (!projectPath) {
      throw new Error(`Project path not found: ${projectId}`);
    }

    // Generate asset ID and determine destination path
    const assetId = `asset_${randomUUID().replace(/-/g, '').slice(0, 8)}`;
    const fileName = basename(assetPath);
    const assetSubDir = this.getAssetSubDir(assetType);
    const destDir = join(projectPath, 'assets', assetSubDir);
    const destPath = join(destDir, `${assetId}_${fileName}`);

    // Ensure destination directory exists
    mkdirSync(destDir, { recursive: true });

    // Copy asset file to project
    copyFileSync(assetPath, destPath);

    // Create asset record
    const asset: Asset = {
      id: assetId,
      type: assetType,
      name: fileName,
      path: `assets/${assetSubDir}/${assetId}_${fileName}`,
    };

    // Update project with new asset
    const updatedAssets = [...project.assets, asset];
    await this.updateProject(projectId, { assets: updatedAssets });

    videoLog.info(`Added ${assetType} asset "${fileName}" to project ${projectId}`);
    return asset;
  }

  /**
   * Remove an asset from a project
   */
  async removeAsset(projectId: string, assetId: string): Promise<boolean> {
    const project = await this.getProject(projectId);
    if (!project) {
      return false;
    }

    const asset = project.assets.find((a: Asset) => a.id === assetId);
    if (!asset) {
      return false;
    }

    // Get project path
    const projectPath = this.findProjectPath(projectId);
    if (!projectPath) {
      return false;
    }

    // Delete asset file
    const assetFilePath = join(projectPath, asset.path);
    if (existsSync(assetFilePath)) {
      try {
        rmSync(assetFilePath);
      } catch (error) {
        videoLog.warn(`Failed to delete asset file ${assetFilePath}:`, error);
      }
    }

    // Update project to remove asset
    const updatedAssets = project.assets.filter((a: Asset) => a.id !== assetId);
    await this.updateProject(projectId, { assets: updatedAssets });

    videoLog.info(`Removed asset ${assetId} from project ${projectId}`);
    return true;
  }

  /**
   * Get the file system path for a project
   */
  getProjectPath(projectId: string): string {
    const path = this.findProjectPath(projectId);
    if (!path) {
      throw new Error(`Project not found: ${projectId}`);
    }
    return path;
  }

  // ============================================================================
  // Private Helper Methods
  // ============================================================================

  /**
   * Find the file system path for a project
   */
  private findProjectPath(projectId: string): string | null {
    // Check cache first
    const workspaceId = this.projectWorkspaceMap.get(projectId);
    if (workspaceId) {
      const workspaceRootPath = getWorkspaceRootPath(workspaceId);
      const projectPath = join(getVideoProjectsBasePath(workspaceRootPath), projectId);
      if (existsSync(projectPath)) {
        return projectPath;
      }
    }

    // Search all workspaces
    const workspacesDir = join(homedir(), '.creator-flow', 'workspaces');
    if (!existsSync(workspacesDir)) {
      return null;
    }

    const workspaces = readdirSync(workspacesDir, { withFileTypes: true });
    for (const ws of workspaces) {
      if (!ws.isDirectory()) continue;

      const workspaceRootPath = join(workspacesDir, ws.name);
      const projectPath = join(getVideoProjectsBasePath(workspaceRootPath), projectId);

      if (existsSync(projectPath)) {
        // Update cache
        this.projectWorkspaceMap.set(projectId, ws.name);
        return projectPath;
      }
    }

    return null;
  }

  /**
   * Save project data to disk
   * @requirements 12.2
   */
  private saveProjectToDisk(projectPath: string, project: VideoProject): void {
    const projectJsonPath = join(projectPath, 'project.json');
    writeFileSync(projectJsonPath, JSON.stringify(project, null, 2), 'utf-8');
  }

  /**
   * Get the subdirectory name for an asset type
   */
  private getAssetSubDir(assetType: AssetType): string {
    switch (assetType) {
      case 'image':
        return 'images';
      case 'video':
        return 'videos';
      case 'audio':
        return 'audio';
      case 'font':
        return 'fonts';
      default:
        return 'other';
    }
  }
}
