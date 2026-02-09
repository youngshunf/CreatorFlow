/**
 * Video Project Manager
 *
 * Manages video projects including CRUD operations, asset management,
 * and project persistence.
 *
 * Project Storage Structure:
 * {workspaceRoot}/
 * └── 视频创作/
 *     └── {项目名称}/
 *         ├── project.json          # 项目配置
 *         ├── 素材/                  # 素材文件
 *         │   ├── images/
 *         │   ├── videos/
 *         │   ├── audio/
 *         │   └── fonts/
 *         ├── 组合/                  # 组合代码
 *         │   └── {compositionId}.tsx
 *         └── 输出/                  # 渲染输出
 *             └── {renderName}.mp4
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
// Import MCP Server path utilities
import {
  getVideoProjectsDir,
  getProjectPath,
  getProjectConfigPath,
  getAssetsPath,
  getAssetTypePath,
  getCompositionsPath,
  getCompositionFilePath,
  getOutputPath,
  VIDEO_PROJECTS_DIR_NAME,
  ASSETS_DIR_NAME,
  COMPOSITIONS_DIR_NAME,
  OUTPUT_DIR_NAME,
  extractProjectName,
  getAssetRelativePath,
} from '@creator-flow/video/mcp-server/utils/paths';
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
  return getVideoProjectsDir(workspaceRootPath);
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

  /** Map of project ID to project name for quick lookup */
  private projectIdToNameCache: Map<string, string> = new Map();

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

    // Create project directory structure using MCP Server path utilities
    const projectPath = getProjectPath(workspaceRootPath, name);
    const compositionsPath = getCompositionsPath(workspaceRootPath, name);
    const assetsPath = getAssetsPath(workspaceRootPath, name);
    const outputPath = getOutputPath(workspaceRootPath, name);

    // Create all directories
    mkdirSync(projectPath, { recursive: true });
    mkdirSync(compositionsPath, { recursive: true });
    mkdirSync(getAssetTypePath(workspaceRootPath, name, 'image'), { recursive: true });
    mkdirSync(getAssetTypePath(workspaceRootPath, name, 'video'), { recursive: true });
    mkdirSync(getAssetTypePath(workspaceRootPath, name, 'audio'), { recursive: true });
    mkdirSync(getAssetTypePath(workspaceRootPath, name, 'font'), { recursive: true });
    mkdirSync(outputPath, { recursive: true });

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
          const compositionFilePath = getCompositionFilePath(workspaceRootPath, name, compositionId);

          // Write composition code to file
          writeFileSync(
            compositionFilePath,
            templateData.compositionCode,
            'utf-8'
          );

          compositions.push({
            id: compositionId,
            name: templateData.name,
            code: `${COMPOSITIONS_DIR_NAME}/${compositionId}.tsx`,
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

    // Update project-workspace mapping and ID-to-name cache
    this.projectWorkspaceMap.set(projectId, workspaceId);
    this.projectIdToNameCache.set(projectId, name);

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

        // Update project-workspace mapping and ID-to-name cache
        this.projectWorkspaceMap.set(validated.id, workspaceId);
        this.projectIdToNameCache.set(validated.id, validated.name);
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

    // Get workspace root path
    const workspaceId = this.projectWorkspaceMap.get(projectId);
    if (!workspaceId) {
      throw new Error(`Workspace not found for project: ${projectId}`);
    }
    const workspaceRootPath = getWorkspaceRootPath(workspaceId);

    // Generate asset ID and determine destination path
    const assetId = `asset_${randomUUID().replace(/-/g, '').slice(0, 8)}`;
    const fileName = basename(assetPath);
    const destDir = getAssetTypePath(workspaceRootPath, project.name, assetType);
    const destPath = join(destDir, `${assetId}_${fileName}`);

    // Ensure destination directory exists
    mkdirSync(destDir, { recursive: true });

    // Copy asset file to project
    copyFileSync(assetPath, destPath);

    // Create asset record with relative path
    const asset: Asset = {
      id: assetId,
      type: assetType,
      name: fileName,
      path: getAssetRelativePath(assetType, `${assetId}_${fileName}`),
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
    const projectName = this.projectIdToNameCache.get(projectId);

    if (workspaceId && projectName) {
      const workspaceRootPath = getWorkspaceRootPath(workspaceId);
      const projectPath = getProjectPath(workspaceRootPath, projectName);
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
      const videoProjectsDir = getVideoProjectsDir(workspaceRootPath);

      if (!existsSync(videoProjectsDir)) continue;

      // Scan all project directories
      const projectDirs = readdirSync(videoProjectsDir, { withFileTypes: true });
      for (const projectDir of projectDirs) {
        if (!projectDir.isDirectory()) continue;

        const projectPath = join(videoProjectsDir, projectDir.name);
        const projectJsonPath = join(projectPath, 'project.json');

        if (!existsSync(projectJsonPath)) continue;

        try {
          const data = readFileSync(projectJsonPath, 'utf-8');
          const parsed = JSON.parse(data);

          if (parsed.id === projectId) {
            // Update cache
            this.projectWorkspaceMap.set(projectId, ws.name);
            this.projectIdToNameCache.set(projectId, projectDir.name);
            return projectPath;
          }
        } catch (error) {
          // Skip invalid project files
          continue;
        }
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
}
