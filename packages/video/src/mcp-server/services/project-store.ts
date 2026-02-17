/**
 * 项目存储服务
 *
 * 负责视频项目的 CRUD 操作和文件系统持久化
 * 项目数据存储在用户可见的目录结构中，方便直接查看和编辑
 *
 * 文件系统结构:
 * {Workspace_Root}/
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
 */

import {
  existsSync,
  mkdirSync,
  readdirSync,
  readFileSync,
  rmSync,
  writeFileSync,
  statSync,
  copyFileSync,
  unlinkSync,
} from "fs";
import { join, basename } from "path";
import { randomUUID } from "crypto";
import {
  type VideoProject,
  type VideoConfig,
  type Scene,
  type Transition,
  type Asset,
  type AssetType,
  type CreateProjectInput,
  type UpdateProjectInput,
  type ProjectSummary,
  VideoProjectSchema,
} from "../types";
import {
  createProjectNotFoundError,
  createAssetNotFoundError,
  createSceneNotFoundError,
  createValidationError,
  MCPError,
} from "../types/errors";
import {
  getVideoProjectsDir,
  getProjectPath,
  getProjectConfigPath,
  getAssetsPath,
  getAssetTypePath,
  getCompositionsPath,
  getOutputPath,
  getAssetRelativePath,
  resolveVideoProjectPath,
  ASSET_SUBDIRS,
} from "../utils/paths";
import {
  assertValidProjectName,
  assertValidWorkspacePath,
  validateAssetFile,
  validateAssetType,
} from "../utils/validators";

// ============================================================================
// 类型定义
// ============================================================================

/**
 * ProjectStore 配置
 */
export interface ProjectStoreConfig {
  /** 工作区根路径 */
  workspacePath: string;
}

/**
 * 添加素材输入
 */
export interface AddAssetInput {
  /** 源文件路径（绝对路径） */
  sourcePath: string;
  /** 素材类型 */
  assetType: AssetType;
  /** 素材名称（可选，默认使用文件名） */
  name?: string;
}

/**
 * 按类型分组的素材列表
 */
export interface AssetsByType {
  image: Asset[];
  video: Asset[];
  audio: Asset[];
  font: Asset[];
}

/**
 * 添加场景输入
 */
export interface AddSceneInput {
  /** 场景名称 */
  name: string;
  /** 内置组合 ID */
  compositionId: string;
  /** 时长（帧数） */
  durationInFrames: number;
  /** 组合参数 */
  props?: Record<string, any>;
  /** 插入位置（默认末尾） */
  insertAt?: number;
}

/**
 * 更新场景输入
 */
export interface UpdateSceneInput {
  /** 组合参数 */
  props?: Record<string, any>;
  /** 时长（帧数） */
  durationInFrames?: number;
  /** 更换组合类型 */
  compositionId?: string;
  /** 场景名称 */
  name?: string;
}

// ============================================================================
// ProjectStore 类
// ============================================================================

/**
 * 项目存储服务
 *
 * 提供视频项目的 CRUD 操作，数据持久化到文件系统
 */
export class ProjectStore {
  private readonly workspacePath: string;

  /**
   * 创建 ProjectStore 实例
   * @param config 配置对象
   */
  constructor(config: ProjectStoreConfig) {
    this.workspacePath = config.workspacePath;
  }

  // ==========================================================================
  // 项目 CRUD 操作
  // ==========================================================================

  /**
   * 创建新项目
   * @param input 创建项目输入
   * @returns 创建的项目
   * @throws MCPError 如果项目名称无效或项目已存在
   */
  async createProject(input: CreateProjectInput): Promise<VideoProject> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 验证项目名称
    assertValidProjectName(input.name);

    // 解析项目路径（根据 creator.db 是否存在选择路径模式）
    const projectPath = resolveVideoProjectPath(this.workspacePath, input.name);
    if (existsSync(projectPath)) {
      throw createValidationError(
        `项目已存在: ${input.name}`,
        "name",
        "不存在的项目名称",
        input.name,
      );
    }

    // 生成项目 ID 和时间戳
    const now = new Date().toISOString();
    const projectId = randomUUID();

    // 计算帧数
    const fps = input.fps ?? 30;
    const durationInFrames = Math.ceil(input.durationInSeconds * fps);

    // 创建项目配置
    const config: VideoConfig = {
      width: input.width ?? 1920,
      height: input.height ?? 1080,
      fps,
      durationInFrames,
    };

    // 创建项目对象
    const project: VideoProject = {
      id: projectId,
      name: input.name,
      description: input.description,
      createdAt: now,
      updatedAt: now,
      config,
      scenes: [],
      transitions: [],
      assets: [],
      renders: [],
    };

    // 创建目录结构
    this.createProjectDirectories(input.name);

    // 保存项目配置
    this.saveProject(project);

    return project;
  }

  /**
   * 获取项目详情
   * @param projectId 项目 ID
   * @returns 项目详情，如果不存在则返回 null
   */
  async getProject(projectId: string): Promise<VideoProject | null> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 遍历所有项目查找匹配的 ID
    const projectsDir = getVideoProjectsDir(this.workspacePath);
    if (!existsSync(projectsDir)) {
      return null;
    }

    const projectDirs = readdirSync(projectsDir, { withFileTypes: true })
      .filter((dirent) => dirent.isDirectory())
      .map((dirent) => dirent.name);

    for (const projectName of projectDirs) {
      const project = this.loadProject(projectName);
      if (project && project.id === projectId) {
        return project;
      }
    }

    return null;
  }

  /**
   * 通过名称获取项目
   * @param projectName 项目名称
   * @returns 项目详情，如果不存在则返回 null
   */
  async getProjectByName(projectName: string): Promise<VideoProject | null> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    return this.loadProject(projectName);
  }

  /**
   * 列出所有项目
   * @returns 项目摘要列表，按更新时间降序排列
   */
  async listProjects(): Promise<ProjectSummary[]> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    const projectsDir = getVideoProjectsDir(this.workspacePath);
    if (!existsSync(projectsDir)) {
      return [];
    }

    const projectDirs = readdirSync(projectsDir, { withFileTypes: true })
      .filter((dirent) => dirent.isDirectory())
      .map((dirent) => dirent.name);

    const summaries: ProjectSummary[] = [];

    for (const projectName of projectDirs) {
      const project = this.loadProject(projectName);
      if (project) {
        summaries.push({
          id: project.id,
          name: project.name,
          description: project.description,
          createdAt: project.createdAt,
          updatedAt: project.updatedAt,
          sceneCount: project.scenes.length,
          assetCount: project.assets.length,
        });
      }
    }

    // 按更新时间降序排列
    summaries.sort((a, b) => {
      return new Date(b.updatedAt).getTime() - new Date(a.updatedAt).getTime();
    });

    return summaries;
  }

  /**
   * 更新项目
   * @param projectId 项目 ID
   * @param updates 更新内容
   * @returns 更新后的项目
   * @throws MCPError 如果项目不存在
   */
  async updateProject(
    projectId: string,
    updates: UpdateProjectInput,
  ): Promise<VideoProject> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取现有项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    const oldName = project.name;
    let needsRename = false;

    // 如果更新名称，需要验证新名称
    if (updates.name && updates.name !== project.name) {
      assertValidProjectName(updates.name);

      // 检查新名称是否已被使用
      const newProjectPath = getProjectPath(this.workspacePath, updates.name);
      if (existsSync(newProjectPath)) {
        throw createValidationError(
          `项目名称已存在: ${updates.name}`,
          "name",
          "不存在的项目名称",
          updates.name,
        );
      }

      needsRename = true;
    }

    // 更新项目字段
    const updatedProject: VideoProject = {
      ...project,
      name: updates.name ?? project.name,
      description:
        updates.description !== undefined
          ? updates.description
          : project.description,
      config: updates.config
        ? {
            ...project.config,
            ...updates.config,
          }
        : project.config,
      updatedAt: new Date().toISOString(),
    };

    // 如果需要重命名，先重命名目录
    if (needsRename && updates.name) {
      const oldPath = getProjectPath(this.workspacePath, oldName);
      const newPath = getProjectPath(this.workspacePath, updates.name);

      // 使用 fs.renameSync 重命名目录
      const { renameSync } = await import("fs");
      renameSync(oldPath, newPath);
    }

    // 保存更新后的项目
    this.saveProject(updatedProject);

    return updatedProject;
  }

  /**
   * 删除项目
   * @param projectId 项目 ID
   * @returns 是否删除成功
   * @throws MCPError 如果项目不存在
   */
  async deleteProject(projectId: string): Promise<boolean> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    // 删除项目目录
    const projectPath = getProjectPath(this.workspacePath, project.name);
    rmSync(projectPath, { recursive: true, force: true });

    return true;
  }

  // ==========================================================================
  // 素材管理操作
  // ==========================================================================

  /**
   * 添加素材到项目
   *
   * 将源文件复制到项目的素材目录，并在项目配置中注册
   *
   * @param projectId 项目 ID
   * @param input 添加素材输入
   * @returns 添加的素材对象
   * @throws MCPError 如果项目不存在、文件不存在或格式不支持
   */
  async addAsset(projectId: string, input: AddAssetInput): Promise<Asset> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    // 验证素材类型
    const typeResult = validateAssetType(input.assetType);
    if (!typeResult.valid) {
      throw typeResult.error;
    }

    // 验证素材文件（存在性和扩展名）
    const fileResult = validateAssetFile(input.sourcePath, input.assetType);
    if (!fileResult.valid) {
      throw fileResult.error;
    }

    // 生成素材 ID
    const assetId = randomUUID();

    // 获取文件名
    const fileName = basename(input.sourcePath);

    // 生成唯一的目标文件名（避免冲突）
    const targetFileName = this.generateUniqueFileName(
      project.name,
      input.assetType,
      fileName,
    );

    // 获取目标路径
    const targetDir = getAssetTypePath(
      this.workspacePath,
      project.name,
      input.assetType,
    );
    const targetPath = join(targetDir, targetFileName);

    // 确保目标目录存在
    if (!existsSync(targetDir)) {
      mkdirSync(targetDir, { recursive: true });
    }

    // 复制文件到项目素材目录
    copyFileSync(input.sourcePath, targetPath);

    // 创建素材对象
    const asset: Asset = {
      id: assetId,
      type: input.assetType,
      name: input.name ?? fileName,
      path: getAssetRelativePath(input.assetType, targetFileName),
    };

    // 更新项目配置
    const updatedProject: VideoProject = {
      ...project,
      assets: [...project.assets, asset],
      updatedAt: new Date().toISOString(),
    };

    // 保存项目
    this.saveProject(updatedProject);

    return asset;
  }

  /**
   * 从项目移除素材
   *
   * 从项目配置中移除素材注册，可选删除实际文件
   *
   * @param projectId 项目 ID
   * @param assetId 素材 ID
   * @param deleteFile 是否删除实际文件（默认 false）
   * @returns 是否移除成功
   * @throws MCPError 如果项目或素材不存在
   */
  async removeAsset(
    projectId: string,
    assetId: string,
    deleteFile: boolean = false,
  ): Promise<boolean> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    // 查找素材
    const assetIndex = project.assets.findIndex((a) => a.id === assetId);
    if (assetIndex === -1) {
      throw createAssetNotFoundError(assetId);
    }

    const asset = project.assets[assetIndex];
    if (!asset) {
      throw createAssetNotFoundError(assetId);
    }

    // 如果需要删除文件
    if (deleteFile) {
      const projectPath = getProjectPath(this.workspacePath, project.name);
      const assetFilePath = join(projectPath, asset.path);

      if (existsSync(assetFilePath)) {
        unlinkSync(assetFilePath);
      }
    }

    // 从项目中移除素材
    const updatedAssets = [...project.assets];
    updatedAssets.splice(assetIndex, 1);

    const updatedProject: VideoProject = {
      ...project,
      assets: updatedAssets,
      updatedAt: new Date().toISOString(),
    };

    // 保存项目
    this.saveProject(updatedProject);

    return true;
  }

  /**
   * 列出项目中的所有素材（按类型分组）
   *
   * @param projectId 项目 ID
   * @returns 按类型分组的素材列表
   * @throws MCPError 如果项目不存在
   */
  async listAssets(projectId: string): Promise<AssetsByType> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    // 按类型分组
    const assetsByType: AssetsByType = {
      image: [],
      video: [],
      audio: [],
      font: [],
    };

    for (const asset of project.assets) {
      assetsByType[asset.type].push(asset);
    }

    return assetsByType;
  }

  /**
   * 获取单个素材详情
   *
   * @param projectId 项目 ID
   * @param assetId 素材 ID
   * @returns 素材对象，如果不存在则返回 null
   * @throws MCPError 如果项目不存在
   */
  async getAsset(projectId: string, assetId: string): Promise<Asset | null> {
    // 验证工作区路径
    assertValidWorkspacePath(this.workspacePath);

    // 获取项目
    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    // 查找素材
    const asset = project.assets.find((a) => a.id === assetId);
    return asset ?? null;
  }

  // ==========================================================================
  // 场景管理操作
  // ==========================================================================

  /**
   * 添加场景到项目
   */
  async addScene(projectId: string, input: AddSceneInput): Promise<Scene> {
    assertValidWorkspacePath(this.workspacePath);

    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    if (!input.name || input.name.trim().length === 0) {
      throw createValidationError(
        "场景名称不能为空",
        "name",
        "非空字符串",
        input.name ?? "",
      );
    }
    if (!input.compositionId || input.compositionId.trim().length === 0) {
      throw createValidationError(
        "组合 ID 不能为空",
        "compositionId",
        "非空字符串",
        input.compositionId ?? "",
      );
    }
    if (!input.durationInFrames || input.durationInFrames <= 0) {
      throw createValidationError(
        "时长必须大于 0",
        "durationInFrames",
        "正整数",
        String(input.durationInFrames),
      );
    }

    const scene: Scene = {
      id: randomUUID(),
      name: input.name.trim(),
      compositionId: input.compositionId,
      durationInFrames: input.durationInFrames,
      props: input.props ?? {},
    };

    const updatedScenes = [...project.scenes];
    if (
      input.insertAt !== undefined &&
      input.insertAt >= 0 &&
      input.insertAt <= updatedScenes.length
    ) {
      updatedScenes.splice(input.insertAt, 0, scene);
    } else {
      updatedScenes.push(scene);
    }

    const updatedProject: VideoProject = {
      ...project,
      scenes: updatedScenes,
      updatedAt: new Date().toISOString(),
    };

    this.saveProject(updatedProject);
    return scene;
  }

  /**
   * 更新场景参数
   */
  async updateScene(
    projectId: string,
    sceneId: string,
    updates: UpdateSceneInput,
  ): Promise<Scene> {
    assertValidWorkspacePath(this.workspacePath);

    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    const sceneIndex = project.scenes.findIndex((s) => s.id === sceneId);
    if (sceneIndex === -1) {
      throw createSceneNotFoundError(sceneId);
    }

    const existing = project.scenes[sceneIndex]!;
    const updatedScene: Scene = {
      ...existing,
      ...(updates.name !== undefined && { name: updates.name.trim() }),
      ...(updates.compositionId !== undefined && {
        compositionId: updates.compositionId,
      }),
      ...(updates.durationInFrames !== undefined && {
        durationInFrames: updates.durationInFrames,
      }),
      ...(updates.props !== undefined && { props: updates.props }),
    };

    const updatedScenes = [...project.scenes];
    updatedScenes[sceneIndex] = updatedScene;

    const updatedProject: VideoProject = {
      ...project,
      scenes: updatedScenes,
      updatedAt: new Date().toISOString(),
    };

    this.saveProject(updatedProject);
    return updatedScene;
  }

  /**
   * 移除场景
   */
  async removeScene(projectId: string, sceneId: string): Promise<boolean> {
    assertValidWorkspacePath(this.workspacePath);

    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    const sceneIndex = project.scenes.findIndex((s) => s.id === sceneId);
    if (sceneIndex === -1) {
      throw createSceneNotFoundError(sceneId);
    }

    const updatedScenes = [...project.scenes];
    updatedScenes.splice(sceneIndex, 1);

    // 同步调整 transitions（长度 = scenes.length - 1）
    const updatedTransitions = [...project.transitions];
    if (updatedTransitions.length > 0) {
      if (sceneIndex === 0) {
        updatedTransitions.splice(0, 1);
      } else if (sceneIndex >= updatedScenes.length) {
        updatedTransitions.splice(updatedTransitions.length - 1, 1);
      } else {
        updatedTransitions.splice(sceneIndex - 1, 1);
      }
    }

    const updatedProject: VideoProject = {
      ...project,
      scenes: updatedScenes,
      transitions: updatedTransitions,
      updatedAt: new Date().toISOString(),
    };

    this.saveProject(updatedProject);
    return true;
  }

  /**
   * 重排场景顺序
   */
  async reorderScenes(
    projectId: string,
    sceneIds: string[],
  ): Promise<VideoProject> {
    assertValidWorkspacePath(this.workspacePath);

    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    if (sceneIds.length !== project.scenes.length) {
      throw createValidationError(
        `场景 ID 数量不匹配: 期望 ${project.scenes.length}，收到 ${sceneIds.length}`,
        "sceneIds",
        `${project.scenes.length} 个场景 ID`,
        `${sceneIds.length} 个场景 ID`,
      );
    }

    const sceneMap = new Map(project.scenes.map((s) => [s.id, s]));
    const reorderedScenes: Scene[] = [];

    for (const id of sceneIds) {
      const scene = sceneMap.get(id);
      if (!scene) {
        throw createSceneNotFoundError(id);
      }
      reorderedScenes.push(scene);
    }

    const updatedProject: VideoProject = {
      ...project,
      scenes: reorderedScenes,
      updatedAt: new Date().toISOString(),
    };

    this.saveProject(updatedProject);
    return updatedProject;
  }

  /**
   * 设置场景间过渡效果
   */
  async setTransitions(
    projectId: string,
    transitions: Transition[],
  ): Promise<VideoProject> {
    assertValidWorkspacePath(this.workspacePath);

    const project = await this.getProject(projectId);
    if (!project) {
      throw createProjectNotFoundError(projectId);
    }

    const expectedLength = Math.max(0, project.scenes.length - 1);
    if (transitions.length !== expectedLength) {
      throw createValidationError(
        `过渡效果数量不匹配: 期望 ${expectedLength}（场景数 - 1），收到 ${transitions.length}`,
        "transitions",
        `${expectedLength} 个过渡效果`,
        `${transitions.length} 个过渡效果`,
      );
    }

    const updatedProject: VideoProject = {
      ...project,
      transitions,
      updatedAt: new Date().toISOString(),
    };

    this.saveProject(updatedProject);
    return updatedProject;
  }

  // ==========================================================================
  // 路径工具方法
  // ==========================================================================

  /**
   * 获取项目目录路径
   * @param projectName 项目名称
   * @returns 项目目录路径
   */
  getProjectPath(projectName: string): string {
    return getProjectPath(this.workspacePath, projectName);
  }

  /**
   * 获取素材目录路径
   * @param projectName 项目名称
   * @returns 素材目录路径
   */
  getAssetsPath(projectName: string): string {
    return getAssetsPath(this.workspacePath, projectName);
  }

  /**
   * 获取输出目录路径
   * @param projectName 项目名称
   * @returns 输出目录路径
   */
  getOutputPath(projectName: string): string {
    return getOutputPath(this.workspacePath, projectName);
  }

  /**
   * 获取组合目录路径
   * @param projectName 项目名称
   * @returns 组合目录路径
   */
  getCompositionsPath(projectName: string): string {
    return getCompositionsPath(this.workspacePath, projectName);
  }

  // ==========================================================================
  // 私有方法
  // ==========================================================================

  /**
   * 生成唯一的文件名
   *
   * 如果目标目录中已存在同名文件，则添加数字后缀
   *
   * @param projectName 项目名称
   * @param assetType 素材类型
   * @param originalFileName 原始文件名
   * @returns 唯一的文件名
   */
  private generateUniqueFileName(
    projectName: string,
    assetType: AssetType,
    originalFileName: string,
  ): string {
    const targetDir = getAssetTypePath(
      this.workspacePath,
      projectName,
      assetType,
    );

    // 如果目标目录不存在，直接返回原始文件名
    if (!existsSync(targetDir)) {
      return originalFileName;
    }

    // 检查文件是否已存在
    let targetPath = join(targetDir, originalFileName);
    if (!existsSync(targetPath)) {
      return originalFileName;
    }

    // 分离文件名和扩展名
    const lastDotIndex = originalFileName.lastIndexOf(".");
    const nameWithoutExt =
      lastDotIndex > 0
        ? originalFileName.slice(0, lastDotIndex)
        : originalFileName;
    const ext = lastDotIndex > 0 ? originalFileName.slice(lastDotIndex) : "";

    // 添加数字后缀直到找到唯一的文件名
    let counter = 1;
    let newFileName = `${nameWithoutExt}_${counter}${ext}`;
    targetPath = join(targetDir, newFileName);

    while (existsSync(targetPath)) {
      counter++;
      newFileName = `${nameWithoutExt}_${counter}${ext}`;
      targetPath = join(targetDir, newFileName);
    }

    return newFileName;
  }

  /**
   * 创建项目目录结构
   * @param projectName 项目名称
   */
  private createProjectDirectories(projectName: string): void {
    // 创建视频项目根目录
    const projectsDir = getVideoProjectsDir(this.workspacePath);
    if (!existsSync(projectsDir)) {
      mkdirSync(projectsDir, { recursive: true });
    }

    // 创建项目目录
    const projectPath = getProjectPath(this.workspacePath, projectName);
    mkdirSync(projectPath, { recursive: true });

    // 创建素材目录及子目录
    const assetsPath = getAssetsPath(this.workspacePath, projectName);
    mkdirSync(assetsPath, { recursive: true });

    for (const subdir of Object.values(ASSET_SUBDIRS)) {
      mkdirSync(join(assetsPath, subdir), { recursive: true });
    }

    // 创建组合目录
    const compositionsPath = getCompositionsPath(
      this.workspacePath,
      projectName,
    );
    mkdirSync(compositionsPath, { recursive: true });

    // 创建输出目录
    const outputPath = getOutputPath(this.workspacePath, projectName);
    mkdirSync(outputPath, { recursive: true });
  }

  /**
   * 加载项目配置
   * @param projectName 项目名称
   * @returns 项目对象，如果不存在或无效则返回 null
   */
  private loadProject(projectName: string): VideoProject | null {
    const configPath = getProjectConfigPath(this.workspacePath, projectName);

    if (!existsSync(configPath)) {
      return null;
    }

    try {
      const content = readFileSync(configPath, "utf-8");
      const data = JSON.parse(content);

      // 使用 Zod schema 验证数据
      const result = VideoProjectSchema.safeParse(data);
      if (!result.success) {
        console.error(
          `[ProjectStore] Invalid project data in ${configPath}:`,
          result.error,
        );
        return null;
      }

      return result.data;
    } catch (error) {
      console.error(
        `[ProjectStore] Failed to load project ${projectName}:`,
        error,
      );
      return null;
    }
  }

  /**
   * 保存项目配置
   * @param project 项目对象
   */
  private saveProject(project: VideoProject): void {
    const configPath = getProjectConfigPath(this.workspacePath, project.name);

    // 序列化为 JSON（格式化以便用户阅读）
    const content = JSON.stringify(project, null, 2);

    writeFileSync(configPath, content, "utf-8");
  }

  // ==========================================================================
  // 静态工厂方法
  // ==========================================================================

  /**
   * 创建 ProjectStore 实例
   * @param workspacePath 工作区路径
   * @returns ProjectStore 实例
   */
  static create(workspacePath: string): ProjectStore {
    return new ProjectStore({ workspacePath });
  }
}

// ============================================================================
// 序列化/反序列化辅助函数
// ============================================================================

/**
 * 将项目对象序列化为 JSON 字符串
 * @param project 项目对象
 * @returns JSON 字符串
 */
export function serializeProject(project: VideoProject): string {
  return JSON.stringify(project, null, 2);
}

/**
 * 从 JSON 字符串反序列化项目对象
 * @param json JSON 字符串
 * @returns 项目对象，如果无效则返回 null
 */
export function deserializeProject(json: string): VideoProject | null {
  try {
    const data = JSON.parse(json);
    const result = VideoProjectSchema.safeParse(data);
    if (!result.success) {
      console.error("[ProjectStore] Invalid project JSON:", result.error);
      return null;
    }
    return result.data;
  } catch (error) {
    console.error("[ProjectStore] Failed to parse project JSON:", error);
    return null;
  }
}

/**
 * 验证项目对象是否有效
 * @param project 项目对象
 * @returns 是否有效
 */
export function isValidProject(project: unknown): project is VideoProject {
  const result = VideoProjectSchema.safeParse(project);
  return result.success;
}
