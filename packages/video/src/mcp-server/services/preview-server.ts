/**
 * 预览服务器管理
 *
 * 管理 Remotion Studio 预览服务器的生命周期
 * 支持启动、停止预览服务器，以及跟踪活跃服务器
 *
 * @requirements 6.1, 6.3, 6.4, 6.5
 */

import { spawn, type ChildProcess } from 'child_process';
import { existsSync } from 'fs';
import { createServer, type Server } from 'net';
import { join, dirname } from 'path';
import type { PreviewInstance } from '../types';
import {
  createPreviewFailedError,
  createProjectNotFoundError,
  MCPError,
} from '../types/errors';

// ============================================================================
// 类型定义
// ============================================================================

/**
 * 预览配置
 */
export interface PreviewConfig {
  /** 项目路径（包含 Remotion 入口文件的目录） */
  projectPath: string;
  /** 项目 ID */
  projectId: string;
  /** 组合 ID（可选，用于指定默认显示的组合） */
  compositionId?: string;
  /** 端口号（可选，默认自动查找可用端口） */
  port?: number;
}

/**
 * 内部预览实例（包含进程引用）
 */
interface InternalPreviewInstance extends PreviewInstance {
  /** 子进程引用 */
  process: ChildProcess;
}

// ============================================================================
// 常量定义
// ============================================================================

/** 默认起始端口 */
const DEFAULT_START_PORT = 3100;

/** 最大端口号 */
const MAX_PORT = 65535;

/** 服务器启动超时时间（毫秒） */
const SERVER_START_TIMEOUT = 30000;

/** 服务器停止超时时间（毫秒） */
const SERVER_STOP_TIMEOUT = 5000;

// ============================================================================
// PreviewServerManager 类
// ============================================================================

/**
 * 预览服务器管理器
 *
 * 提供预览服务器的生命周期管理：
 * - 启动 Remotion Studio 预览服务器
 * - 停止预览服务器
 * - 跟踪活跃服务器
 * - 防止同一项目重复启动服务器
 */
export class PreviewServerManager {
  /** 活跃的预览服务器映射（projectId -> PreviewInstance） */
  private activeServers: Map<string, InternalPreviewInstance> = new Map();

  constructor() {
    // 注册进程退出时的清理
    this.setupCleanup();
  }

  // ==========================================================================
  // 公共方法
  // ==========================================================================

  /**
   * 启动预览服务器
   *
   * @param config 预览配置
   * @returns 预览实例
   * @throws MCPError 如果启动失败或项目路径无效
   *
   * @requirements 6.1 - 启动预览服务器并返回 URL
   * @requirements 6.5 - 防止同一项目重复启动服务器
   */
  async start(config: PreviewConfig): Promise<PreviewInstance> {
    const { projectPath, projectId, compositionId, port } = config;

    console.error(`[PreviewServerManager] Starting preview server for project: ${projectId}`);

    // 检查是否已有活跃的服务器
    // @requirements 6.5 - 防止重复启动
    const existingServer = this.activeServers.get(projectId);
    if (existingServer) {
      console.error(`[PreviewServerManager] Returning existing server for project: ${projectId}`);
      return this.toPublicInstance(existingServer);
    }

    // 验证项目路径
    if (!existsSync(projectPath)) {
      throw createProjectNotFoundError(projectId);
    }

    // 查找可用端口
    const serverPort = port ?? await this.findAvailablePort(DEFAULT_START_PORT);
    console.error(`[PreviewServerManager] Using port: ${serverPort}`);

    try {
      // 启动 Remotion Studio
      const instance = await this.startRemotionStudio(
        projectPath,
        projectId,
        serverPort,
        compositionId
      );

      // 保存到活跃服务器映射
      this.activeServers.set(projectId, instance);

      console.error(`[PreviewServerManager] Preview server started: ${instance.url}`);

      return this.toPublicInstance(instance);
    } catch (error) {
      // @requirements 6.4 - 返回启动失败的错误
      if (error instanceof MCPError) {
        throw error;
      }

      const message = error instanceof Error ? error.message : String(error);
      throw createPreviewFailedError(message, { path: projectPath });
    }
  }

  /**
   * 停止预览服务器
   *
   * @param projectId 项目 ID
   * @returns 是否成功停止
   *
   * @requirements 6.3 - 优雅关闭服务器并释放端口
   */
  async stop(projectId: string): Promise<boolean> {
    console.error(`[PreviewServerManager] Stopping preview server for project: ${projectId}`);

    const instance = this.activeServers.get(projectId);
    if (!instance) {
      console.error(`[PreviewServerManager] No active server found for project: ${projectId}`);
      return false;
    }

    try {
      // 优雅关闭进程
      await this.stopProcess(instance.process);

      // 从活跃服务器映射中移除
      this.activeServers.delete(projectId);

      console.error(`[PreviewServerManager] Preview server stopped for project: ${projectId}`);
      return true;
    } catch (error) {
      console.error(`[PreviewServerManager] Error stopping server:`, error);

      // 强制终止进程
      try {
        instance.process.kill('SIGKILL');
      } catch {
        // 忽略强制终止的错误
      }

      // 从映射中移除
      this.activeServers.delete(projectId);

      return true;
    }
  }

  /**
   * 获取活跃的预览服务器
   *
   * @param projectId 项目 ID
   * @returns 预览实例，如果不存在则返回 null
   */
  getActiveServer(projectId: string): PreviewInstance | null {
    const instance = this.activeServers.get(projectId);
    return instance ? this.toPublicInstance(instance) : null;
  }

  /**
   * 列出所有活跃服务器
   *
   * @returns 所有活跃的预览实例列表
   */
  listActiveServers(): PreviewInstance[] {
    return Array.from(this.activeServers.values()).map(instance =>
      this.toPublicInstance(instance)
    );
  }

  /**
   * 停止所有活跃服务器
   *
   * 用于清理资源
   */
  async stopAll(): Promise<void> {
    console.error(`[PreviewServerManager] Stopping all preview servers...`);

    const projectIds = Array.from(this.activeServers.keys());
    await Promise.all(projectIds.map(id => this.stop(id)));

    console.error(`[PreviewServerManager] All preview servers stopped`);
  }

  // ==========================================================================
  // 私有方法
  // ==========================================================================

  /**
   * 查找可用端口
   *
   * @param startPort 起始端口
   * @returns 可用的端口号
   * @throws MCPError 如果找不到可用端口
   */
  private async findAvailablePort(startPort: number): Promise<number> {
    let port = startPort;

    while (port <= MAX_PORT) {
      const isAvailable = await this.isPortAvailable(port);
      if (isAvailable) {
        return port;
      }
      port++;
    }

    throw createPreviewFailedError('无法找到可用端口');
  }

  /**
   * 检查端口是否可用
   *
   * @param port 端口号
   * @returns 端口是否可用
   */
  private isPortAvailable(port: number): Promise<boolean> {
    return new Promise((resolve) => {
      const server: Server = createServer();

      server.once('error', () => {
        resolve(false);
      });

      server.once('listening', () => {
        server.close(() => {
          resolve(true);
        });
      });

      server.listen(port, '127.0.0.1');
    });
  }

  /**
   * 启动 Remotion Studio
   *
   * @param projectPath 项目路径
   * @param projectId 项目 ID
   * @param port 端口号
   * @param compositionId 组合 ID（可选）
   * @returns 内部预览实例
   */
  private async startRemotionStudio(
    projectPath: string,
    projectId: string,
    port: number,
    compositionId?: string
  ): Promise<InternalPreviewInstance> {
    return new Promise((resolve, reject) => {
      // 查找入口文件
      const entryPoint = this.findEntryPoint(projectPath);
      console.error(`[PreviewServerManager] Using entry point: ${entryPoint}`);

      // 构建命令参数
      const args = [
        'remotion',
        'studio',
        entryPoint,
        '--port',
        port.toString(),
        '--log',
        'info',
      ];

      // 如果指定了组合 ID，添加到参数
      if (compositionId) {
        args.push('--props', JSON.stringify({ compositionId }));
      }

      console.error(`[PreviewServerManager] Starting: npx ${args.join(' ')}`);

      // 启动子进程
      const childProcess = spawn('npx', args, {
        cwd: projectPath,
        stdio: ['ignore', 'pipe', 'pipe'],
        env: {
          ...process.env,
          // 确保使用正确的 Node 环境
          NODE_ENV: 'development',
        },
        // 在 Windows 上需要 shell
        shell: process.platform === 'win32',
      });

      let started = false;
      let output = '';

      // 设置超时
      const timeout = setTimeout(() => {
        if (!started) {
          childProcess.kill();
          reject(createPreviewFailedError('服务器启动超时'));
        }
      }, SERVER_START_TIMEOUT);

      // 监听 stdout
      childProcess.stdout?.on('data', (data: Buffer) => {
        const text = data.toString();
        output += text;
        console.error(`[Remotion Studio] ${text.trim()}`);

        // 检查是否启动成功
        // Remotion Studio 启动后会输出类似 "Server listening on port XXXX" 的信息
        if (!started && (
          text.includes('Server listening') ||
          text.includes('http://localhost:') ||
          text.includes(`http://127.0.0.1:${port}`) ||
          text.includes('Ready')
        )) {
          started = true;
          clearTimeout(timeout);

          const instance: InternalPreviewInstance = {
            projectId,
            port,
            url: `http://localhost:${port}`,
            pid: childProcess.pid ?? 0,
            process: childProcess,
          };

          resolve(instance);
        }
      });

      // 监听 stderr
      childProcess.stderr?.on('data', (data: Buffer) => {
        const text = data.toString();
        output += text;
        console.error(`[Remotion Studio Error] ${text.trim()}`);
      });

      // 监听进程退出
      childProcess.on('exit', (code, signal) => {
        if (!started) {
          clearTimeout(timeout);
          reject(createPreviewFailedError(
            `服务器进程意外退出 (code: ${code}, signal: ${signal})`,
            { path: projectPath }
          ));
        }
      });

      // 监听错误
      childProcess.on('error', (error) => {
        if (!started) {
          clearTimeout(timeout);
          reject(createPreviewFailedError(error.message, { path: projectPath }));
        }
      });

      // 如果 30 秒后仍未检测到启动成功的信号，但进程仍在运行，则假设启动成功
      // 这是一个后备机制，因为不同版本的 Remotion 可能输出不同的启动信息
      setTimeout(() => {
        if (!started && childProcess.pid && !childProcess.killed) {
          started = true;
          clearTimeout(timeout);

          const instance: InternalPreviewInstance = {
            projectId,
            port,
            url: `http://localhost:${port}`,
            pid: childProcess.pid,
            process: childProcess,
          };

          console.error(`[PreviewServerManager] Assuming server started (fallback)`);
          resolve(instance);
        }
      }, 10000);
    });
  }

  /**
   * 查找 Remotion 入口文件
   *
   * @param projectPath 项目路径
   * @returns 入口文件路径
   */
  private findEntryPoint(projectPath: string): string {
    // 尝试常见的入口文件位置
    const possibleEntryPoints = [
      join(projectPath, 'src', 'Root.tsx'),
      join(projectPath, 'src', 'index.tsx'),
      join(projectPath, 'Root.tsx'),
      join(projectPath, 'index.tsx'),
      join(projectPath, 'remotion', 'Root.tsx'),
      join(projectPath, 'remotion', 'index.tsx'),
    ];

    for (const entryPoint of possibleEntryPoints) {
      if (existsSync(entryPoint)) {
        return entryPoint;
      }
    }

    // 如果项目中没有找到入口文件，尝试使用 @sprouty-ai/video 包的 Root
    try {
      const videoPackageRoot = require.resolve('@sprouty-ai/video');
      const videoPackageDir = dirname(videoPackageRoot);
      const defaultEntryPoint = join(videoPackageDir, 'Root.tsx');
      if (existsSync(defaultEntryPoint)) {
        return defaultEntryPoint;
      }
    } catch {
      // 忽略 video 包未找到的情况
    }

    // 回退到项目路径本身
    return projectPath;
  }

  /**
   * 优雅停止进程
   *
   * @param process 子进程
   */
  private async stopProcess(process: ChildProcess): Promise<void> {
    return new Promise((resolve, reject) => {
      if (!process.pid || process.killed) {
        resolve();
        return;
      }

      // 设置超时
      const timeout = setTimeout(() => {
        // 超时后强制终止
        try {
          process.kill('SIGKILL');
        } catch {
          // 忽略错误
        }
        resolve();
      }, SERVER_STOP_TIMEOUT);

      // 监听退出事件
      process.once('exit', () => {
        clearTimeout(timeout);
        resolve();
      });

      // 发送 SIGTERM 信号
      try {
        process.kill('SIGTERM');
      } catch (error) {
        clearTimeout(timeout);
        reject(error);
      }
    });
  }

  /**
   * 转换为公共预览实例（不包含进程引用）
   *
   * @param instance 内部预览实例
   * @returns 公共预览实例
   */
  private toPublicInstance(instance: InternalPreviewInstance): PreviewInstance {
    return {
      projectId: instance.projectId,
      port: instance.port,
      url: instance.url,
      pid: instance.pid,
    };
  }

  /**
   * 设置进程退出时的清理
   */
  private setupCleanup(): void {
    // 在进程退出时停止所有服务器
    const cleanup = async () => {
      console.error('[PreviewServerManager] Cleaning up...');
      await this.stopAll();
    };

    // 注册退出处理器
    process.on('exit', () => {
      // 同步清理（exit 事件中不能使用异步）
      const instances = Array.from(this.activeServers.values());
      for (const instance of instances) {
        try {
          instance.process.kill('SIGKILL');
        } catch {
          // 忽略错误
        }
      }
    });

    // 注册信号处理器
    process.on('SIGINT', async () => {
      await cleanup();
      process.exit(0);
    });

    process.on('SIGTERM', async () => {
      await cleanup();
      process.exit(0);
    });
  }

  // ==========================================================================
  // 静态工厂方法
  // ==========================================================================

  /**
   * 创建 PreviewServerManager 实例
   * @returns PreviewServerManager 实例
   */
  static create(): PreviewServerManager {
    return new PreviewServerManager();
  }
}

// ============================================================================
// 导出单例实例
// ============================================================================

/**
 * 默认的预览服务器管理器实例
 */
export const previewServerManager = PreviewServerManager.create();
