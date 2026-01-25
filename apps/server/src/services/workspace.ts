import { join } from 'path'
import { mkdir, access, readdir, rm } from 'fs/promises'
import { existsSync } from 'fs'

const DATA_DIR = process.env.DATA_DIR || '/data/users'

export class WorkspaceService {
  /**
   * Get user's data root directory
   */
  static getUserRoot(userId: string): string {
    return join(DATA_DIR, userId)
  }
  
  /**
   * Get user's workspace path
   */
  static getWorkspacePath(userId: string): string {
    return join(this.getUserRoot(userId), 'workspace')
  }
  
  /**
   * Ensure user workspace exists and is initialized
   */
  static async ensureUserWorkspace(userId: string): Promise<string> {
    const workspacePath = this.getWorkspacePath(userId)
    
    // Create workspace directory
    await mkdir(workspacePath, { recursive: true })
    
    // Check if workspace is already initialized
    const configPath = join(workspacePath, 'config.json')
    if (!existsSync(configPath)) {
      // Initialize with default structure
      await this.initializeWorkspace(workspacePath, userId)
    }
    
    return workspacePath
  }
  
  /**
   * Initialize workspace with default structure
   */
  private static async initializeWorkspace(workspacePath: string, userId: string): Promise<void> {
    // Create subdirectories
    const dirs = ['sessions', 'sources', 'skills', 'statuses', 'permissions']
    await Promise.all(
      dirs.map(dir => mkdir(join(workspacePath, dir), { recursive: true }))
    )
    
    // Create default config
    const config = {
      id: userId,
      name: 'Default Workspace',
      rootPath: workspacePath,
      createdAt: new Date().toISOString(),
    }
    
    await Bun.write(
      join(workspacePath, 'config.json'),
      JSON.stringify(config, null, 2)
    )
    
    // Create default permissions (cloud safe mode)
    const permissions = {
      version: '1.0',
      blockedTools: [],
      readOnlyBashPatterns: [
        { pattern: '^ls( |$)', description: 'List files' },
        { pattern: '^cat ', description: 'Read file' },
        { pattern: '^head ', description: 'Read file head' },
        { pattern: '^tail ', description: 'Read file tail' },
        { pattern: '^grep ', description: 'Search in files' },
        { pattern: '^find ', description: 'Find files' },
        { pattern: '^pwd$', description: 'Print working directory' },
        { pattern: '^echo ', description: 'Echo text' },
        { pattern: '^which ', description: 'Locate command' },
      ],
    }
    
    await Bun.write(
      join(workspacePath, 'permissions', 'default.json'),
      JSON.stringify(permissions, null, 2)
    )
  }
  
  /**
   * Check if workspace exists
   */
  static async exists(userId: string): Promise<boolean> {
    try {
      await access(this.getWorkspacePath(userId))
      return true
    } catch {
      return false
    }
  }
  
  /**
   * Delete user workspace (admin only, for account deletion)
   */
  static async deleteWorkspace(userId: string): Promise<void> {
    const userRoot = this.getUserRoot(userId)
    
    if (existsSync(userRoot)) {
      await rm(userRoot, { recursive: true, force: true })
    }
  }
  
  /**
   * Get workspace disk usage
   */
  static async getDiskUsage(userId: string): Promise<number> {
    const workspacePath = this.getWorkspacePath(userId)
    
    if (!existsSync(workspacePath)) {
      return 0
    }
    
    let totalSize = 0
    
    async function calculateSize(dirPath: string): Promise<void> {
      const entries = await readdir(dirPath, { withFileTypes: true })
      
      for (const entry of entries) {
        const fullPath = join(dirPath, entry.name)
        
        if (entry.isDirectory()) {
          await calculateSize(fullPath)
        } else {
          const file = Bun.file(fullPath)
          totalSize += file.size
        }
      }
    }
    
    await calculateSize(workspacePath)
    
    return totalSize
  }
}
