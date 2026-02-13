#!/usr/bin/env bun

/**
 * 测试 MCP Source 路径解析
 */

import { join } from 'path';
import { existsSync } from 'fs';

// 模拟 resolveCwd 函数
function resolveCwd(cwd: string | undefined, workspaceRootPath: string): string | undefined {
  if (!cwd) return undefined;

  // app: 前缀 - 相对于应用安装目录
  if (cwd.startsWith('app:')) {
    const relativePath = cwd.slice(4); // Remove 'app:' prefix
    // 在开发/测试环境中，使用 dist 目录模拟打包后的结构
    const appRoot = join(__dirname, 'dist');
    const resolved = join(appRoot, relativePath);
    console.log(`[resolveCwd] Resolved 'app:${relativePath}' → '${resolved}'`);
    return resolved;
  }

  // 绝对路径 - 直接使用
  if (cwd.startsWith('/')) {
    return cwd;
  }

  // 相对路径 - 相对于工作区根目录
  const resolved = join(workspaceRootPath, cwd);
  console.log(`[resolveCwd] Resolved relative path '${cwd}' → '${resolved}'`);
  return resolved;
}

// 模拟 resolveCommand 函数
function resolveCommand(command: string): string {
  console.log(`[resolveCommand] Resolving command: ${command}`);

  // 如果是绝对路径，直接返回
  if (command.startsWith('/')) {
    console.log(`[resolveCommand] Using absolute path: ${command}`);
    return command;
  }

  // 尝试使用 which 查找
  try {
    const { execSync } = require('child_process');
    const { homedir } = require('os');

    const extendedPath = [
      process.env.PATH || '',
      join(homedir(), '.bun', 'bin'),
      join(homedir(), '.deno', 'bin'),
      '/usr/local/bin',
      '/opt/homebrew/bin',
    ].join(':');

    const result = execSync(`which ${command}`, {
      env: { ...process.env, PATH: extendedPath },
      encoding: 'utf-8',
    }).trim();

    if (result && existsSync(result)) {
      console.log(`[resolveCommand] Found via which: ${result}`);
      return result;
    }
  } catch (err) {
    console.log(`[resolveCommand] which failed, trying well-known paths...`);
  }

  // 检查 well-known 安装位置
  const { homedir } = require('os');
  const wellKnownPaths = [
    join('/opt/homebrew/bin', command),
    join(homedir(), '.bun', 'bin', command),
    join(homedir(), '.deno', 'bin', command),
    join('/usr/local/bin', command),
  ];

  for (const path of wellKnownPaths) {
    if (existsSync(path)) {
      console.log(`[resolveCommand] Found at well-known path: ${path}`);
      return path;
    }
  }

  console.log(`[resolveCommand] Using original command: ${command}`);
  return command;
}

// 测试场景
console.log('=== MCP Source 路径解析测试 ===\n');

// 场景 1: 打包应用中的 app: 前缀
console.log('场景 1: 打包应用 (app: 前缀)');
const appCwd = 'app:resources/video';
const workspaceRoot = '/Users/mac/.sprouty-ai/workspaces/test-workspace';
const resolvedAppPath = resolveCwd(appCwd, workspaceRoot);
console.log(`输入: ${appCwd}`);
console.log(`输出: ${resolvedAppPath}`);
console.log(`存在: ${existsSync(resolvedAppPath || '')}`);
console.log();

// 场景 2: 开发环境中的相对路径
console.log('场景 2: 开发环境 (相对路径)');
const devCwd = '../../packages/video';
const resolvedDevPath = resolveCwd(devCwd, workspaceRoot);
console.log(`输入: ${devCwd}`);
console.log(`输出: ${resolvedDevPath}`);
console.log();

// 场景 3: 绝对路径
console.log('场景 3: 绝对路径');
const absCwd = '/usr/local/lib/video-mcp';
const resolvedAbsPath = resolveCwd(absCwd, workspaceRoot);
console.log(`输入: ${absCwd}`);
console.log(`输出: ${resolvedAbsPath}`);
console.log();

// 场景 4: 命令路径解析
console.log('场景 4: 命令路径解析');
const bunCommand = resolveCommand('bun');
console.log(`命令: bun`);
console.log(`解析: ${bunCommand}`);
console.log(`存在: ${existsSync(bunCommand)}`);
console.log();

// 验证 MCP 服务器入口文件
console.log('场景 5: 验证 MCP 服务器入口');
if (resolvedAppPath) {
  const mcpEntryPoint = join(resolvedAppPath, 'src/mcp-server/index.ts');
  console.log(`入口文件: ${mcpEntryPoint}`);
  console.log(`存在: ${existsSync(mcpEntryPoint)}`);

  if (existsSync(mcpEntryPoint)) {
    console.log('✅ MCP 服务器入口文件存在');
  } else {
    console.log('❌ MCP 服务器入口文件不存在');
  }
}

console.log('\n=== 测试完成 ===');
