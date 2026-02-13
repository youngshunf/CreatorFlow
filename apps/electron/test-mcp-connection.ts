#!/usr/bin/env bun

/**
 * 测试 Video MCP 服务器连接
 */

import { spawn } from 'child_process';
import { join } from 'path';

const bunPath = '/Users/mac/.bun/bin/bun';
const videoCwd = join(__dirname, 'dist/resources/video');
const mcpScript = 'src/mcp-server/index.ts';

console.log('=== Video MCP 服务器连接测试 ===\n');
console.log(`Bun 路径: ${bunPath}`);
console.log(`工作目录: ${videoCwd}`);
console.log(`脚本路径: ${mcpScript}`);
console.log();

// 启动 MCP 服务器
console.log('启动 MCP 服务器...');
const child = spawn(bunPath, ['run', mcpScript], {
  cwd: videoCwd,
  stdio: ['pipe', 'pipe', 'pipe'],
});

let initTimeout: NodeJS.Timeout;
let hasReceivedResponse = false;

// 监听 stdout（JSON-RPC 响应）
child.stdout.on('data', (data) => {
  const output = data.toString();
  console.log('收到响应:', output);
  hasReceivedResponse = true;

  try {
    const response = JSON.parse(output);
    if (response.result?.serverInfo) {
      console.log('\n✅ MCP 服务器连接成功！');
      console.log('服务器信息:', JSON.stringify(response.result.serverInfo, null, 2));

      if (response.result.capabilities?.tools) {
        console.log('\n可用工具数量:', response.result.capabilities.tools.length || '未知');
      }

      cleanup();
    }
  } catch (err) {
    // 可能是非 JSON 输出，继续等待
  }
});

// 监听 stderr（错误日志）
child.stderr.on('data', (data) => {
  const error = data.toString();
  if (!error.includes('[MCP Video Server]')) {
    console.error('错误输出:', error);
  }
});

// 监听进程退出
child.on('exit', (code, signal) => {
  console.log(`\n进程退出: code=${code}, signal=${signal}`);
  if (!hasReceivedResponse) {
    console.log('❌ 未收到服务器响应');
  }
  process.exit(code || 0);
});

// 监听错误
child.on('error', (err) => {
  console.error('❌ 启动失败:', err);
  process.exit(1);
});

// 发送初始化请求
console.log('发送初始化请求...\n');
const initRequest = {
  jsonrpc: '2.0',
  id: 1,
  method: 'initialize',
  params: {
    protocolVersion: '2024-11-05',
    capabilities: {},
    clientInfo: {
      name: 'test-client',
      version: '1.0.0',
    },
  },
};

child.stdin.write(JSON.stringify(initRequest) + '\n');

// 设置超时
initTimeout = setTimeout(() => {
  if (!hasReceivedResponse) {
    console.log('❌ 初始化超时（10秒）');
    cleanup();
  }
}, 10000);

function cleanup() {
  clearTimeout(initTimeout);
  child.kill();
  setTimeout(() => process.exit(0), 100);
}

// 处理 Ctrl+C
process.on('SIGINT', () => {
  console.log('\n中断测试...');
  cleanup();
});
