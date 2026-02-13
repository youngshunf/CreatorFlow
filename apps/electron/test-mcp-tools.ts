#!/usr/bin/env bun

/**
 * 简单测试：验证 Video MCP 服务器工具列表
 */

import { spawn } from 'child_process';
import { join } from 'path';

const bunPath = '/Users/mac/.bun/bin/bun';
const videoCwd = join(__dirname, 'dist/resources/video');
const mcpScript = 'src/mcp-server/index.ts';

console.log('=== Video MCP 工具列表验证 ===\n');

// 启动 MCP 服务器
const child = spawn(bunPath, ['run', mcpScript], {
  cwd: videoCwd,
  stdio: ['pipe', 'pipe', 'pipe'],
});

let toolCount = 0;
const tools: string[] = [];

// 监听 stdout
child.stdout.on('data', (data) => {
  const output = data.toString();

  // 提取工具名称
  const toolMatch = output.match(/- (video_\w+):/g);
  if (toolMatch) {
    toolMatch.forEach(match => {
      const toolName = match.replace(/- /, '').replace(/:/, '');
      if (!tools.includes(toolName)) {
        tools.push(toolName);
      }
    });
  }

  // 提取工具总数
  const totalMatch = output.match(/Total: (\d+) tools/);
  if (totalMatch) {
    toolCount = parseInt(totalMatch[1], 10);
  }
});

// 监听 stderr
child.stderr.on('data', (data) => {
  const error = data.toString();
  // 忽略日志信息
  if (!error.includes('[MCP Video Server]')) {
    console.error('错误:', error);
  }
});

// 等待 2 秒后检查结果
setTimeout(() => {
  child.kill();

  console.log(`\n检测到的工具数量: ${toolCount}`);
  console.log(`\n工具列表 (${tools.length} 个):`);
  tools.forEach((tool, index) => {
    console.log(`  ${index + 1}. ${tool}`);
  });

  if (toolCount === 19 && tools.length >= 18) {
    console.log('\n✅ Video MCP 服务器工具验证成功！');
    console.log(`   - 预期工具数: 18-19`);
    console.log(`   - 实际工具数: ${toolCount}`);
    console.log(`   - 检测到工具: ${tools.length}`);
  } else {
    console.log('\n⚠️  工具数量不匹配');
    console.log(`   - 预期: 18-19`);
    console.log(`   - 实际: ${toolCount}`);
  }

  process.exit(0);
}, 2000);

// 处理 Ctrl+C
process.on('SIGINT', () => {
  child.kill();
  process.exit(0);
});
