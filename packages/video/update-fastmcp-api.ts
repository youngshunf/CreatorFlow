#!/usr/bin/env bun
/**
 * 批量更新 FastMCP 工具注册代码
 * 从旧 API: mcp.tool(name, description, schema, handler)
 * 到新 API: mcp.addTool({ name, description, parameters, execute })
 */

import { readFileSync, writeFileSync } from 'fs';
import { glob } from 'glob';

const toolFiles = glob.sync('src/mcp-server/tools/*.ts', {
  ignore: ['**/*.test.ts', '**/index.ts']
});

console.log(`Found ${toolFiles.length} tool files to update:\n`);

for (const file of toolFiles) {
  console.log(`Processing: ${file}`);

  let content = readFileSync(file, 'utf-8');
  let modified = false;

  // 匹配模式: mcp.tool('name', 'description', Schema, handler)
  const toolPattern = /mcp\.tool\(\s*['"]([^'"]+)['"]\s*,\s*['"]([^'"]+)['"]\s*,\s*(\w+)\s*,\s*(\w+)\s*\)/g;

  const matches = [...content.matchAll(toolPattern)];

  if (matches.length > 0) {
    for (const match of matches) {
      const [fullMatch, name, description, schema, handler] = match;

      const newCode = `mcp.addTool({
    name: '${name}',
    description: '${description}',
    parameters: ${schema},
    execute: ${handler},
  })`;

      content = content.replace(fullMatch, newCode);
      modified = true;

      console.log(`  ✓ Updated: ${name}`);
    }

    if (modified) {
      writeFileSync(file, content, 'utf-8');
      console.log(`  ✅ Saved: ${file}\n`);
    }
  } else {
    console.log(`  ⏭️  No changes needed\n`);
  }
}

console.log('✅ All files processed!');
