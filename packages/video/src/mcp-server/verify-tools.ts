#!/usr/bin/env bun

/**
 * 验证 MCP Server 工具注册
 *
 * 这个脚本验证所有 MCP 工具是否正确注册
 */

import { TOOL_LIST } from "./tools/index";

console.log("=".repeat(60));
console.log("MCP Video Server 工具注册验证");
console.log("=".repeat(60));
console.log();

console.log(`总工具数: ${TOOL_LIST.length}`);
console.log();

// 按类别分组
const categories = new Map<string, Array<(typeof TOOL_LIST)[number]>>();
for (const tool of TOOL_LIST) {
  const category = tool.category;
  if (!categories.has(category)) {
    categories.set(category, []);
  }
  categories.get(category)!.push(tool);
}

// 显示每个类别的工具
for (const [category, tools] of categories.entries()) {
  console.log(`\n${category.toUpperCase()} (${tools.length} 个工具):`);
  console.log("-".repeat(60));
  for (const tool of tools) {
    console.log(`  ✓ ${tool.name}`);
    console.log(`    ${tool.description}`);
  }
}

console.log();
console.log("=".repeat(60));
console.log("验证完成！");
console.log("=".repeat(60));

// 验证新增的工具
const newTools = [
  "video_list_available_assets",
  "video_validate_composition",
  "video_get_render_status",
];

console.log();
console.log("新增工具验证:");
console.log("-".repeat(60));

for (const toolName of newTools) {
  const found = TOOL_LIST.find((t) => t.name === toolName);
  if (found) {
    console.log(`  ✓ ${toolName} - 已注册`);
  } else {
    console.log(`  ✗ ${toolName} - 未找到`);
    process.exit(1);
  }
}

console.log();
console.log("所有新工具已成功注册！");
