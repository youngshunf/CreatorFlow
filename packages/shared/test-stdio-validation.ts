#!/usr/bin/env bun
/**
 * Test script to verify validateStdioMcpConnection() uses resolveCommand()
 */

import { validateStdioMcpConnection } from './src/mcp/validation.ts';

console.log('Testing validateStdioMcpConnection() with command resolution...\n');

// Test with 'bun' command (should auto-resolve to absolute path)
console.log('Testing video-mcp with "bun" command...');
const result = await validateStdioMcpConnection({
  command: 'bun',
  args: ['run', '../../packages/video/src/mcp-server/index.ts'],
  timeout: 10000,
});

console.log('\nValidation Result:');
console.log(`  Success: ${result.success}`);
console.log(`  Error: ${result.error || 'none'}`);
console.log(`  Tools: ${result.tools?.length || 0} tools found`);

if (result.success) {
  console.log('\n✅ Validation passed! Command resolution is working.');
  console.log(`\nAvailable tools (${result.tools?.length}):`);
  result.tools?.forEach((tool, i) => {
    console.log(`  ${i + 1}. ${tool}`);
  });
} else {
  console.log('\n❌ Validation failed!');
  console.log(`Error: ${result.error}`);
}
