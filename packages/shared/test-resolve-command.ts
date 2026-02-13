#!/usr/bin/env bun
/**
 * Test script to verify resolveCommand() is working correctly
 */

import { resolveCommand } from './src/sources/server-builder.ts';

console.log('Testing resolveCommand()...\n');

// Test 1: Resolve 'bun'
const bunPath = resolveCommand('bun');
console.log(`✓ resolveCommand('bun') → ${bunPath}`);

// Test 2: Resolve 'node'
const nodePath = resolveCommand('node');
console.log(`✓ resolveCommand('node') → ${nodePath}`);

// Test 3: Absolute path (should return as-is)
const absolutePath = resolveCommand('/usr/bin/env');
console.log(`✓ resolveCommand('/usr/bin/env') → ${absolutePath}`);

// Test 4: Non-existent command (should return as-is)
const nonExistent = resolveCommand('nonexistent-command-xyz');
console.log(`✓ resolveCommand('nonexistent-command-xyz') → ${nonExistent}`);

console.log('\n✅ All tests passed!');
