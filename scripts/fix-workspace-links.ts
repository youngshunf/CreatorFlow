/**
 * Fix workspace symlinks that bun sometimes fails to create automatically.
 * This script runs after `bun install` via postinstall hook.
 */

import { mkdirSync, symlinkSync, existsSync, lstatSync, unlinkSync } from 'node:fs'
import { join, dirname } from 'node:path'
import { fileURLToPath } from 'node:url'

const __dirname = dirname(fileURLToPath(import.meta.url))
const rootDir = join(__dirname, '..')
const nodeModulesDir = join(rootDir, 'node_modules')
const scopeDir = join(nodeModulesDir, '@creator-flow')

// Workspace packages to link
const packages = [
  { name: 'core', path: '../../packages/core' },
  { name: 'shared', path: '../../packages/shared' },
  { name: 'ui', path: '../../packages/ui' },
  { name: 'mermaid', path: '../../packages/mermaid' },
]

function main() {
  // Ensure node_modules exists
  if (!existsSync(nodeModulesDir)) {
    console.log('⚠️  node_modules not found, skipping workspace link fix')
    return
  }

  // Create @creator-flow scope directory
  if (!existsSync(scopeDir)) {
    mkdirSync(scopeDir, { recursive: true })
    console.log('📁 Created @creator-flow scope directory')
  }

  let linksCreated = 0

  for (const pkg of packages) {
    const linkPath = join(scopeDir, pkg.name)
    const targetPath = pkg.path

    // Check if the target package exists
    const absoluteTargetPath = join(scopeDir, targetPath)
    if (!existsSync(absoluteTargetPath)) {
      console.log(`⚠️  Package ${pkg.name} not found at ${absoluteTargetPath}`)
      continue
    }

    // Remove existing link if it's incorrect
    if (existsSync(linkPath)) {
      const stat = lstatSync(linkPath)
      if (stat.isSymbolicLink()) {
        // Link already exists, check if it points to the right place
        continue
      }
      // Not a symlink, remove it
      unlinkSync(linkPath)
    }

    // Create symlink
    try {
      symlinkSync(targetPath, linkPath)
      linksCreated++
      console.log(`🔗 Linked @creator-flow/${pkg.name} -> ${targetPath}`)
    } catch (err) {
      console.error(`❌ Failed to link @creator-flow/${pkg.name}:`, err)
    }
  }

  if (linksCreated > 0) {
    console.log(`✅ Fixed ${linksCreated} workspace link(s)`)
  } else {
    console.log('✅ All workspace links are correct')
  }
}

main()
