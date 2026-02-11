/**
 * Orchestrates the full production build for the Electron app.
 * Runs each build step sequentially.
 */

import { spawn } from "bun";
import { join } from "path";

const ROOT_DIR = join(import.meta.dir, "..");

const steps = [
  { name: "Main Process", script: "scripts/electron-build-main.ts" },
  { name: "Preload", script: "scripts/electron-build-preload.ts" },
  { name: "Renderer", script: "scripts/electron-build-renderer.ts" },
  { name: "Assets", script: "scripts/electron-build-assets.ts" },
];

async function build() {
  console.log("🚀 Starting full build...\n");

  for (const step of steps) {
    console.log(`📦 Building ${step.name}...`);
    const proc = spawn({
      cmd: ["bun", "run", join(ROOT_DIR, step.script)],
      cwd: ROOT_DIR,
      stdout: "inherit",
      stderr: "inherit",
    });

    const exitCode = await proc.exited;
    if (exitCode !== 0) {
      console.error(`\n❌ ${step.name} build failed (exit code ${exitCode})`);
      process.exit(exitCode);
    }
    console.log(`✅ ${step.name} done\n`);
  }

  console.log("🎉 Full build completed successfully!");
}

build();
