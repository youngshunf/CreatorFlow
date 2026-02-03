/**
 * Cross-platform renderer build script
 */

import { spawn } from "bun";
import { existsSync, rmSync } from "fs";
import { join } from "path";

const ROOT_DIR = join(import.meta.dir, "..");
const ELECTRON_DIR = join(ROOT_DIR, "apps/electron");

// Get environment from parent process (prioritize VITE_APP_ENV for consistency with build scripts)
const appEnv = process.env.VITE_APP_ENV || process.env.APP_ENV || "development";
console.log(`📦 Building renderer for environment: ${appEnv}`);

// Clean renderer dist first
const rendererDir = join(ELECTRON_DIR, "dist/renderer");
if (existsSync(rendererDir)) {
  rmSync(rendererDir, { recursive: true, force: true });
}

const proc = spawn({
  cmd: ["bun", "run", "vite", "build", "--config", "apps/electron/vite.config.ts"],
  cwd: ROOT_DIR,
  stdout: "inherit",
  stderr: "inherit",
  env: {
    ...process.env,
    // Pass both for compatibility
    VITE_APP_ENV: appEnv,
    APP_ENV: appEnv,
  },
});

const exitCode = await proc.exited;
process.exit(exitCode);
