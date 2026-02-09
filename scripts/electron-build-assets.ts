/**
 * Cross-platform assets copy script
 *
 * Copies shared assets into the Electron dist directory:
 * 1. Documentation files (packages/shared/assets/docs/)
 * 2. Theme files (apps/electron/resources/themes/)
 *
 * At runtime, these files are installed to ~/.sprouty-ai/
 */

import { existsSync, cpSync, mkdirSync } from "fs";
import { join } from "path";

const ROOT_DIR = join(import.meta.dir, "..");
const ELECTRON_DIR = join(ROOT_DIR, "apps/electron");

// Copy documentation assets
const docsSourceDir = join(ROOT_DIR, "packages/shared/assets/docs");
const docsDestDir = join(ELECTRON_DIR, "dist/assets/docs");

if (existsSync(docsSourceDir)) {
  mkdirSync(join(ELECTRON_DIR, "dist/assets"), { recursive: true });
  cpSync(docsSourceDir, docsDestDir, { recursive: true, force: true });
  console.log("📦 Copied doc assets to dist");
} else {
  console.log("⚠️ No shared assets/docs directory found");
}

// Copy theme files
const themesSourceDir = join(ELECTRON_DIR, "resources/themes");
const themesDestDir = join(ELECTRON_DIR, "dist/assets/themes");

if (existsSync(themesSourceDir)) {
  mkdirSync(join(ELECTRON_DIR, "dist/assets"), { recursive: true });
  cpSync(themesSourceDir, themesDestDir, { recursive: true, force: true });
  console.log("🎨 Copied theme assets to dist");
} else {
  console.log("⚠️ No resources/themes directory found");
}
