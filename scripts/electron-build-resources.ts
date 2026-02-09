/**
 * Cross-platform resources copy script
 */

import { existsSync, cpSync } from "fs";
import { join } from "path";

const ROOT_DIR = join(import.meta.dir, "..");
const ELECTRON_DIR = join(ROOT_DIR, "apps/electron");
const SHARED_DIR = join(ROOT_DIR, "packages/shared");

const srcDir = join(ELECTRON_DIR, "resources");
const destDir = join(ELECTRON_DIR, "dist/resources");

if (existsSync(srcDir)) {
  cpSync(srcDir, destDir, { recursive: true, force: true });
  console.log("📦 Copied resources to dist");
} else {
  console.log("⚠️ No resources directory found");
}

// Copy bundled-skills for app initialization
const bundledSkillsSrc = join(SHARED_DIR, "src/apps/bundled-skills");
const bundledSkillsDest = join(ELECTRON_DIR, "dist/resources/bundled-skills");

if (existsSync(bundledSkillsSrc)) {
  cpSync(bundledSkillsSrc, bundledSkillsDest, { recursive: true, force: true });
  console.log("📦 Copied bundled-skills to dist");
} else {
  console.log("⚠️ No bundled-skills directory found");
}

// Copy bundled-plugins for app initialization
const bundledPluginsSrc = join(SHARED_DIR, "src/apps/bundled-plugins");
const bundledPluginsDest = join(ELECTRON_DIR, "dist/resources/bundled-plugins");

if (existsSync(bundledPluginsSrc)) {
  cpSync(bundledPluginsSrc, bundledPluginsDest, { recursive: true, force: true });
  console.log("📦 Copied bundled-plugins to dist");
} else {
  console.log("⚠️ No bundled-plugins directory found");
}
