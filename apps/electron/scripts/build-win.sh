#!/bin/bash
# Build script for Windows NSIS installer (cross-compile from macOS/Linux)
# Usage: bash scripts/build-win.sh [development|staging|production]
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ELECTRON_DIR="$(dirname "$SCRIPT_DIR")"
ROOT_DIR="$(dirname "$(dirname "$ELECTRON_DIR")")"

# Helper function to check required file/directory exists
require_path() {
    local path="$1"
    local description="$2"
    local hint="$3"

    if [ ! -e "$path" ]; then
        echo "ERROR: $description not found at $path"
        [ -n "$hint" ] && echo "$hint"
        exit 1
    fi
}

# Parse arguments
MODE="${1:-development}"
case "$MODE" in
    development|staging|production)
        ;;
    *)
        echo "Invalid mode: $MODE"
        echo "Usage: build-win.sh [development|staging|production]"
        exit 1
        ;;
esac

# Set environment variable for Vite build mode
export VITE_APP_ENV="$MODE"
export APP_ENV="$MODE"

# Configuration
BUN_VERSION="bun-v1.3.5"  # Pinned version for reproducible builds

# Disable SSL verification if behind corporate proxy/VPN (common in China)
# Remove this if you don't need it
export NODE_TLS_REJECT_UNAUTHORIZED=0

echo "=== Building CreatorFlow Windows Installer ($MODE) using electron-builder ==="

# 1. Clean previous build artifacts
echo "Cleaning previous builds..."
rm -rf "$ELECTRON_DIR/vendor"
rm -rf "$ELECTRON_DIR/node_modules/@anthropic-ai"
rm -rf "$ELECTRON_DIR/packages"
rm -rf "$ELECTRON_DIR/release"

# 2. Install dependencies
echo "Installing dependencies..."
cd "$ROOT_DIR"
bun install

# 3. Download Bun binary for Windows (with cache support)
# Use baseline build - works on all x64 CPUs (no AVX2 requirement)
echo "Preparing Bun ${BUN_VERSION} for Windows x64 (baseline)..."
mkdir -p "$ELECTRON_DIR/vendor/bun"
BUN_DOWNLOAD="bun-windows-x64-baseline"
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Cache directory
CACHE_DIR="$HOME/.cache/sprouty-ai/bun"
CACHED_ZIP="$CACHE_DIR/${BUN_VERSION}-${BUN_DOWNLOAD}.zip"
mkdir -p "$CACHE_DIR"

DOWNLOAD_SUCCESS=false

# Check if cached version exists
if [ -f "$CACHED_ZIP" ]; then
    echo "Found cached Bun binary at $CACHED_ZIP"
    echo "Extracting from cache..."
    unzip -o "$CACHED_ZIP" -d "$TEMP_DIR"
    if [ -f "$TEMP_DIR/${BUN_DOWNLOAD}/bun.exe" ]; then
        cp "$TEMP_DIR/${BUN_DOWNLOAD}/bun.exe" "$ELECTRON_DIR/vendor/bun/"
        DOWNLOAD_SUCCESS=true
        echo "Using cached Bun ${BUN_VERSION} for Windows"
    else
        echo "Cached file is corrupted, will re-download..."
        rm -f "$CACHED_ZIP"
    fi
fi

# Download if not cached or cache is corrupted
if [ "$DOWNLOAD_SUCCESS" = false ]; then
    ZIP_URL="https://github.com/oven-sh/bun/releases/download/${BUN_VERSION}/${BUN_DOWNLOAD}.zip"
    CHECKSUM_URL="https://github.com/oven-sh/bun/releases/download/${BUN_VERSION}/SHASUMS256.txt"

    echo "Downloading from $ZIP_URL..."
    if curl -fSL --connect-timeout 30 --max-time 180 "$ZIP_URL" -o "$TEMP_DIR/${BUN_DOWNLOAD}.zip" && \
       curl -fSL --connect-timeout 30 --max-time 60 "$CHECKSUM_URL" -o "$TEMP_DIR/SHASUMS256.txt"; then
        # Verify checksum
        echo "Verifying checksum..."
        cd "$TEMP_DIR"
        if grep "${BUN_DOWNLOAD}.zip" SHASUMS256.txt | shasum -a 256 -c -; then
            cd - > /dev/null
            # Extract and install
            echo "Extracting Bun..."
            unzip -o "$TEMP_DIR/${BUN_DOWNLOAD}.zip" -d "$TEMP_DIR"
            cp "$TEMP_DIR/${BUN_DOWNLOAD}/bun.exe" "$ELECTRON_DIR/vendor/bun/"
            # Save to cache
            echo "Saving to cache..."
            cp "$TEMP_DIR/${BUN_DOWNLOAD}.zip" "$CACHED_ZIP"
            DOWNLOAD_SUCCESS=true
            echo "Downloaded and cached Bun ${BUN_VERSION} for Windows successfully"
        else
            cd - > /dev/null
            echo "Checksum verification failed!"
        fi
    else
        echo "Download failed!"
    fi
fi

if [ "$DOWNLOAD_SUCCESS" = false ]; then
    echo "ERROR: Could not download Windows Bun binary"
    echo "Please check your network connection and try again."
    echo "Or manually download from: https://github.com/oven-sh/bun/releases/download/${BUN_VERSION}/${BUN_DOWNLOAD}.zip"
    echo "And place it at: $CACHED_ZIP"
    exit 1
fi

# Verify bun.exe exists
BUN_EXE="$ELECTRON_DIR/vendor/bun/bun.exe"
if [ ! -f "$BUN_EXE" ]; then
    echo "ERROR: bun.exe not found at $BUN_EXE"
    exit 1
fi
echo "Bun extracted to: $BUN_EXE"
ls -la "$BUN_EXE"

# 4. Copy SDK from root node_modules (monorepo hoisting)
SDK_SOURCE="$ROOT_DIR/node_modules/@anthropic-ai/claude-agent-sdk"
SDK_DEST="$ELECTRON_DIR/node_modules/@anthropic-ai/claude-agent-sdk"
require_path "$SDK_SOURCE" "SDK" "Run 'bun install' from the repository root first."
echo "Copying SDK..."
mkdir -p "$ELECTRON_DIR/node_modules/@anthropic-ai"
# Remove existing symlink or directory first to avoid "identical" error
rm -rf "$SDK_DEST"
# Use -L to follow symlinks and copy actual content
cp -rL "$SDK_SOURCE" "$SDK_DEST"

# 5. Copy interceptor
INTERCEPTOR_SOURCE="$ROOT_DIR/packages/shared/src/network-interceptor.ts"
require_path "$INTERCEPTOR_SOURCE" "Interceptor" "Ensure packages/shared/src/network-interceptor.ts exists."
echo "Copying interceptor..."
mkdir -p "$ELECTRON_DIR/packages/shared/src"
cp "$INTERCEPTOR_SOURCE" "$ELECTRON_DIR/packages/shared/src/"

# 6. Build Electron app
echo "Building Electron app for $MODE environment..."
cd "$ROOT_DIR"

# Use appropriate build command based on MODE
if [ "$MODE" = "staging" ]; then
    echo "Building for STAGING environment..."
    bun run electron:build:staging
elif [ "$MODE" = "production" ]; then
    echo "Building for PRODUCTION environment..."
    bun run electron:build:production
else
    echo "Building for DEVELOPMENT environment..."
    bun run electron:build
fi

# 8. Package with electron-builder
echo "Packaging app with electron-builder..."
cd "$ELECTRON_DIR"

# Verify bun.exe still exists before packaging (sanity check)
if [ ! -f "$BUN_EXE" ]; then
    echo "ERROR: bun.exe disappeared before packaging!"
    exit 1
fi

npx electron-builder --win --x64

# 9. Verify the installer was built
# Read version from package.json
ELECTRON_VERSION=$(cat "$ELECTRON_DIR/package.json" | grep '"version"' | head -1 | sed 's/.*"version": *"\([^"]*\)".*/\1/')
# electron-builder.yml uses artifactName: 智小芽_v${version}-x64.exe
INSTALLER_NAME="智小芽_v${ELECTRON_VERSION}-x64.exe"
INSTALLER_PATH="$ELECTRON_DIR/release/$INSTALLER_NAME"

if [ ! -f "$INSTALLER_PATH" ]; then
    echo "ERROR: Expected installer not found at $INSTALLER_PATH"
    echo "Contents of release directory:"
    ls -la "$ELECTRON_DIR/release/"
    exit 1
fi

echo ""
echo "=== Build Complete ==="
echo "Installer: $INSTALLER_PATH"
echo "Size: $(du -h "$INSTALLER_PATH" | cut -f1)"
