// Load user's shell environment first (before other imports that may use env)
// This ensures tools like Homebrew, nvm, etc. are available to the agent
import { loadShellEnv } from './shell-env'
loadShellEnv()

import { app, BrowserWindow, session } from 'electron'
import { createHash } from 'crypto'
import { hostname, homedir } from 'os'
import * as Sentry from '@sentry/electron/main'

// Initialize Sentry error tracking as early as possible after app import.
// Only enabled in production (packaged) builds to avoid noise during development.
// DSN is baked in at build time via esbuild --define (same pattern as OAuth secrets).
//
// NOTE: Source map upload is intentionally disabled. Stack traces in Sentry will show
// bundled/minified code. To enable source map upload in the future:
//   1. Add SENTRY_AUTH_TOKEN, SENTRY_ORG, SENTRY_PROJECT to CI secrets
//   2. Re-enable the @sentry/vite-plugin in vite.config.ts (handles renderer maps)
//   3. Add @sentry/esbuild-plugin to scripts/electron-build-main.ts (handles main process maps)
Sentry.init({
  dsn: process.env.SENTRY_ELECTRON_INGEST_URL,
  environment: app.isPackaged ? 'production' : 'development',
  release: app.getVersion(),
  // Enabled whenever the ingest URL is available — works in both production (baked via CI)
  // and development (injected via .env / 1Password). Filter by environment in Sentry dashboard.
  enabled: !!process.env.SENTRY_ELECTRON_INGEST_URL,

  // Scrub sensitive data before sending to Sentry.
  // Removes authorization headers, API keys/tokens, and credential-like values.
  beforeSend(event) {
    // Scrub request headers (authorization, cookies)
    if (event.request?.headers) {
      const sensitiveHeaders = ['authorization', 'cookie', 'x-api-key']
      for (const header of sensitiveHeaders) {
        if (event.request.headers[header]) {
          event.request.headers[header] = '[REDACTED]'
        }
      }
    }

    // Scrub breadcrumb data that may contain sensitive values
    if (event.breadcrumbs) {
      for (const breadcrumb of event.breadcrumbs) {
        if (breadcrumb.data) {
          for (const key of Object.keys(breadcrumb.data)) {
            const lowerKey = key.toLowerCase()
            if (lowerKey.includes('token') || lowerKey.includes('key') ||
                lowerKey.includes('secret') || lowerKey.includes('password') ||
                lowerKey.includes('credential') || lowerKey.includes('auth')) {
              breadcrumb.data[key] = '[REDACTED]'
            }
          }
        }
      }
    }

    return event
  },
})

// Set anonymous machine ID for Sentry user tracking (no PII — just a hash).
// Uses hostname + homedir to produce a stable per-machine identifier.
const machineId = createHash('sha256').update(hostname() + homedir()).digest('hex').slice(0, 16)
Sentry.setUser({ id: machineId })

import { join } from 'path'
import { existsSync } from 'fs'
import { SessionManager } from './sessions'
import { registerIpcHandlers } from './ipc'
import { createApplicationMenu } from './menu'
import { WindowManager } from './window-manager'
import { loadWindowState, saveWindowState } from './window-state'
import { getWorkspaces, loadStoredConfig } from '@creator-flow/shared/config'
import { initializeDocs } from '@creator-flow/shared/docs'
import { ensureDefaultPermissions } from '@creator-flow/shared/agent/permissions-config'
import { ensureToolIcons } from '@creator-flow/shared/config'
import { setBundledAssetsRoot } from '@creator-flow/shared/utils'
import { handleDeepLink } from './deep-link'
import { registerThumbnailScheme, registerThumbnailHandler } from './thumbnail-protocol'
import log, { isDebugMode, mainLog, getLogFilePath } from './logger'
import { setPerfEnabled, enableDebug } from '@creator-flow/shared/utils'
import { initNotificationService, clearBadgeCount, initBadgeIcon, initInstanceBadge } from './notifications'
import { checkForUpdatesOnLaunch, setWindowManager as setAutoUpdateWindowManager, isUpdating } from './auto-update'
import { registerBundledApps } from '@creator-flow/shared/apps'

// Initialize electron-log for renderer process support
log.initialize()

// Enable debug/perf in dev mode (running from source)
if (isDebugMode) {
  process.env.CRAFT_DEBUG = '1'
  enableDebug()
  setPerfEnabled(true)
}

// Custom URL scheme for deeplinks (e.g., craftagents://auth-complete)
// Supports multi-instance dev: CRAFT_DEEPLINK_SCHEME env var (craftagents1, craftagents2, etc.)
const DEEPLINK_SCHEME = process.env.CRAFT_DEEPLINK_SCHEME || 'craftagents'

let windowManager: WindowManager | null = null
let sessionManager: SessionManager | null = null

// Store pending deep link if app not ready yet (cold start)
let pendingDeepLink: string | null = null

// Set app name early (before app.whenReady) to ensure correct macOS menu bar title
// Supports multi-instance dev: CRAFT_APP_NAME env var (e.g., "CreatorFlow [1]")
app.setName(process.env.CRAFT_APP_NAME || '智小芽')

// Register as default protocol client for craftagents:// URLs
// This must be done before app.whenReady() on some platforms
if (process.defaultApp) {
  // Development mode: need to pass the app path
  if (process.argv.length >= 2) {
    app.setAsDefaultProtocolClient(DEEPLINK_SCHEME, process.execPath, [process.argv[1]])
  }
} else {
  // Production mode
  app.setAsDefaultProtocolClient(DEEPLINK_SCHEME)
}

// Register thumbnail:// custom protocol for file preview thumbnails in the sidebar.
// Must happen before app.whenReady() — Electron requires early scheme registration.
registerThumbnailScheme()

// Handle deeplink on macOS (when app is already running)
app.on('open-url', (event, url) => {
  event.preventDefault()
  mainLog.info('Received deeplink:', url)

  if (windowManager) {
    handleDeepLink(url, windowManager).catch(err => {
      mainLog.error('Failed to handle deep link:', err)
    })
  } else {
    // App not ready - store for later
    pendingDeepLink = url
  }
})

// Handle deeplink on Windows/Linux (single instance check)
const gotTheLock = app.requestSingleInstanceLock()
if (!gotTheLock) {
  app.quit()
} else {
  app.on('second-instance', (_event, commandLine, _workingDirectory) => {
    // Someone tried to run a second instance, we should focus our window.
    // On Windows/Linux, the deeplink is in commandLine
    const url = commandLine.find(arg => arg.startsWith(`${DEEPLINK_SCHEME}://`))
    if (url && windowManager) {
      mainLog.info('Received deeplink from second instance:', url)
      handleDeepLink(url, windowManager).catch(err => {
        mainLog.error('Failed to handle deep link:', err)
      })
    } else if (windowManager) {
      // No deep link - just focus the first window
      const windows = windowManager.getAllWindows()
      if (windows.length > 0) {
        const win = windows[0].window
        if (win.isMinimized()) win.restore()
        win.focus()
      }
    }
  })
}

// Helper to create initial windows on startup
async function createInitialWindows(): Promise<void> {
  if (!windowManager) return

  // Load saved window state
  const savedState = loadWindowState()
  const workspaces = getWorkspaces()
  const validWorkspaceIds = workspaces.map(ws => ws.id)

  if (workspaces.length === 0) {
    // No workspaces configured - create window without workspace (will show onboarding)
    windowManager.createWindow({ workspaceId: '' })
    return
  }

  if (savedState?.windows.length) {
    // Restore windows from saved state
    let restoredCount = 0

    for (const saved of savedState.windows) {
      // Skip invalid workspaces
      if (!validWorkspaceIds.includes(saved.workspaceId)) continue

      // Restore main window with focused mode if it was saved
      mainLog.info(`Restoring window: workspaceId=${saved.workspaceId}, focused=${saved.focused ?? false}, url=${saved.url ?? 'none'}`)
      const win = windowManager.createWindow({
        workspaceId: saved.workspaceId,
        focused: saved.focused,
        restoreUrl: saved.url,
      })
      win.setBounds(saved.bounds)

      restoredCount++
    }

    if (restoredCount > 0) {
      mainLog.info(`Restored ${restoredCount} window(s) from saved state`)
      return
    }
  }

  // Default: open window for first workspace
  windowManager.createWindow({ workspaceId: workspaces[0].id })
  mainLog.info(`Created window for first workspace: ${workspaces[0].name}`)
}

// ============================================================
// Certificate Error Handler for Staging (self-signed certs)
// ============================================================
// In staging/development environments, we may connect to internal servers
// with self-signed certificates. This handler allows those connections.
// In production, this only affects our whitelisted internal domains.
// ============================================================
app.on('certificate-error', (event, _webContents, url, _error, _certificate, callback) => {
  // Allow self-signed certs for internal staging/dev servers
  // TEMPORARY: Also allow api.ai.dcfuture.cn until SSL cert is fixed (CN mismatch)
  const allowedHosts = ['192.168.1.92', 'localhost', '127.0.0.1', 'api.ai.dcfuture.cn']
  try {
    const parsedUrl = new URL(url)
    if (allowedHosts.some(host => parsedUrl.hostname === host)) {
      event.preventDefault()
      callback(true) // Trust this certificate
      return
    }
  } catch {
    // Invalid URL, reject
  }
  callback(false) // Use default behavior (reject invalid certs)
})

app.whenReady().then(async () => {
  // ============================================================
  // CORS Bypass for Cloud API
  // ============================================================
  // Electron's renderer process is subject to CORS restrictions when making
  // fetch requests to external APIs. We intercept responses and add the
  // necessary CORS headers to allow requests to our cloud API endpoints.
  //
  // This is safe because:
  // 1. We only modify responses from our own API domains
  // 2. The main process has full control over what domains are whitelisted
  // 3. contextIsolation is enabled, preventing renderer from accessing Node
  // ============================================================
  const cloudApiDomains = [
    'localhost',               // Development
    '127.0.0.1',               // Development (IP)
    '192.168.1.92',            // Staging (internal network)
    'api.ai.dcfuture.cn',      // Production
  ]

  session.defaultSession.webRequest.onHeadersReceived((details, callback) => {
    // Skip if URL is not provided
    if (!details.url) {
      callback({ responseHeaders: details.responseHeaders })
      return
    }

    try {
      const url = new URL(details.url)
      const isCloudApi = cloudApiDomains.some(domain => url.hostname.includes(domain))

      if (isCloudApi) {
        const headers = details.responseHeaders || {}
        // 检查是否已经有 CORS 头，避免重复添加导致 "multiple values" 错误
        const hasAllowOrigin = Object.keys(headers).some(
          key => key.toLowerCase() === 'access-control-allow-origin'
        )

        if (hasAllowOrigin) {
          // 后端已经返回了 CORS 头，不需要再添加
          callback({ responseHeaders: headers })
        } else {
          // 后端没有返回 CORS 头，添加它们
          callback({
            responseHeaders: {
              ...headers,
              'Access-Control-Allow-Origin': ['*'],
              'Access-Control-Allow-Methods': ['GET, POST, PUT, DELETE, OPTIONS'],
              'Access-Control-Allow-Headers': ['Content-Type, Authorization'],
            },
          })
        }
      } else {
        callback({ responseHeaders: details.responseHeaders })
      }
    } catch {
      // Invalid URL, pass through without modification
      callback({ responseHeaders: details.responseHeaders })
    }
  })

  mainLog.info('CORS bypass configured for cloud API domains:', cloudApiDomains)

  // Register bundled assets root so all seeding functions can find their files
  // (docs, permissions, themes, tool-icons resolve via getBundledAssetsDir)
  setBundledAssetsRoot(__dirname)

  // Initialize bundled docs
  initializeDocs()

  // Register bundled apps (must happen early so apps are available for workspace creation)
  registerBundledApps()

  // Ensure default permissions file exists (copies bundled default.json on first run)
  ensureDefaultPermissions()

  // Seed tool icons to ~/.craft-agent/tool-icons/ (copies bundled SVGs on first run)
  ensureToolIcons()

  // Register thumbnail:// protocol handler (scheme was registered earlier, before app.whenReady)
  registerThumbnailHandler()

  // Note: electron-updater handles pending updates internally via autoInstallOnAppQuit

  // Application menu is created after windowManager initialization (see below)

  // Initialize notification icon for all platforms
  const notificationIconPath = join(__dirname, '../resources/icon.png')
  const notificationIcoPath = join(__dirname, '../resources/icon.ico')
  if (existsSync(notificationIconPath)) {
    initBadgeIcon(notificationIconPath, notificationIcoPath)
  }

  // Set dock icon on macOS (required for dev mode, bundled apps use Info.plist)
  if (process.platform === 'darwin' && app.dock) {
    const dockIconPath = join(__dirname, '../resources/icon.png')
    if (existsSync(dockIconPath)) {
      app.dock.setIcon(dockIconPath)
    }

    // Multi-instance dev: show instance number badge on dock icon
    // CRAFT_INSTANCE_NUMBER is set by detect-instance.sh for numbered folders
    const instanceNum = process.env.CRAFT_INSTANCE_NUMBER
    if (instanceNum) {
      const num = parseInt(instanceNum, 10)
      if (!isNaN(num) && num > 0) {
        initInstanceBadge(num)
      }
    }
  }

  try {
    // Initialize window manager
    windowManager = new WindowManager()

    // Create the application menu (needs windowManager for New Window action)
    createApplicationMenu(windowManager)

    // Initialize session manager
    sessionManager = new SessionManager()
    sessionManager.setWindowManager(windowManager)

    // Initialize notification service
    initNotificationService(windowManager)

    // Register IPC handlers (must happen before window creation)
    registerIpcHandlers(sessionManager, windowManager)

    // Create initial windows (restores from saved state or opens first workspace)
    await createInitialWindows()

    // Initialize auth (must happen after window creation for error reporting)
    await sessionManager.initialize()

    // Set Sentry context tags for error grouping (no PII — just config classification).
    // Runs after init so config and auth state are available.
    try {
      const config = loadStoredConfig()
      const workspaces = getWorkspaces()
      Sentry.setTag('authType', config?.authType ?? 'unknown')
      Sentry.setTag('hasCustomEndpoint', String(!!config?.anthropicBaseUrl))
      Sentry.setTag('model', config?.model ?? 'default')
      Sentry.setTag('customModel', config?.customModel ?? 'none')
      Sentry.setTag('workspaceCount', String(workspaces.length))
    } catch (err) {
      mainLog.warn('Failed to set Sentry context tags:', err)
    }

    // Initialize auto-update (check immediately on launch)
    // Skip in dev mode to avoid replacing /Applications app and launching it instead
    setAutoUpdateWindowManager(windowManager)
    if (app.isPackaged) {
      checkForUpdatesOnLaunch().catch(err => {
        mainLog.error('[auto-update] Launch check failed:', err)
      })
    } else {
      mainLog.info('[auto-update] Skipping auto-update in dev mode')
    }

    // Process pending deep link from cold start
    if (pendingDeepLink) {
      mainLog.info('Processing pending deep link:', pendingDeepLink)
      await handleDeepLink(pendingDeepLink, windowManager)
      pendingDeepLink = null
    }

    mainLog.info('App initialized successfully')
    mainLog.info('Logs at:', getLogFilePath() || 'file logging enabled')
    
    // Log environment configuration for debugging
    try {
      const { getCurrentEnv, getCloudApiUrl, getLlmGatewayUrl } = await import('@creator-flow/shared/config/environments')
      mainLog.info('Environment config:', {
        env: getCurrentEnv(),
        cloudApiUrl: getCloudApiUrl(),
        llmGatewayUrl: getLlmGatewayUrl(),
        isPackaged: app.isPackaged,
      })
    } catch (err) {
      mainLog.warn('Failed to log environment config:', err)
    }
  } catch (error) {
    // Log full error details (error objects may not serialize properly)
    const errorMessage = error instanceof Error 
      ? `${error.message}\n${error.stack}` 
      : String(error)
    mainLog.error('Failed to initialize app:', errorMessage)
    // Continue anyway - the app will show errors in the UI
  }

  // macOS: Re-create window when dock icon is clicked
  app.on('activate', () => {
    if (!windowManager?.hasWindows()) {
      // Open first workspace or last focused
      const workspaces = getWorkspaces()
      if (workspaces.length > 0 && windowManager) {
        const savedState = loadWindowState()
        const wsId = savedState?.lastFocusedWorkspaceId || workspaces[0].id
        // Verify workspace still exists
        if (workspaces.some(ws => ws.id === wsId)) {
          windowManager.createWindow({ workspaceId: wsId })
        } else {
          windowManager.createWindow({ workspaceId: workspaces[0].id })
        }
      }
    }
  })
})

app.on('window-all-closed', () => {
  // On macOS, apps typically stay active until explicitly quit
  if (process.platform !== 'darwin') {
    app.quit()
  }
})

// Track if we're in the process of quitting (to avoid re-entry)
let isQuitting = false

// Save window state and clean up resources before quitting
app.on('before-quit', async (event) => {
  // Avoid re-entry when we call app.exit()
  if (isQuitting) return
  isQuitting = true

  if (windowManager) {
    // Get full window states (includes bounds, type, and query)
    const windows = windowManager.getWindowStates()
    // Get the focused window's workspace as last focused
    const focusedWindow = BrowserWindow.getFocusedWindow()
    let lastFocusedWorkspaceId: string | undefined
    if (focusedWindow) {
      lastFocusedWorkspaceId = windowManager.getWorkspaceForWindow(focusedWindow.webContents.id) ?? undefined
    }

    saveWindowState({
      windows,
      lastFocusedWorkspaceId,
    })
    mainLog.info('Saved window state:', windows.length, 'windows')
  }

  // Flush all pending session writes before quitting
  if (sessionManager) {
    // Prevent quit until sessions are flushed
    event.preventDefault()
    try {
      await sessionManager.flushAllSessions()
      mainLog.info('Flushed all pending session writes')
    } catch (error) {
      mainLog.error('Failed to flush sessions:', error)
    }
    // Clean up SessionManager resources (file watchers, timers, etc.)
    sessionManager.cleanup()

    // If update is in progress, let electron-updater handle the quit flow
    // Force exit breaks the NSIS installer on Windows
    if (isUpdating()) {
      mainLog.info('Update in progress, letting electron-updater handle quit')
      app.quit()
      return
    }

    // Now actually quit
    app.exit(0)
  }
})

// Handle uncaught exceptions — forward to Sentry explicitly since registering
// a custom handler can interfere with @sentry/electron's automatic capture.
process.on('uncaughtException', (error) => {
  mainLog.error('Uncaught exception:', error)
  Sentry.captureException(error)
})

process.on('unhandledRejection', (reason, promise) => {
  mainLog.error('Unhandled rejection at:', promise, 'reason:', reason)
  Sentry.captureException(reason instanceof Error ? reason : new Error(String(reason)))
})
