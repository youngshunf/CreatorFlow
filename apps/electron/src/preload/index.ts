// Capture errors in the isolated preload context and forward to Sentry
import '@sentry/electron/preload'
import { contextBridge, ipcRenderer } from 'electron'
import { IPC_CHANNELS, type SessionEvent, type ElectronAPI, type FileAttachment, type AuthType, type FMDirectoryChangeEvent } from '../shared/types'

const api: ElectronAPI = {
  // Session management
  getSessions: () => ipcRenderer.invoke(IPC_CHANNELS.GET_SESSIONS),
  getSessionMessages: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.GET_SESSION_MESSAGES, sessionId),
  createSession: (workspaceId: string, options?: import('../shared/types').CreateSessionOptions) => ipcRenderer.invoke(IPC_CHANNELS.CREATE_SESSION, workspaceId, options),
  deleteSession: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.DELETE_SESSION, sessionId),
  sendMessage: (sessionId: string, message: string, attachments?: FileAttachment[], storedAttachments?: import('../shared/types').StoredAttachment[], options?: import('../shared/types').SendMessageOptions) => ipcRenderer.invoke(IPC_CHANNELS.SEND_MESSAGE, sessionId, message, attachments, storedAttachments, options),
  cancelProcessing: (sessionId: string, silent?: boolean) => ipcRenderer.invoke(IPC_CHANNELS.CANCEL_PROCESSING, sessionId, silent),
  killShell: (sessionId: string, shellId: string) => ipcRenderer.invoke(IPC_CHANNELS.KILL_SHELL, sessionId, shellId),
  getTaskOutput: (taskId: string) => ipcRenderer.invoke(IPC_CHANNELS.GET_TASK_OUTPUT, taskId),
  respondToPermission: (sessionId: string, requestId: string, allowed: boolean, alwaysAllow: boolean) =>
    ipcRenderer.invoke(IPC_CHANNELS.RESPOND_TO_PERMISSION, sessionId, requestId, allowed, alwaysAllow),
  respondToCredential: (sessionId: string, requestId: string, response: import('../shared/types').CredentialResponse) =>
    ipcRenderer.invoke(IPC_CHANNELS.RESPOND_TO_CREDENTIAL, sessionId, requestId, response),
  respondToInteractive: (sessionId: string, requestId: string, response: import('@sprouty-ai/shared/interactive-ui').InteractiveResponse) =>
    ipcRenderer.invoke(IPC_CHANNELS.RESPOND_TO_INTERACTIVE, sessionId, requestId, response),

  // Consolidated session command handler
  sessionCommand: (sessionId: string, command: import('../shared/types').SessionCommand) =>
    ipcRenderer.invoke(IPC_CHANNELS.SESSION_COMMAND, sessionId, command),

  // Pending plan execution (for reload recovery)
  getPendingPlanExecution: (sessionId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.GET_PENDING_PLAN_EXECUTION, sessionId),

  // Workspace management
  getWorkspaces: () => ipcRenderer.invoke(IPC_CHANNELS.GET_WORKSPACES),
  createWorkspace: (folderPath: string, name: string, appId?: string, appSource?: 'bundled' | 'marketplace', installMode?: 'force' | 'merge') =>
    ipcRenderer.invoke(IPC_CHANNELS.CREATE_WORKSPACE, folderPath, name, appId, appSource, installMode),
  deleteWorkspace: (workspaceId: string, mode: 'delete' | 'backup') =>
    ipcRenderer.invoke(IPC_CHANNELS.DELETE_WORKSPACE, workspaceId, mode),
  checkWorkspaceSlug: (slug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.CHECK_WORKSPACE_SLUG, slug),

  // Window management
  getWindowWorkspace: () => ipcRenderer.invoke(IPC_CHANNELS.GET_WINDOW_WORKSPACE),
  getWindowMode: () => ipcRenderer.invoke(IPC_CHANNELS.GET_WINDOW_MODE),
  openWorkspace: (workspaceId: string) => ipcRenderer.invoke(IPC_CHANNELS.OPEN_WORKSPACE, workspaceId),
  openSessionInNewWindow: (workspaceId: string, sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.OPEN_SESSION_IN_NEW_WINDOW, workspaceId, sessionId),
  switchWorkspace: (workspaceId: string) => ipcRenderer.invoke(IPC_CHANNELS.SWITCH_WORKSPACE, workspaceId),
  closeWindow: () => ipcRenderer.invoke(IPC_CHANNELS.CLOSE_WINDOW),
  confirmCloseWindow: () => ipcRenderer.invoke(IPC_CHANNELS.WINDOW_CONFIRM_CLOSE),
  onCloseRequested: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.WINDOW_CLOSE_REQUESTED, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.WINDOW_CLOSE_REQUESTED, handler)
  },
  setTrafficLightsVisible: (visible: boolean) => ipcRenderer.invoke(IPC_CHANNELS.WINDOW_SET_TRAFFIC_LIGHTS, visible),
  adjustWindowWidth: (delta: number) => ipcRenderer.invoke(IPC_CHANNELS.WINDOW_ADJUST_WIDTH, delta),

  // Event listeners
  onSessionEvent: (callback: (event: SessionEvent) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, sessionEvent: SessionEvent) => {
      callback(sessionEvent)
    }
    ipcRenderer.on(IPC_CHANNELS.SESSION_EVENT, handler)
    // Return cleanup function
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.SESSION_EVENT, handler)
    }
  },

  // File operations
  readFile: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.READ_FILE, path),
  readFileDataUrl: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.READ_FILE_DATA_URL, path),
  readFileBinary: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.READ_FILE_BINARY, path) as Promise<Uint8Array>,
  openFileDialog: () => ipcRenderer.invoke(IPC_CHANNELS.OPEN_FILE_DIALOG),
  readFileAttachment: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.READ_FILE_ATTACHMENT, path),
  storeAttachment: (sessionId: string, attachment: FileAttachment) => ipcRenderer.invoke(IPC_CHANNELS.STORE_ATTACHMENT, sessionId, attachment),
  generateThumbnail: (base64: string, mimeType: string) => ipcRenderer.invoke(IPC_CHANNELS.GENERATE_THUMBNAIL, base64, mimeType),

  // Theme
  getSystemTheme: () => ipcRenderer.invoke(IPC_CHANNELS.GET_SYSTEM_THEME),
  onSystemThemeChange: (callback: (isDark: boolean) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, isDark: boolean) => {
      callback(isDark)
    }
    ipcRenderer.on(IPC_CHANNELS.SYSTEM_THEME_CHANGED, handler)
    // Return cleanup function
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.SYSTEM_THEME_CHANGED, handler)
    }
  },

  // System
  getVersions: () => ({
    node: process.versions.node,
    chrome: process.versions.chrome,
    electron: process.versions.electron
  }),
  getHomeDir: () => ipcRenderer.invoke(IPC_CHANNELS.GET_HOME_DIR),
  isDebugMode: () => ipcRenderer.invoke(IPC_CHANNELS.IS_DEBUG_MODE),

  // Auto-update
  checkForUpdates: () => ipcRenderer.invoke(IPC_CHANNELS.UPDATE_CHECK),
  getUpdateInfo: () => ipcRenderer.invoke(IPC_CHANNELS.UPDATE_GET_INFO),
  installUpdate: () => ipcRenderer.invoke(IPC_CHANNELS.UPDATE_INSTALL),
  dismissUpdate: (version: string) => ipcRenderer.invoke(IPC_CHANNELS.UPDATE_DISMISS, version),
  getDismissedUpdateVersion: () => ipcRenderer.invoke(IPC_CHANNELS.UPDATE_GET_DISMISSED),
  onUpdateAvailable: (callback: (info: import('../shared/types').UpdateInfo) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, info: import('../shared/types').UpdateInfo) => {
      callback(info)
    }
    ipcRenderer.on(IPC_CHANNELS.UPDATE_AVAILABLE, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.UPDATE_AVAILABLE, handler)
  },
  onUpdateDownloadProgress: (callback: (progress: number) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, progress: number) => {
      callback(progress)
    }
    ipcRenderer.on(IPC_CHANNELS.UPDATE_DOWNLOAD_PROGRESS, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.UPDATE_DOWNLOAD_PROGRESS, handler)
  },

  // Shell operations
  openUrl: (url: string) => ipcRenderer.invoke(IPC_CHANNELS.OPEN_URL, url),
  openFile: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.OPEN_FILE, path),
  showInFolder: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.SHOW_IN_FOLDER, path),

  // Menu event listeners
  onMenuNewChat: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.MENU_NEW_CHAT, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.MENU_NEW_CHAT, handler)
  },
  onMenuOpenSettings: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.MENU_OPEN_SETTINGS, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.MENU_OPEN_SETTINGS, handler)
  },
  onMenuKeyboardShortcuts: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.MENU_KEYBOARD_SHORTCUTS, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.MENU_KEYBOARD_SHORTCUTS, handler)
  },
  onMenuToggleFocusMode: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.MENU_TOGGLE_FOCUS_MODE, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.MENU_TOGGLE_FOCUS_MODE, handler)
  },
  onMenuToggleSidebar: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.MENU_TOGGLE_SIDEBAR, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.MENU_TOGGLE_SIDEBAR, handler)
  },

  // Deep link navigation listener (for external sproutyai:// URLs)
  onDeepLinkNavigate: (callback: (nav: import('../shared/types').DeepLinkNavigation) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, nav: import('../shared/types').DeepLinkNavigation) => {
      callback(nav)
    }
    ipcRenderer.on(IPC_CHANNELS.DEEP_LINK_NAVIGATE, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.DEEP_LINK_NAVIGATE, handler)
  },

  // Auth
  showLogoutConfirmation: () => ipcRenderer.invoke(IPC_CHANNELS.SHOW_LOGOUT_CONFIRMATION),
  showDeleteSessionConfirmation: (name: string) => ipcRenderer.invoke(IPC_CHANNELS.SHOW_DELETE_SESSION_CONFIRMATION, name),
  logout: () => ipcRenderer.invoke(IPC_CHANNELS.LOGOUT),

  // Cloud LLM Gateway
  setCloudConfig: (config: { gatewayUrl: string; llmToken: string }) =>
    ipcRenderer.invoke(IPC_CHANNELS.CLOUD_SET_CONFIG, config),
  getCloudConfig: () =>
    ipcRenderer.invoke(IPC_CHANNELS.CLOUD_GET_CONFIG) as Promise<{ gatewayUrl: string; llmToken: string } | null>,
  clearCloudConfig: () =>
    ipcRenderer.invoke(IPC_CHANNELS.CLOUD_CLEAR_CONFIG),

  // Onboarding
  getAuthState: () => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_GET_AUTH_STATE).then(r => r.authState),
  getSetupNeeds: () => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_GET_AUTH_STATE).then(r => r.setupNeeds),
  startWorkspaceMcpOAuth: (mcpUrl: string) => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_START_MCP_OAUTH, mcpUrl),
  saveOnboardingConfig: (config: {
    authType?: AuthType
    workspace?: { name: string; iconUrl?: string; mcpUrl?: string }
    credential?: string
    mcpCredentials?: { accessToken: string; clientId?: string }
    anthropicBaseUrl?: string | null
    customModel?: string | null
  }) => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_SAVE_CONFIG, config),
  // Claude OAuth (two-step flow)
  startClaudeOAuth: () => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_START_CLAUDE_OAUTH),
  exchangeClaudeCode: (code: string) => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_EXCHANGE_CLAUDE_CODE, code),
  hasClaudeOAuthState: () => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_HAS_CLAUDE_OAUTH_STATE),
  clearClaudeOAuthState: () => ipcRenderer.invoke(IPC_CHANNELS.ONBOARDING_CLEAR_CLAUDE_OAUTH_STATE),

  // Settings - API Setup
  getApiSetup: () => ipcRenderer.invoke(IPC_CHANNELS.SETTINGS_GET_API_SETUP),
  updateApiSetup: (authType: AuthType, credential?: string, anthropicBaseUrl?: string | null, customModel?: string | null) =>
    ipcRenderer.invoke(IPC_CHANNELS.SETTINGS_UPDATE_API_SETUP, authType, credential, anthropicBaseUrl, customModel),
  testApiConnection: (apiKey: string, baseUrl?: string, modelName?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SETTINGS_TEST_API_CONNECTION, apiKey, baseUrl, modelName),

  // Settings - Model (global default)
  getModel: () => ipcRenderer.invoke(IPC_CHANNELS.SETTINGS_GET_MODEL),
  setModel: (model: string) => ipcRenderer.invoke(IPC_CHANNELS.SETTINGS_SET_MODEL, model),
  // Session-specific model (overrides global)
  getSessionModel: (sessionId: string, workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SESSION_GET_MODEL, sessionId, workspaceId),
  setSessionModel: (sessionId: string, workspaceId: string, model: string | null) =>
    ipcRenderer.invoke(IPC_CHANNELS.SESSION_SET_MODEL, sessionId, workspaceId, model),

  // Workspace Settings (per-workspace configuration)
  getWorkspaceSettings: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.WORKSPACE_SETTINGS_GET, workspaceId),
  updateWorkspaceSetting: <K extends string>(workspaceId: string, key: K, value: unknown) =>
    ipcRenderer.invoke(IPC_CHANNELS.WORKSPACE_SETTINGS_UPDATE, workspaceId, key, value),

  // Folder dialog
  openFolderDialog: (defaultPath?: string) => ipcRenderer.invoke(IPC_CHANNELS.OPEN_FOLDER_DIALOG, defaultPath),

  // Filesystem search (for @ mention file selection)
  searchFiles: (basePath: string, query: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.FS_SEARCH, basePath, query),
  // Debug: send renderer logs to main process log file
  debugLog: (...args: unknown[]) =>
    ipcRenderer.send(IPC_CHANNELS.DEBUG_LOG, ...args),

  // User Preferences
  readPreferences: () => ipcRenderer.invoke(IPC_CHANNELS.PREFERENCES_READ),
  writePreferences: (content: string) => ipcRenderer.invoke(IPC_CHANNELS.PREFERENCES_WRITE, content),

  // Session Drafts (persisted input text)
  getDraft: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.DRAFTS_GET, sessionId),
  setDraft: (sessionId: string, text: string) => ipcRenderer.invoke(IPC_CHANNELS.DRAFTS_SET, sessionId, text),
  deleteDraft: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.DRAFTS_DELETE, sessionId),
  getAllDrafts: () => ipcRenderer.invoke(IPC_CHANNELS.DRAFTS_GET_ALL),

  // Session Info Panel
  getSessionFiles: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.GET_SESSION_FILES, sessionId),
  getSessionNotes: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.GET_SESSION_NOTES, sessionId),
  setSessionNotes: (sessionId: string, content: string) => ipcRenderer.invoke(IPC_CHANNELS.SET_SESSION_NOTES, sessionId, content),
  watchSessionFiles: (sessionId: string) => ipcRenderer.invoke(IPC_CHANNELS.WATCH_SESSION_FILES, sessionId),
  unwatchSessionFiles: () => ipcRenderer.invoke(IPC_CHANNELS.UNWATCH_SESSION_FILES),
  onSessionFilesChanged: (callback: (sessionId: string) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, sessionId: string) => callback(sessionId)
    ipcRenderer.on(IPC_CHANNELS.SESSION_FILES_CHANGED, handler)
    return () => ipcRenderer.removeListener(IPC_CHANNELS.SESSION_FILES_CHANGED, handler)
  },

  // Sources
  getSources: (workspaceId: string) => ipcRenderer.invoke(IPC_CHANNELS.SOURCES_GET, workspaceId),
  createSource: (workspaceId: string, config: Partial<import('@sprouty-ai/shared/sources').FolderSourceConfig>) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_CREATE, workspaceId, config),
  deleteSource: (workspaceId: string, sourceSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_DELETE, workspaceId, sourceSlug),
  startSourceOAuth: (workspaceId: string, sourceSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_START_OAUTH, workspaceId, sourceSlug),
  saveSourceCredentials: (workspaceId: string, sourceSlug: string, credential: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_SAVE_CREDENTIALS, workspaceId, sourceSlug, credential),
  getSourcePermissionsConfig: (workspaceId: string, sourceSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_GET_PERMISSIONS, workspaceId, sourceSlug),
  getWorkspacePermissionsConfig: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.WORKSPACE_GET_PERMISSIONS, workspaceId),
  getDefaultPermissionsConfig: () =>
    ipcRenderer.invoke(IPC_CHANNELS.DEFAULT_PERMISSIONS_GET),
  // Default permissions change listener (live updates when default.json changes)
  onDefaultPermissionsChanged: (callback: () => void) => {
    const handler = () => callback()
    ipcRenderer.on(IPC_CHANNELS.DEFAULT_PERMISSIONS_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.DEFAULT_PERMISSIONS_CHANGED, handler)
    }
  },
  getMcpTools: (workspaceId: string, sourceSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SOURCES_GET_MCP_TOOLS, workspaceId, sourceSlug),

  // Session content search (full-text search via ripgrep)
  searchSessionContent: (workspaceId: string, query: string, searchId?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SEARCH_SESSIONS, workspaceId, query, searchId),

  // Status management
  listStatuses: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.STATUSES_LIST, workspaceId),
  reorderStatuses: (workspaceId: string, orderedIds: string[]) =>
    ipcRenderer.invoke(IPC_CHANNELS.STATUSES_REORDER, workspaceId, orderedIds),

  // Generic workspace image loading/saving
  readWorkspaceImage: (workspaceId: string, relativePath: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.WORKSPACE_READ_IMAGE, workspaceId, relativePath),
  writeWorkspaceImage: (workspaceId: string, relativePath: string, base64: string, mimeType: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.WORKSPACE_WRITE_IMAGE, workspaceId, relativePath, base64, mimeType),

  // Sources change listener (live updates when sources are added/removed)
  onSourcesChanged: (callback: (sources: import('@sprouty-ai/shared/sources').LoadedSource[]) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, sources: import('@sprouty-ai/shared/sources').LoadedSource[]) => {
      callback(sources)
    }
    ipcRenderer.on(IPC_CHANNELS.SOURCES_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.SOURCES_CHANGED, handler)
    }
  },

  // Skills
  getSkills: (workspaceId: string, workingDirectory?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SKILLS_GET, workspaceId, workingDirectory),
  getSkillFiles: (workspaceId: string, skillSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SKILLS_GET_FILES, workspaceId, skillSlug),
  deleteSkill: (workspaceId: string, skillSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SKILLS_DELETE, workspaceId, skillSlug),
  openSkillInEditor: (workspaceId: string, skillSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SKILLS_OPEN_EDITOR, workspaceId, skillSlug),
  openSkillInFinder: (workspaceId: string, skillSlug: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.SKILLS_OPEN_FINDER, workspaceId, skillSlug),

  // Skills change listener (live updates when skills are added/removed/modified)
  onSkillsChanged: (callback: (skills: import('@sprouty-ai/shared/skills').LoadedSkill[]) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, skills: import('@sprouty-ai/shared/skills').LoadedSkill[]) => {
      callback(skills)
    }
    ipcRenderer.on(IPC_CHANNELS.SKILLS_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.SKILLS_CHANGED, handler)
    }
  },

  // Apps (local bundled apps)
  listBundledApps: () =>
    ipcRenderer.invoke(IPC_CHANNELS.APPS_LIST_BUNDLED),

  // Marketplace
  marketplaceListSkills: (options?: import('@sprouty-ai/shared/marketplace').ListSkillsParams) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_LIST_SKILLS, options),
  marketplaceGetSkill: (skillId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_GET_SKILL, skillId),
  marketplaceListApps: (options?: import('@sprouty-ai/shared/marketplace').ListAppsParams) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_LIST_APPS, options),
  marketplaceGetApp: (appId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_GET_APP, appId),
  marketplaceGetAppSkills: (appId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_GET_APP_SKILLS, appId),
  marketplaceListCategories: () =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_LIST_CATEGORIES),
  marketplaceSearch: (query: string, options?: import('@sprouty-ai/shared/marketplace').SearchParams) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_SEARCH, query, options),
  marketplaceInstallSkill: (workspaceId: string, skillId: string, version?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_INSTALL_SKILL, workspaceId, skillId, version),
  marketplaceUpdateSkill: (workspaceId: string, skillId: string, version?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_UPDATE_SKILL, workspaceId, skillId, version),
  marketplaceCheckUpdates: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_CHECK_UPDATES, workspaceId),
  marketplaceGetInstalled: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.MARKETPLACE_GET_INSTALLED, workspaceId),
  onMarketplaceInstallProgress: (callback: (progress: import('@sprouty-ai/shared/marketplace').InstallProgress) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, progress: import('@sprouty-ai/shared/marketplace').InstallProgress) => {
      callback(progress)
    }
    ipcRenderer.on(IPC_CHANNELS.MARKETPLACE_INSTALL_PROGRESS, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.MARKETPLACE_INSTALL_PROGRESS, handler)
    }
  },

  // Statuses change listener (live updates when statuses config or icon files change)
  onStatusesChanged: (callback: (workspaceId: string) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, workspaceId: string) => {
      callback(workspaceId)
    }
    ipcRenderer.on(IPC_CHANNELS.STATUSES_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.STATUSES_CHANGED, handler)
    }
  },

  // Label management
  listLabels: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.LABELS_LIST, workspaceId),
  createLabel: (workspaceId: string, input: any) =>
    ipcRenderer.invoke(IPC_CHANNELS.LABELS_CREATE, workspaceId, input),
  deleteLabel: (workspaceId: string, labelId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.LABELS_DELETE, workspaceId, labelId),

  // Labels change listener (live updates when labels config changes)
  onLabelsChanged: (callback: (workspaceId: string) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, workspaceId: string) => {
      callback(workspaceId)
    }
    ipcRenderer.on(IPC_CHANNELS.LABELS_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.LABELS_CHANGED, handler)
    }
  },

  // Views (dynamic, expression-based filters stored in views.json)
  listViews: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.VIEWS_LIST, workspaceId),
  saveViews: (workspaceId: string, views: any[]) =>
    ipcRenderer.invoke(IPC_CHANNELS.VIEWS_SAVE, workspaceId, views),

  // Tool icon mappings (for Appearance settings page)
  getToolIconMappings: () => ipcRenderer.invoke(IPC_CHANNELS.TOOL_ICONS_GET_MAPPINGS),

  // Theme (app-level default)
  getAppTheme: () => ipcRenderer.invoke(IPC_CHANNELS.THEME_GET_APP),
  // Preset themes (app-level)
  loadPresetThemes: () => ipcRenderer.invoke(IPC_CHANNELS.THEME_GET_PRESETS),
  loadPresetTheme: (themeId: string) => ipcRenderer.invoke(IPC_CHANNELS.THEME_LOAD_PRESET, themeId),
  getColorTheme: () => ipcRenderer.invoke(IPC_CHANNELS.THEME_GET_COLOR_THEME),
  setColorTheme: (themeId: string) => ipcRenderer.invoke(IPC_CHANNELS.THEME_SET_COLOR_THEME, themeId),
  // Workspace-level theme overrides
  getWorkspaceColorTheme: (workspaceId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.THEME_GET_WORKSPACE_COLOR_THEME, workspaceId),
  setWorkspaceColorTheme: (workspaceId: string, themeId: string | null) =>
    ipcRenderer.invoke(IPC_CHANNELS.THEME_SET_WORKSPACE_COLOR_THEME, workspaceId, themeId),
  getAllWorkspaceThemes: () =>
    ipcRenderer.invoke(IPC_CHANNELS.THEME_GET_ALL_WORKSPACE_THEMES),

  // Logo URL resolution (uses Node.js filesystem cache for provider domains)
  getLogoUrl: (serviceUrl: string, provider?: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.LOGO_GET_URL, serviceUrl, provider),

  // Theme change listeners (live updates when theme.json files change)
  onAppThemeChange: (callback: (theme: import('@sprouty-ai/shared/config').ThemeOverrides | null) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, theme: import('@sprouty-ai/shared/config').ThemeOverrides | null) => {
      callback(theme)
    }
    ipcRenderer.on(IPC_CHANNELS.THEME_APP_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.THEME_APP_CHANGED, handler)
    }
  },
  // Theme preferences sync across windows (mode, colorTheme, font)
  broadcastThemePreferences: (preferences: { mode: string; colorTheme: string; font: string }) =>
    ipcRenderer.invoke(IPC_CHANNELS.THEME_BROADCAST_PREFERENCES, preferences),
  onThemePreferencesChange: (callback: (preferences: { mode: string; colorTheme: string; font: string }) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, preferences: { mode: string; colorTheme: string; font: string }) => {
      callback(preferences)
    }
    ipcRenderer.on(IPC_CHANNELS.THEME_PREFERENCES_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.THEME_PREFERENCES_CHANGED, handler)
    }
  },

  // Workspace theme sync across windows
  broadcastWorkspaceThemeChange: (workspaceId: string, themeId: string | null) =>
    ipcRenderer.invoke(IPC_CHANNELS.THEME_BROADCAST_WORKSPACE_THEME, workspaceId, themeId),
  onWorkspaceThemeChange: (callback: (data: { workspaceId: string; themeId: string | null }) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, data: { workspaceId: string; themeId: string | null }) => {
      callback(data)
    }
    ipcRenderer.on(IPC_CHANNELS.THEME_WORKSPACE_THEME_CHANGED, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.THEME_WORKSPACE_THEME_CHANGED, handler)
    }
  },

  // Notifications
  showNotification: (title: string, body: string, workspaceId: string, sessionId: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.NOTIFICATION_SHOW, title, body, workspaceId, sessionId),
  getNotificationsEnabled: () =>
    ipcRenderer.invoke(IPC_CHANNELS.NOTIFICATION_GET_ENABLED) as Promise<boolean>,
  setNotificationsEnabled: (enabled: boolean) =>
    ipcRenderer.invoke(IPC_CHANNELS.NOTIFICATION_SET_ENABLED, enabled),

  // Input settings
  getAutoCapitalisation: () =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_GET_AUTO_CAPITALISATION) as Promise<boolean>,
  setAutoCapitalisation: (enabled: boolean) =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_SET_AUTO_CAPITALISATION, enabled),
  getSendMessageKey: () =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_GET_SEND_MESSAGE_KEY) as Promise<'enter' | 'cmd-enter'>,
  setSendMessageKey: (key: 'enter' | 'cmd-enter') =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_SET_SEND_MESSAGE_KEY, key),
  getSpellCheck: () =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_GET_SPELL_CHECK) as Promise<boolean>,
  setSpellCheck: (enabled: boolean) =>
    ipcRenderer.invoke(IPC_CHANNELS.INPUT_SET_SPELL_CHECK, enabled),

  updateBadgeCount: (count: number) =>
    ipcRenderer.invoke(IPC_CHANNELS.BADGE_UPDATE, count),
  clearBadgeCount: () =>
    ipcRenderer.invoke(IPC_CHANNELS.BADGE_CLEAR),
  setDockIconWithBadge: (dataUrl: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.BADGE_SET_ICON, dataUrl),
  onBadgeDraw: (callback: (data: { count: number; iconDataUrl: string }) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, data: { count: number; iconDataUrl: string }) => {
      callback(data)
    }
    ipcRenderer.on(IPC_CHANNELS.BADGE_DRAW, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.BADGE_DRAW, handler)
    }
  },
  getWindowFocusState: () =>
    ipcRenderer.invoke(IPC_CHANNELS.WINDOW_GET_FOCUS_STATE),
  onWindowFocusChange: (callback: (isFocused: boolean) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, isFocused: boolean) => {
      callback(isFocused)
    }
    ipcRenderer.on(IPC_CHANNELS.WINDOW_FOCUS_STATE, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.WINDOW_FOCUS_STATE, handler)
    }
  },
  onNotificationNavigate: (callback: (data: { workspaceId: string; sessionId: string }) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, data: { workspaceId: string; sessionId: string }) => {
      callback(data)
    }
    ipcRenderer.on(IPC_CHANNELS.NOTIFICATION_NAVIGATE, handler)
    return () => {
      ipcRenderer.removeListener(IPC_CHANNELS.NOTIFICATION_NAVIGATE, handler)
    }
  },
  getGitBranch: (dirPath: string) =>
    ipcRenderer.invoke(IPC_CHANNELS.GET_GIT_BRANCH, dirPath),

  // Git Bash (Windows)
  checkGitBash: () => ipcRenderer.invoke(IPC_CHANNELS.GITBASH_CHECK),
  browseForGitBash: () => ipcRenderer.invoke(IPC_CHANNELS.GITBASH_BROWSE),
  setGitBashPath: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.GITBASH_SET_PATH, path),

  // File Manager
  fm: {
    listDirectory: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_LIST_DIRECTORY, path),
    createFolder: (fullPath: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_CREATE_FOLDER, fullPath),
    delete: (paths: string[]) => ipcRenderer.invoke(IPC_CHANNELS.FM_DELETE, paths),
    rename: (oldPath: string, newName: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_RENAME, oldPath, newName),
    move: (srcPaths: string[], destDir: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_MOVE, srcPaths, destDir),
    copy: (srcPaths: string[], destDir: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_COPY, srcPaths, destDir),
    getFileInfo: (path: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_GET_FILE_INFO, path),
    readFileBase64: (path: string, maxSize?: number) => ipcRenderer.invoke(IPC_CHANNELS.FM_READ_FILE_BASE64, path, maxSize),
    writeFile: (path: string, content: string) => ipcRenderer.invoke(IPC_CHANNELS.FM_WRITE_FILE, path, content),
    watchDirectory: (path: string) => ipcRenderer.send(IPC_CHANNELS.FM_WATCH_DIRECTORY, path),
    unwatchDirectory: (path: string) => ipcRenderer.send(IPC_CHANNELS.FM_UNWATCH_DIRECTORY, path),
    onDirectoryChanged: (callback: (event: FMDirectoryChangeEvent) => void) => {
      const handler = (_event: Electron.IpcRendererEvent, fmEvent: FMDirectoryChangeEvent) => {
        callback(fmEvent)
      }
      ipcRenderer.on(IPC_CHANNELS.FM_DIRECTORY_CHANGED, handler)
      return () => {
        ipcRenderer.removeListener(IPC_CHANNELS.FM_DIRECTORY_CHANGED, handler)
      }
    },
  },

  // Menu actions (for unified Craft menu)
  menuQuit: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_QUIT),
  menuNewWindow: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_NEW_WINDOW),
  menuMinimize: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_MINIMIZE),
  menuMaximize: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_MAXIMIZE),
  menuZoomIn: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_ZOOM_IN),
  menuZoomOut: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_ZOOM_OUT),
  menuZoomReset: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_ZOOM_RESET),
  menuToggleDevTools: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_TOGGLE_DEVTOOLS),
  menuUndo: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_UNDO),
  menuRedo: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_REDO),
  menuCut: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_CUT),
  menuCopy: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_COPY),
  menuPaste: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_PASTE),
  menuSelectAll: () => ipcRenderer.invoke(IPC_CHANNELS.MENU_SELECT_ALL),

  // Creator Media (自媒体创作 APP v2.0)
  creatorMedia: {
    projects: {
      list: (workspaceId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_LIST, workspaceId),
      get: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_GET, workspaceId, projectId),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_CREATE, workspaceId, data),
      update: (workspaceId: string, projectId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_UPDATE, workspaceId, projectId, data),
      delete: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_DELETE, workspaceId, projectId),
      setActive: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_SET_ACTIVE, workspaceId, projectId),
      getActive: (workspaceId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROJECTS_GET_ACTIVE, workspaceId),
    },
    profiles: {
      get: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROFILES_GET, workspaceId, projectId),
      upsert: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PROFILES_UPSERT, workspaceId, data),
    },
    platformAccounts: {
      list: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_LIST, workspaceId, projectId),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_CREATE, workspaceId, data),
      update: (workspaceId: string, id: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_UPDATE, workspaceId, id, data),
      delete: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PLATFORM_ACCOUNTS_DELETE, workspaceId, id),
    },
    competitors: {
      list: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_LIST, workspaceId, projectId),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_CREATE, workspaceId, data),
      update: (workspaceId: string, id: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_UPDATE, workspaceId, id, data),
      delete: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_COMPETITORS_DELETE, workspaceId, id),
    },
    contents: {
      list: (workspaceId: string, projectId: string, filters?: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_LIST, workspaceId, projectId, filters),
      get: (workspaceId: string, contentId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_GET, workspaceId, contentId),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_CREATE, workspaceId, data),
      update: (workspaceId: string, contentId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_UPDATE, workspaceId, contentId, data),
      updateStatus: (workspaceId: string, contentId: string, status: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_UPDATE_STATUS, workspaceId, contentId, status),
      delete: (workspaceId: string, contentId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTENTS_DELETE, workspaceId, contentId),
    },
    publishRecords: {
      list: (workspaceId: string, contentId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_LIST, workspaceId, contentId),
      get: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_GET, workspaceId, id),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_CREATE, workspaceId, data),
      update: (workspaceId: string, id: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_UPDATE, workspaceId, id, data),
      delete: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_PUBLISH_RECORDS_DELETE, workspaceId, id),
    },
    viralPatterns: {
      list: (workspaceId: string, filters?: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_LIST, workspaceId, filters),
      get: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_GET, workspaceId, id),
      create: (workspaceId: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_CREATE, workspaceId, data),
      update: (workspaceId: string, id: string, data: unknown) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_UPDATE, workspaceId, id, data),
      delete: (workspaceId: string, id: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_VIRAL_PATTERNS_DELETE, workspaceId, id),
    },
    context: {
      get: (workspaceId: string, projectId: string) => ipcRenderer.invoke(IPC_CHANNELS.CREATOR_MEDIA_CONTEXT_GET, workspaceId, projectId),
    },
  },
}

contextBridge.exposeInMainWorld('electronAPI', api)
