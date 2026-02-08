# AGENTS.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Development commands

All commands below assume the project root (`creator-flow` monorepo).

### Install and bootstrapping
- Install dependencies:
  - `bun install`

### Electron desktop app (primary UI)
- Hot-reload development (recommended when iterating on the app or shared logic):
  - `bun run electron:dev`
    - Spins up Vite for the renderer, esbuild in `--watch` mode for main and preload, and launches Electron pointed at the dev server.
- Build and run the app from a fresh state:
  - `bun run electron:start`
    - Runs the full `electron:build` pipeline (main, preload, renderer, resources) and then starts Electron on the built assets.
- Build artifacts only (no launch):
  - `bun run electron:build`
- Electron app linting:
  - From repo root (shortcut): `bun run lint:electron`
  - Equivalent from app directory: `cd apps/electron && bun run lint`

### Type checking
- Shared business-logic package only:
  - `bun run typecheck`  # runs `tsc --noEmit` in `packages/shared`
- All core packages (`@sprouty-ai/core` and `@sprouty-ai/shared`):
  - `bun run typecheck:all`

### Tests
- Monorepo-level tests (where configured for Bun):
  - `bun run test`  # invokes `bun test`
- Shared business-logic package tests (where most non-UI logic lives):
  - `cd packages/shared && bun test`
  - Watch mode for shared tests: `cd packages/shared && bun test --watch`
- Running a specific test file with Bun (useful for focused debugging):
  - `bun test path/to/file.test.ts`

### Supporting tooling
- Sync OAuth secrets into a local `.env` (required for Gmail/Slack/Microsoft integrations to work smoothly in Electron):
  - `bun run sync-secrets`
    - Requires 1Password CLI and configured vault entries as described in `apps/electron/README.md`.
- Build packaged binaries (Electron Builder) when needed:
  - Cross-platform default: `bun run electron:dist`
  - macOS only: `bun run electron:dist:mac`
  - Windows only: `bun run electron:dist:win`
  - Linux only: `bun run electron:dist:linux`
- Viewer app (log/docs viewer, if you need to work on it):
  - Dev server: `bun run viewer:dev`
  - Build: `bun run viewer:build`
  - Preview built artifacts: `bun run viewer:preview`

## Architecture overview

### Monorepo layout and roles

This repository is a Bun-based workspaces monorepo with two primary layers:

- **Apps (`apps/*`)**
  - `apps/electron`: Main desktop Electron + React UI for CreatorFlow.
    - Uses `esbuild` for the main and preload processes and Vite (React + Tailwind v4 + shadcn/ui) for the renderer.
    - Wraps the business logic in `@sprouty-ai/shared` to drive the UI via IPC.
  - `apps/viewer`: Vite-based viewer application referenced by `viewer:*` scripts (e.g., for documentation/log viewing); less central than the Electron app.

- **Packages (`packages/*`)**
  - `@sprouty-ai/core` (`packages/core`)
    - Shared TypeScript **types** and a minimal debug utility; it is intentionally light-weight.
    - Key type groups (from `CLAUDE.md` and `README.md`):
      - Workspace/auth/config (`Workspace`, `McpAuthType`, `AuthType`, `StoredConfig`, `OAuthCredentials`, `CumulativeUsage`).
      - Sessions (`Session`, `StoredSession`, `SessionMetadata`).
      - Messages and events (`Message`, `StoredMessage`, `MessageRole`, `ToolStatus`, `TokenUsage`, `AgentEvent`, `TypedError`, `Question`).
    - Design note: sessions are the primary isolation boundary; each `Session` is 1:1 with an SDK session and belongs to exactly one workspace.
  - `@sprouty-ai/shared` (`packages/shared`)
    - Core **business logic** for CreatorFlow and the main integration surface for agents:
      - `src/agent/`: `SproutyAgent`, permission modes, session-scoped tools, and permission configuration.
      - `src/auth/`: OAuth, token handling for Craft/Claude, and persisted auth state.
      - `src/config/`: Application/workspace configuration, themes, preferences, and file-watcher for live updates.
      - `src/credentials/`: AES‑256‑GCM encrypted credential storage.
      - `src/mcp/`: MCP client and connection validation.
      - `src/sessions/`: Session listing, persistence queue, and on-disk storage.
      - `src/sources/`: Definition and storage of external sources (MCP, APIs, local, gmail).
      - `src/statuses/`: Dynamic status system for session workflows.
      - `src/headless/`: Headless (non-UI) execution mode.
      - `src/prompts/`, `src/version/`, `src/workspaces/`, `src/utils/`, `src/network-interceptor.ts`: system prompt generation, version and install logic, workspace storage, shared utilities, and HTTP interception for API/MCP tooling.
    - The Electron app and any other consumers import from `@sprouty-ai/shared` rather than reaching into these directories directly.

### Electron app structure (apps/electron)

The Electron app is the main entry point for end users and is split into three layers:

- **Main process (`apps/electron/src/main`)**
  - `index.ts`: Bootstraps the app, creates windows, and wires up dev tools.
  - `sessions.ts`: Wraps `SproutyAgent`, manages session lifecycle, handles event streaming from the Claude Agent SDK, and integrates external sources.
  - `ipc.ts`: Declares IPC channels for sessions, files, shell actions, etc.
  - `menu.ts`: Application menus and accelerators.
  - `deep-link.ts`: Handles `creatorflow://…` deep links and routes them into the navigation system.
  - `agent-service.ts` / `sources-service.ts`: Discover and cache agents, validate auth, and manage sources.
  - **Critical SDK wiring:**
    - `setPathToClaudeCodeExecutable` must be pointed at `node_modules/@anthropic-ai/claude-agent-sdk/cli.js` before creating agents, because esbuild’s bundling breaks the SDK’s default `import.meta.url` resolution.
    - Auth environment variables (`CLAUDE_CODE_OAUTH_TOKEN` or `ANTHROPIC_API_KEY`) are set based on the current billing/auth state before any agents are constructed.

- **Preload (`apps/electron/src/preload`)**
  - `index.ts`: Context bridge that exposes a typed `electronAPI` surface to the renderer (navigation, session ops, file dialogs, etc.).

- **Renderer (`apps/electron/src/renderer`)**
  - React/Vite UI, using shadcn/ui and Tailwind v4.
  - Key concepts:
    - `App.tsx`: Top-level shell and wiring into IPC/state.
    - `components/chat/*`: Chat UI (input, display, attachment preview, multi-file diff view).
    - `components/ui/*`: Shared UI primitives including `source-avatar` for consistent source icons.
    - `contexts/NavigationContext.tsx` and `lib/navigate.ts`: Type-safe routing for tabs, actions, sidebars, and deep links.
    - `hooks/useAgentState.ts`, `useBackgroundTasks.ts`, `useStatuses.ts`, `useTheme.ts`: State machines and derived data for agent activation, background tasks, workflow statuses, and theme resolution.
    - `playground/`: Component playground for iterating on UI pieces in isolation.

### Configuration, storage, and runtime behavior

Understanding config layout is important when you need to adjust behavior without changing code:

- **Config and data root:** `~/.creator-flow/`
  - `config.json`: Global app configuration (workspaces list, auth types, some preferences).
  - `credentials.enc`: Encrypted credentials (AES‑256‑GCM) managed via `@sprouty-ai/shared/credentials`.
  - `preferences.json`: User preferences.
  - `theme.json`: App-level theme.
  - `workspaces/{id}/`:
    - `config.json`: Workspace-level configuration, including MCP server definitions.
    - `theme.json`: Workspace theme, merged on top of global theme.
    - `sessions/`: Session history (JSONL/JSON data used by the sessions module).
    - `sources/`: Per-source configuration and guides.
    - `skills/`: Saved agent skills.
    - `statuses/`: Dynamic status configuration used by the status system.

- **Permission modes and safety rules** (from `packages/shared/CLAUDE.md`):
  - Modes are per-session: `'safe'` (Explore), `'ask'` (Ask to Edit), `'allow-all'` (Auto).
  - Keyboard shortcut: **SHIFT+TAB** cycles modes in the UI.
  - Permission configuration is file-backed and layered:
    - Workspace-level: `~/.creator-flow/workspaces/{id}/permissions.json`.
    - Source-level: `~/.creator-flow/workspaces/{id}/sources/{slug}/permissions.json`.
  - Config fields include `blockedTools`, `allowedBashPatterns`, `allowedMcpPatterns`, `allowedApiEndpoints`, and `allowedWritePaths`.
  - When changing behavior around tool/shell access, prefer editing these JSON configs or their parsing logic rather than hard-coding special cases in the UI.

- **Dynamic statuses**
  - The workflow/status system is entirely data-driven and backed by files under `~/.creator-flow/workspaces/{id}/statuses/`.
  - Default statuses (e.g., Todo, In Progress, Needs Review, Done, Cancelled) can be customized via `createStatus`, `updateStatus`, `deleteStatus`, and `reorderStatuses` in `@sprouty-ai/shared/statuses`.

- **MCP auth separation (important invariant from core CLAUDE rules)**
  - Craft OAuth (`craft_oauth::global`) is used strictly for Craft API operations (spaces, MCP links). It is **not** reused for MCP server authentication.
  - Each MCP server has its own OAuth credentials under `workspace_oauth::{workspaceId}`.
  - When modifying or extending auth flows, preserve this separation to avoid leaking privileges across boundaries.

### Where to hook agent logic

For Warp-based automation or new features, the primary extension points are:

- `@sprouty-ai/shared/agent` and related exports for adjusting how `SproutyAgent` behaves (tools, permission hooks, summarization rules).
- `@sprouty-ai/shared/config`, `credentials`, `sources`, and `statuses` for modifying storage schemas or runtime configuration behavior.
- `apps/electron/src/main/sessions.ts` and `renderer` hooks/components for wiring new capabilities into the desktop UI.

Future agents working in this repo should prefer using these existing layers and conventions rather than introducing parallel mechanisms for configuration, permissions, or session management.

## Internationalization (i18n)

### Overview
The project uses a custom i18n system with Chinese text as translation keys. All user-facing text should be wrapped with the `t()` function.

### Usage patterns

**In `packages/ui` and `packages/shared`:**
```typescript
import { t } from '@sprouty-ai/shared/locale'

// Usage
{t('中文文本')}
```

**In `apps/electron` renderer (React components):**
```typescript
import { useT } from '../../hooks/useT'  // adjust path as needed

const MyComponent = () => {
  const t = useT()
  return <div>{t('中文文本')}</div>
}
```

### Translation files
- **Source locale (Chinese):** `apps/electron/src/renderer/locales/zh.json` - Chinese keys map to themselves
- **Target locales:** `apps/electron/src/renderer/locales/en.json`, etc. - Chinese keys map to translated values

### Adding new translations
1. Use `t('新的中文文本')` in your code
2. Run the translation script to auto-extract and translate:
   ```bash
   bun run scripts/i18n-translate.ts --translate --lang=en
   ```

### Translation script
- **Location:** `scripts/i18n-translate.ts`
- **Scanned directories:** `apps/electron/src`, `packages/shared/src`, `packages/ui/src`
- If translations are missing, check that the source directory is included in the `scanDirs` array in the script

## Git branch merge guidelines

### Merging main into creator-flow branch

**CRITICAL:** When merging the `main` branch into `creator-flow`, you MUST preserve the internationalization (i18n) changes in the `creator-flow` branch.

#### Merge strategy
1. **Before merging:** Commit all i18n changes in `creator-flow` branch first
2. **During merge:** When conflicts occur in i18n-related files, prefer the `creator-flow` version:
   - `apps/electron/src/renderer/locales/*.json` - Keep creator-flow translations
   - Components with `useT()` and `t()` calls - Keep the Chinese text wrapped with `t()`
   - `packages/shared/src/locale/*` - Keep creator-flow locale utilities
3. **After merging:** Run `bun run i18n:scan` to ensure all new text from main is captured
4. **Verification:** Check the UI to ensure all visible text displays in Chinese

#### Files requiring special attention during merge
- `apps/electron/src/renderer/locales/zh-cn.json` - Translation dictionary
- `apps/electron/src/renderer/components/app-shell/*.tsx` - UI components with i18n
- `apps/electron/src/renderer/components/ui/slash-command-menu.tsx` - Permission mode translations
- `packages/shared/src/agent/mode-types.ts` - Mode config (descriptions are translated at render time)

#### Conflict resolution priority
1. Keep i18n wrapper functions (`t()`, `useT()`) from creator-flow
2. Accept new features/logic from main
3. Merge both when possible (new feature code + i18n wrappers)
