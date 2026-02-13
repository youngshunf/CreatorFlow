/**
 * @sprouty-ai/shared
 *
 * Shared business logic for Sprouty AI.
 * Used by the Electron app.
 *
 * Import specific modules via subpath exports:
 *   import { SproutyAgent } from '@sprouty-ai/shared/agent';
 *   import { loadStoredConfig } from '@sprouty-ai/shared/config';
 *   import { getCredentialManager } from '@sprouty-ai/shared/credentials';
 *   import { CraftMcpClient } from '@sprouty-ai/shared/mcp';
 *   import { debug } from '@sprouty-ai/shared/utils';
 *   import { loadSource, createSource, getSourceCredentialManager } from '@sprouty-ai/shared/sources';
 *   import { createWorkspace, loadWorkspace } from '@sprouty-ai/shared/workspaces';
 *
 * Available modules:
 *   - agent: SproutyAgent SDK wrapper, plan tools
 *   - auth: OAuth, token management, auth state
 *   - clients: Craft API client
 *   - config: Storage, models, preferences
 *   - credentials: Encrypted credential storage
 *   - mcp: MCP client, connection validation
 *   - prompts: System prompt generation
 *   - sources: Workspace-scoped source management (MCP, API, local)
 *   - utils: Debug logging, file handling, summarization
 *   - validation: URL validation
 *   - version: Version and installation management
 *   - workspaces: Workspace management (top-level organizational unit)
 */

// Export branding (standalone, no dependencies)
export * from './branding.ts';
