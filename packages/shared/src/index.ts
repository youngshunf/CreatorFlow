/**
 * @creator-flow/shared
 *
 * Shared business logic for CreatorFlow.
 * Used by the Electron app.
 *
 * Import specific modules via subpath exports:
 *   import { CreatorFlowAgent } from '@creator-flow/shared/agent';
 *   import { loadStoredConfig } from '@creator-flow/shared/config';
 *   import { getCredentialManager } from '@creator-flow/shared/credentials';
 *   import { CraftMcpClient } from '@creator-flow/shared/mcp';
 *   import { debug } from '@creator-flow/shared/utils';
 *   import { loadSource, createSource, getSourceCredentialManager } from '@creator-flow/shared/sources';
 *   import { createWorkspace, loadWorkspace } from '@creator-flow/shared/workspaces';
 *
 * Available modules:
 *   - agent: CreatorFlowAgent SDK wrapper, plan tools
 *   - auth: OAuth, token management, auth state
 *   - clients: Craft API client
 *   - config: Storage, models, preferences
 *   - credentials: Encrypted credential storage
 *   - headless: Non-interactive execution mode
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
