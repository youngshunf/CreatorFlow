import type { PermissionRequest, CredentialRequest, CredentialResponse } from '../../../../../shared/types'
import type { InteractiveRequest, InteractiveResponse } from '@sprouty-ai/shared/interactive-ui'

/**
 * Input mode determines which component is rendered in InputContainer
 */
export type InputMode = 'freeform' | 'structured'

/**
 * Types of structured input UIs
 */
export type StructuredInputType = 'permission' | 'credential' | 'interactive'

/**
 * Union type for structured input data
 */
export type StructuredInputData =
  | { type: 'permission'; data: PermissionRequest }
  | { type: 'credential'; data: CredentialRequest }
  | { type: 'interactive'; data: InteractiveRequest }

/**
 * State for structured input
 */
export interface StructuredInputState {
  type: StructuredInputType
  data: PermissionRequest | CredentialRequest | InteractiveRequest
}

/**
 * Response from permission request
 */
export interface PermissionResponse {
  type: 'permission'
  allowed: boolean
  alwaysAllow: boolean
}

/**
 * Union type for all structured responses
 */
export type StructuredResponse = PermissionResponse | CredentialResponse | InteractiveResponse

// Re-export for convenience
export type { CredentialResponse }
export type { InteractiveRequest, InteractiveResponse }
