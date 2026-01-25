import type { PermissionRequest as PermissionRequestType, CredentialRequest as CredentialRequestType } from '../../../../shared/types'
import type { InteractiveRequest as InteractiveRequestType } from '@creator-flow/shared/interactive-ui'
import { PermissionRequest } from './structured/PermissionRequest'
import { CredentialRequest } from './structured/CredentialRequest'
import { InteractiveRequest } from './structured/InteractiveRequest'
import type { StructuredInputState, StructuredResponse } from './structured/types'

interface StructuredInputProps {
  state: StructuredInputState
  onResponse: (response: StructuredResponse) => void
  /** When true, removes container styling (shadow, bg, rounded) - used when wrapped by InputContainer */
  unstyled?: boolean
}

/**
 * StructuredInput - Router component for structured input UIs
 *
 * Routes to the appropriate component based on the input type:
 * - permission: PermissionRequest (bash command approval)
 * - credential: CredentialRequest (secure auth input)
 * - interactive: InteractiveRequest (agent interactive UI)
 */
export function StructuredInput({ state, onResponse, unstyled = false }: StructuredInputProps) {
  switch (state.type) {
    case 'permission':
      return (
        <PermissionRequest
          request={state.data as PermissionRequestType}
          onResponse={onResponse}
          unstyled={unstyled}
        />
      )
    case 'credential':
      return (
        <CredentialRequest
          request={state.data as CredentialRequestType}
          onResponse={onResponse}
          unstyled={unstyled}
        />
      )
    case 'interactive':
      return (
        <InteractiveRequest
          request={state.data as InteractiveRequestType}
          onResponse={onResponse}
          unstyled={unstyled}
        />
      )
    default:
      return null
  }
}
