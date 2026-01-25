/**
 * Interactive UI Helpers
 *
 * Factory functions for creating interactive UI elements.
 * These make it easy for Agent code to construct valid UI trees.
 */

import { randomUUID } from 'crypto';
import type {
  InteractiveUIElement,
  InteractiveRequest,
  SingleChoiceProps,
  MultiChoiceProps,
  ConfirmProps,
  FormProps,
  CardProps,
  ListProps,
  DataTableProps,
  ButtonProps,
  ButtonGroupProps,
  ChoiceOption,
  FormField,
  TableColumn,
  Action,
} from './types';

// ============================================
// Element Factory Functions
// ============================================

/**
 * Create a SingleChoice element
 */
export function createSingleChoice(props: SingleChoiceProps): InteractiveUIElement {
  return {
    key: `single-choice-${randomUUID().slice(0, 8)}`,
    type: 'single-choice',
    props,
  };
}

/**
 * Create a MultiChoice element
 */
export function createMultiChoice(props: MultiChoiceProps): InteractiveUIElement {
  return {
    key: `multi-choice-${randomUUID().slice(0, 8)}`,
    type: 'multi-choice',
    props,
  };
}

/**
 * Create a Confirm element
 */
export function createConfirm(props: ConfirmProps): InteractiveUIElement {
  return {
    key: `confirm-${randomUUID().slice(0, 8)}`,
    type: 'confirm',
    props,
  };
}

/**
 * Create a Form element
 */
export function createForm(props: FormProps): InteractiveUIElement {
  return {
    key: `form-${randomUUID().slice(0, 8)}`,
    type: 'form',
    props,
  };
}

/**
 * Create a Card element
 */
export function createCard(props: CardProps): InteractiveUIElement {
  return {
    key: `card-${randomUUID().slice(0, 8)}`,
    type: 'card',
    props,
  };
}

/**
 * Create a List element
 */
export function createList(props: ListProps): InteractiveUIElement {
  return {
    key: `list-${randomUUID().slice(0, 8)}`,
    type: 'list',
    props,
  };
}

/**
 * Create a DataTable element
 */
export function createDataTable(props: DataTableProps): InteractiveUIElement {
  return {
    key: `data-table-${randomUUID().slice(0, 8)}`,
    type: 'data-table',
    props,
  };
}

/**
 * Create a Button element
 */
export function createButton(props: ButtonProps): InteractiveUIElement {
  return {
    key: `button-${randomUUID().slice(0, 8)}`,
    type: 'button',
    props,
  };
}

/**
 * Create a ButtonGroup element
 */
export function createButtonGroup(props: ButtonGroupProps): InteractiveUIElement {
  return {
    key: `button-group-${randomUUID().slice(0, 8)}`,
    type: 'button-group',
    props,
  };
}

// ============================================
// Request Factory Functions
// ============================================

/**
 * Create an InteractiveRequest
 */
export function createInteractiveRequest(
  tree: InteractiveUIElement,
  options?: {
    prompt?: string;
    dismissible?: boolean;
    timeout?: number;
  }
): InteractiveRequest {
  return {
    requestId: `interactive-${randomUUID()}`,
    tree,
    prompt: options?.prompt,
    dismissible: options?.dismissible ?? true,
    timeout: options?.timeout,
  };
}

// ============================================
// Option/Field Builder Helpers
// ============================================

/**
 * Create a choice option
 */
export function option(id: string, label: string, description?: string, icon?: string): ChoiceOption {
  return { id, label, description, icon };
}

/**
 * Create multiple choice options from a simple array
 */
export function optionsFromArray(items: string[]): ChoiceOption[] {
  return items.map((item, index) => ({
    id: `option-${index}`,
    label: item,
  }));
}

/**
 * Create a form field
 */
export function field(
  id: string,
  type: FormField['type'],
  label: string,
  options?: Partial<Omit<FormField, 'id' | 'type' | 'label'>>
): FormField {
  return { id, type, label, ...options };
}

/**
 * Create a table column
 */
export function column(
  id: string,
  label: string,
  options?: Partial<Omit<TableColumn, 'id' | 'label'>>
): TableColumn {
  return { id, label, ...options };
}

/**
 * Create an action
 */
export function action(
  id: string,
  label: string,
  options?: Partial<Omit<Action, 'id' | 'label'>>
): Action {
  return { id, label, ...options };
}

// ============================================
// Common UI Patterns (Shorthand)
// ============================================

/**
 * Quick yes/no confirmation
 */
export function yesNoConfirm(title: string, message: string): InteractiveRequest {
  return createInteractiveRequest(
    createConfirm({
      title,
      message,
      confirmLabel: 'Yes',
      cancelLabel: 'No',
    })
  );
}

/**
 * Quick single choice from array of strings
 */
export function quickChoice(label: string, options: string[]): InteractiveRequest {
  return createInteractiveRequest(
    createSingleChoice({
      label,
      options: optionsFromArray(options),
    })
  );
}

/**
 * Quick multi choice from array of strings
 */
export function quickMultiChoice(label: string, options: string[], min?: number, max?: number): InteractiveRequest {
  return createInteractiveRequest(
    createMultiChoice({
      label,
      options: optionsFromArray(options),
      min,
      max,
    })
  );
}

/**
 * Create a multi-question form from multiple interactive elements.
 * All questions are rendered as a single form with one submit button.
 *
 * @example
 * ```ts
 * const request = multiQuestionForm(
 *   [
 *     createSingleChoice({ label: 'Preferred language?', options: optionsFromArray(['TypeScript', 'JavaScript', 'Python']) }),
 *     createMultiChoice({ label: 'Which frameworks?', options: optionsFromArray(['React', 'Vue', 'Angular']), min: 1 }),
 *     createConfirm({ title: 'Use strict mode?', message: 'Enable strict TypeScript mode?' }),
 *   ],
 *   { prompt: 'Please answer these questions to configure your project:' }
 * );
 * ```
 */
export function multiQuestionForm(
  elements: InteractiveUIElement[],
  options?: {
    prompt?: string;
    dismissible?: boolean;
    timeout?: number;
  }
): InteractiveRequest {
  // Create a container element with children
  const container: InteractiveUIElement = {
    key: `form-container-${randomUUID().slice(0, 8)}`,
    type: 'form', // Container type
    props: {}, // No props needed for container
    children: elements,
  };

  return {
    requestId: `interactive-${randomUUID()}`,
    tree: container,
    prompt: options?.prompt,
    dismissible: options?.dismissible ?? true,
    timeout: options?.timeout,
  };
}

/**
 * Simple text input form
 */
export function textInputForm(
  label: string,
  placeholder?: string,
  options?: { required?: boolean; multiline?: boolean }
): InteractiveRequest {
  return createInteractiveRequest(
    createForm({
      fields: [
        field('value', options?.multiline ? 'textarea' : 'text', label, {
          placeholder,
          required: options?.required,
        }),
      ],
      submitLabel: 'Submit',
    })
  );
}

/**
 * Action buttons row
 */
export function actionButtons(...buttons: Array<{ id: string; label: string; variant?: ButtonProps['variant'] }>): InteractiveRequest {
  return createInteractiveRequest(
    createButtonGroup({
      buttons: buttons.map(b => ({
        label: b.label,
        action: b.id,
        variant: b.variant,
      })),
      layout: 'horizontal',
    })
  );
}

// ============================================
// Validation Helpers
// ============================================

/**
 * Check if a response is a submit response
 */
export function isSubmitResponse(response: { type: string }): boolean {
  return response.type === 'submit';
}

/**
 * Check if a response is a cancel response
 */
export function isCancelResponse(response: { type: string }): boolean {
  return response.type === 'cancel';
}

/**
 * Check if a response is an action response
 */
export function isActionResponse(response: { type: string }): boolean {
  return response.type === 'action';
}

/**
 * Extract value from SingleChoice response
 */
export function getSingleChoiceValue(response: { data?: unknown }): string | undefined {
  const data = response.data as { value?: string } | undefined;
  return data?.value;
}

/**
 * Extract values from MultiChoice response
 */
export function getMultiChoiceValues(response: { data?: unknown }): string[] {
  const data = response.data as { values?: string[] } | undefined;
  return data?.values ?? [];
}

/**
 * Extract confirmed from Confirm response
 */
export function getConfirmResult(response: { data?: unknown }): boolean {
  const data = response.data as { confirmed?: boolean } | undefined;
  return data?.confirmed ?? false;
}

/**
 * Extract form data from Form response
 */
export function getFormData<T extends Record<string, unknown>>(response: { data?: unknown }): T | undefined {
  const data = response.data as { data?: T } | undefined;
  return data?.data;
}

/**
 * Extract action ID from action response
 */
export function getActionId(response: { actionId?: string }): string | undefined {
  return response.actionId;
}
