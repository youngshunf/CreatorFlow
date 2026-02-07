/**
 * Interactive UI Types
 *
 * Defines the type system for Agent-generated interactive UI components.
 * Uses a JSON tree structure similar to json-render for LLM-friendly output.
 */

import { z } from 'zod';

// ============================================
// Option Schema (for choice components)
// ============================================

export const ChoiceOptionSchema = z.object({
  /** Unique identifier for this option */
  id: z.string(),
  /** Display label */
  label: z.string(),
  /** Optional description */
  description: z.string().optional(),
  /** Optional icon (emoji or lucide icon name) */
  icon: z.string().optional(),
  /** Whether this option is disabled */
  disabled: z.boolean().optional(),
});

export type ChoiceOption = z.infer<typeof ChoiceOptionSchema>;

// ============================================
// Form Field Schema
// ============================================

export const FormFieldTypeSchema = z.enum([
  'text',
  'textarea',
  'number',
  'email',
  'password',
  'date',
  'select',
  'checkbox',
]);

export type FormFieldType = z.infer<typeof FormFieldTypeSchema>;

export const FormFieldSchema = z.object({
  /** Unique field ID (used as key in response) */
  id: z.string(),
  /** Field type */
  type: FormFieldTypeSchema,
  /** Display label */
  label: z.string(),
  /** Placeholder text */
  placeholder: z.string().optional(),
  /** Whether field is required */
  required: z.boolean().optional(),
  /** Default value */
  defaultValue: z.union([z.string(), z.number(), z.boolean()]).optional(),
  /** Options for select type */
  options: z.array(ChoiceOptionSchema).optional(),
  /** Validation rules */
  validation: z.object({
    min: z.number().optional(),
    max: z.number().optional(),
    minLength: z.number().optional(),
    maxLength: z.number().optional(),
    pattern: z.string().optional(),
    /** Custom validation error message */
    message: z.string().optional(),
  }).optional(),
  /** Help text shown below field */
  helpText: z.string().optional(),
  /** Whether field is disabled */
  disabled: z.boolean().optional(),
  /** Hint text shown below field (alternative to helpText) */
  hint: z.string().optional(),
});

export type FormField = z.infer<typeof FormFieldSchema>;

// ============================================
// Table Column Schema
// ============================================

export const TableColumnSchema = z.object({
  /** Column ID (matches row data key) */
  id: z.string(),
  /** Column header label */
  label: z.string(),
  /** Column width (e.g., '100px', '30%', 'auto') */
  width: z.string().optional(),
  /** Text alignment */
  align: z.enum(['left', 'center', 'right']).optional(),
  /** Column data type for rendering */
  type: z.enum(['text', 'number', 'date', 'boolean', 'badge']).optional(),
});

export type TableColumn = z.infer<typeof TableColumnSchema>;

// ============================================
// Action Schema
// ============================================

export const ActionSchema = z.object({
  /** Action ID (sent to agent) */
  id: z.string(),
  /** Display label */
  label: z.string(),
  /** Button variant */
  variant: z.enum(['default', 'primary', 'secondary', 'destructive', 'ghost', 'outline']).optional(),
  /** Icon (lucide icon name) */
  icon: z.string().optional(),
});

export type Action = z.infer<typeof ActionSchema>;

// ============================================
// Component Props Schemas
// ============================================

/** SingleChoice - Radio button selection */
export const SingleChoicePropsSchema = z.object({
  /** Question/prompt text */
  label: z.string(),
  /** Available options */
  options: z.array(ChoiceOptionSchema),
  /** Currently selected value */
  value: z.string().optional(),
  /** Layout direction */
  layout: z.enum(['vertical', 'horizontal']).optional(),
  /** Whether selection is required (default: true) */
  required: z.boolean().optional(),
});

export type SingleChoiceProps = z.infer<typeof SingleChoicePropsSchema>;

/** MultiChoice - Checkbox selection */
export const MultiChoicePropsSchema = z.object({
  /** Question/prompt text */
  label: z.string(),
  /** Available options */
  options: z.array(ChoiceOptionSchema),
  /** Currently selected values */
  values: z.array(z.string()).optional(),
  /** Minimum selections required */
  min: z.number().optional(),
  /** Maximum selections allowed */
  max: z.number().optional(),
  /** Layout direction */
  layout: z.enum(['vertical', 'horizontal']).optional(),
});

export type MultiChoiceProps = z.infer<typeof MultiChoicePropsSchema>;

/** Confirm - Confirmation dialog */
export const ConfirmPropsSchema = z.object({
  /** Dialog title */
  title: z.string(),
  /** Dialog message */
  message: z.string(),
  /** Confirm button label */
  confirmLabel: z.string().optional(),
  /** Cancel button label */
  cancelLabel: z.string().optional(),
  /** Confirm button variant */
  confirmVariant: z.enum(['default', 'primary', 'destructive']).optional(),
  /** Dialog variant (affects icon and styling) */
  variant: z.enum(['default', 'info', 'warning', 'danger', 'success']).optional(),
});

export type ConfirmProps = z.infer<typeof ConfirmPropsSchema>;

/** Form - Multi-field form */
export const FormPropsSchema = z.object({
  /** Form title */
  title: z.string().optional(),
  /** Form description */
  description: z.string().optional(),
  /** Form fields */
  fields: z.array(FormFieldSchema),
  /** Submit button label */
  submitLabel: z.string().optional(),
  /** Cancel button label */
  cancelLabel: z.string().optional(),
});

export type FormProps = z.infer<typeof FormPropsSchema>;

/** Card - Information card */
export const CardActionSchema = z.object({
  /** Action ID (sent to agent) */
  id: z.string(),
  /** Display label */
  label: z.string(),
  /** Button variant */
  variant: z.enum(['default', 'primary', 'secondary', 'destructive', 'ghost', 'outline']).optional(),
  /** Icon (lucide icon name) */
  icon: z.string().optional(),
  /** Whether action is disabled */
  disabled: z.boolean().optional(),
});

export type CardAction = z.infer<typeof CardActionSchema>;

export const CardPropsSchema = z.object({
  /** Card title */
  title: z.string(),
  /** Card description/content */
  description: z.string().optional(),
  /** Image URL */
  image: z.string().optional(),
  /** Footer actions */
  actions: z.array(CardActionSchema).optional(),
  /** Card variant */
  variant: z.enum(['default', 'outlined', 'elevated']).optional(),
  /** Footer text */
  footer: z.string().optional(),
  /** Metadata items */
  metadata: z.array(z.object({
    /** Icon (lucide icon name) */
    icon: z.string().optional(),
    /** Label text */
    label: z.string().optional(),
    /** Value text */
    value: z.string(),
  })).optional(),
});

export type CardProps = z.infer<typeof CardPropsSchema>;

/** List - Ordered or unordered list */
export const ListItemSchema = z.object({
  /** Unique item ID */
  id: z.string(),
  /** Item label */
  label: z.string(),
  /** Item description */
  description: z.string().optional(),
  /** Icon (lucide icon name) */
  icon: z.string().optional(),
  /** Badge text */
  badge: z.string().optional(),
  /** Secondary text */
  secondary: z.string().optional(),
});

export type ListItem = z.infer<typeof ListItemSchema>;

export const ListPropsSchema = z.object({
  /** List title */
  title: z.string().optional(),
  /** List items */
  items: z.array(ListItemSchema),
  /** Whether list is ordered (numbered) */
  ordered: z.boolean().optional(),
  /** Whether items are selectable */
  selectable: z.boolean().optional(),
  /** Whether to use compact layout */
  compact: z.boolean().optional(),
});

export type ListProps = z.infer<typeof ListPropsSchema>;

/** DataTable - Tabular data display */
export const TableRowSchema = z.object({
  /** Unique row ID */
  id: z.string(),
  /** Row data keyed by column ID */
  data: z.record(z.string(), z.unknown()),
});

export type TableRow = z.infer<typeof TableRowSchema>;

export const DataTablePropsSchema = z.object({
  /** Table title */
  title: z.string().optional(),
  /** Column definitions */
  columns: z.array(TableColumnSchema),
  /** Row data (array of TableRow objects) */
  rows: z.array(TableRowSchema),
  /** Whether rows are selectable */
  selectable: z.boolean().optional(),
  /** Whether to allow multiple selection */
  multiSelect: z.boolean().optional(),
});

export type DataTableProps = z.infer<typeof DataTablePropsSchema>;

/** Button - Standalone action button */
export const ButtonPropsSchema = z.object({
  /** Button label */
  label: z.string(),
  /** Action ID (sent to agent) */
  action: z.string(),
  /** Button variant */
  variant: z.enum(['default', 'primary', 'secondary', 'destructive', 'ghost', 'outline']).optional(),
  /** Icon (lucide icon name) */
  icon: z.string().optional(),
  /** Whether button is full width */
  fullWidth: z.boolean().optional(),
  /** Whether button is disabled */
  disabled: z.boolean().optional(),
  /** Button size */
  size: z.enum(['sm', 'md', 'lg', 'default', 'small', 'large']).optional(),
});

export type ButtonProps = z.infer<typeof ButtonPropsSchema>;

/** ButtonGroup - Multiple buttons in a row */
export const ButtonGroupPropsSchema = z.object({
  /** Buttons to display */
  buttons: z.array(ButtonPropsSchema),
  /** Layout direction */
  layout: z.enum(['horizontal', 'vertical']).optional(),
  /** Alignment */
  align: z.enum(['start', 'center', 'end', 'stretch']).optional(),
});

export type ButtonGroupProps = z.infer<typeof ButtonGroupPropsSchema>;

// ============================================
// Component Type Union
// ============================================

export const InteractiveComponentTypeSchema = z.enum([
  'single-choice',
  'multi-choice',
  'confirm',
  'form',
  'card',
  'list',
  'data-table',
  'button',
  'button-group',
]);

export type InteractiveComponentType = z.infer<typeof InteractiveComponentTypeSchema>;

// ============================================
// Interactive UI Element (Tree Node)
// ============================================

export const InteractiveUIElementSchema: z.ZodType<InteractiveUIElement> = z.object({
  /** Unique key for this element */
  key: z.string(),
  /** Component type */
  type: InteractiveComponentTypeSchema,
  /** Component props (validated per-type at runtime) */
  props: z.record(z.string(), z.unknown()),
  /** Child elements (for container components) */
  children: z.array(z.lazy(() => InteractiveUIElementSchema)).optional(),
});

export interface InteractiveUIElement {
  key: string;
  type: InteractiveComponentType;
  props: Record<string, unknown>;
  children?: InteractiveUIElement[];
}

// ============================================
// Interactive Request (Agent → UI)
// ============================================

export const InteractiveRequestSchema = z.object({
  /** Unique request ID */
  requestId: z.string(),
  /** UI tree to render */
  tree: InteractiveUIElementSchema,
  /** Optional prompt/context shown above the UI */
  prompt: z.string().optional(),
  /** Whether user can dismiss/cancel */
  dismissible: z.boolean().optional(),
  /** Timeout in seconds (auto-cancel after timeout) */
  timeout: z.number().optional(),
});

export type InteractiveRequest = z.infer<typeof InteractiveRequestSchema>;

// ============================================
// Interactive Response (UI → Agent)
// ============================================

export const InteractiveResponseSchema = z.object({
  /** Request ID being responded to */
  requestId: z.string(),
  /** Response type */
  type: z.enum(['submit', 'cancel', 'action', 'timeout']),
  /** Response data (varies by component type) */
  data: z.unknown().optional(),
  /** Action ID (for action responses) */
  actionId: z.string().optional(),
});

export type InteractiveResponse = z.infer<typeof InteractiveResponseSchema>;

// ============================================
// Response Data Types (per component)
// ============================================

export interface SingleChoiceResponse {
  value: string;
}

export interface MultiChoiceResponse {
  values: string[];
}

export interface ConfirmResponse {
  confirmed: boolean;
}

export interface FormResponse {
  data: Record<string, unknown>;
}

export interface ListResponse {
  selectedIndex?: number;
  selectedIndices?: number[];
}

export interface DataTableResponse {
  selectedRows?: number[];
}

export interface ActionResponse {
  actionId: string;
}

// ============================================
// Interactive Status
// ============================================

export type InteractiveStatus = 'pending' | 'completed' | 'cancelled' | 'timeout';

// ============================================
// Interactive Type (for Message role)
// ============================================

export type InteractiveType = 'choice' | 'form' | 'display' | 'confirm' | 'action';

/**
 * Determine the interactive type from component type
 */
export function getInteractiveType(componentType: InteractiveComponentType): InteractiveType {
  switch (componentType) {
    case 'single-choice':
    case 'multi-choice':
      return 'choice';
    case 'confirm':
      return 'confirm';
    case 'form':
      return 'form';
    case 'button':
    case 'button-group':
      return 'action';
    case 'card':
    case 'list':
    case 'data-table':
      return 'display';
  }
}
