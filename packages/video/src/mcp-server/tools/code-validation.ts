/**
 * 代码验证工具
 *
 * MCP 工具：video_validate_composition
 *
 * 验证 Remotion 组合代码的语法和基本正确性
 */

import { z } from 'zod';
import type { FastMCP } from 'fastmcp';
import { createSuccessResponse } from '../types';
import { toErrorResponse } from '../types/errors';

// ============================================================================
// Zod Schemas
// ============================================================================

/**
 * video_validate_composition 输入 Schema
 */
export const ValidateCompositionInputSchema = z.object({
  workspacePath: z.string().describe('工作区根路径'),
  code: z.string().describe('要验证的 Remotion 组件代码'),
  props: z.record(z.string(), z.any()).optional().describe('测试用的属性'),
});

/**
 * 验证错误接口
 */
interface ValidationError {
  line: number;
  column: number;
  message: string;
  severity: 'error' | 'warning';
}

/**
 * 验证结果接口
 */
interface ValidationResult {
  valid: boolean;
  errors: ValidationError[];
  warnings: ValidationError[];
  suggestions: string[];
}

// ============================================================================
// 工具处理函数
// ============================================================================

/**
 * 验证 Remotion 组合代码
 */
async function handleValidateComposition(
  input: z.infer<typeof ValidateCompositionInputSchema>
): Promise<string> {
  try {
    const { code, props } = input;

    const result: ValidationResult = {
      valid: true,
      errors: [],
      warnings: [],
      suggestions: [],
    };

    // 1. 基本语法检查
    validateSyntax(code, result);

    // 2. Remotion API 使用检查
    validateRemotionUsage(code, result);

    // 3. React 组件结构检查
    validateReactStructure(code, result);

    // 4. 常见错误检查
    validateCommonMistakes(code, result);

    // 5. Props 使用检查
    if (props) {
      validatePropsUsage(code, props, result);
    }

    // 如果有错误，标记为无效
    if (result.errors.length > 0) {
      result.valid = false;
    }

    return JSON.stringify(createSuccessResponse(result));
  } catch (error) {
    return JSON.stringify(toErrorResponse(error));
  }
}

/**
 * 基本语法检查
 */
function validateSyntax(code: string, result: ValidationResult): void {
  // 检查括号匹配（不包括 < > 因为 JSX 标签不是简单的括号匹配）
  const brackets = { '(': ')', '[': ']', '{': '}' };
  const stack: string[] = [];
  const lines = code.split('\n');

  for (let lineNum = 0; lineNum < lines.length; lineNum++) {
    const line = lines[lineNum];
    // 跳过字符串内容（简单处理）
    let inString = false;
    let stringChar = '';

    for (let col = 0; col < line.length; col++) {
      const char = line[col];

      // 处理字符串
      if ((char === '"' || char === "'" || char === '`') && (col === 0 || line[col - 1] !== '\\')) {
        if (!inString) {
          inString = true;
          stringChar = char;
        } else if (char === stringChar) {
          inString = false;
          stringChar = '';
        }
        continue;
      }

      // 在字符串内部，跳过括号检查
      if (inString) {
        continue;
      }

      if (char in brackets) {
        stack.push(char);
      } else if (Object.values(brackets).includes(char)) {
        const last = stack.pop();
        if (!last || brackets[last as keyof typeof brackets] !== char) {
          result.errors.push({
            line: lineNum + 1,
            column: col + 1,
            message: `括号不匹配: 期望 ${last ? brackets[last as keyof typeof brackets] : '无'}, 实际 ${char}`,
            severity: 'error',
          });
        }
      }
    }
  }

  // 检查未闭合的括号
  if (stack.length > 0) {
    result.errors.push({
      line: lines.length,
      column: 1,
      message: `有 ${stack.length} 个未闭合的括号`,
      severity: 'error',
    });
  }

  // 检查常见语法错误
  if (code.includes('function(')) {
    result.warnings.push({
      line: 0,
      column: 0,
      message: '使用了 function 关键字，建议使用箭头函数',
      severity: 'warning',
    });
  }
}

/**
 * Remotion API 使用检查
 */
function validateRemotionUsage(code: string, result: ValidationResult): void {
  // 检查必需的导入（支持多行）
  const hasReactImport = /import\s+.*React.*from\s+['"]react['"]/s.test(code);
  const hasRemotionImport = /import\s+.*from\s+['"]remotion['"]/s.test(code);

  if (!hasReactImport) {
    result.errors.push({
      line: 1,
      column: 1,
      message: '缺少 React 导入: import React from "react"',
      severity: 'error',
    });
  }

  if (!hasRemotionImport) {
    result.errors.push({
      line: 1,
      column: 1,
      message: '缺少 Remotion 导入: import { ... } from "remotion"',
      severity: 'error',
    });
  }

  // 检查是否使用了 useCurrentFrame 但没有导入
  if (code.includes('useCurrentFrame') && !/import\s*\{[^}]*useCurrentFrame[^}]*\}\s*from\s+['"]remotion['"]/s.test(code)) {
    result.errors.push({
      line: 0,
      column: 0,
      message: '使用了 useCurrentFrame 但未导入',
      severity: 'error',
    });
  }

  // 检查是否使用了 useVideoConfig 但没有导入
  if (code.includes('useVideoConfig') && !/import\s*\{[^}]*useVideoConfig[^}]*\}\s*from\s+['"]remotion['"]/s.test(code)) {
    result.errors.push({
      line: 0,
      column: 0,
      message: '使用了 useVideoConfig 但未导入',
      severity: 'error',
    });
  }

  // 检查 Sequence 组件的使用
  if (code.includes('<Sequence')) {
    if (!code.includes('from=') || !code.includes('durationInFrames=')) {
      result.warnings.push({
        line: 0,
        column: 0,
        message: 'Sequence 组件应该包含 from 和 durationInFrames 属性',
        severity: 'warning',
      });
    }
  }
}

/**
 * React 组件结构检查
 */
function validateReactStructure(code: string, result: ValidationResult): void {
  // 检查是否导出了组件
  const hasExport = /export\s+(const|function|default)/.test(code);
  if (!hasExport) {
    result.errors.push({
      line: 0,
      column: 0,
      message: '组件必须被导出 (export)',
      severity: 'error',
    });
  }

  // 检查组件是否返回 JSX
  const hasReturn = /return\s*\(/.test(code) || /return\s*</.test(code);
  if (!hasReturn) {
    result.warnings.push({
      line: 0,
      column: 0,
      message: '组件应该返回 JSX 元素',
      severity: 'warning',
    });
  }

  // 检查是否使用了 AbsoluteFill
  if (!code.includes('AbsoluteFill')) {
    result.suggestions.push('建议使用 <AbsoluteFill> 作为根容器以确保全屏显示');
  }
}

/**
 * 常见错误检查
 */
function validateCommonMistakes(code: string, result: ValidationResult): void {
  // 检查是否直接使用了秒数而不是帧数
  if (/durationInFrames=\{[1-9]\d?\}/.test(code)) {
    result.warnings.push({
      line: 0,
      column: 0,
      message: 'durationInFrames 的值看起来很小，确认是否应该乘以 fps',
      severity: 'warning',
    });
  }

  // 检查是否忘记使用 frame 变量
  if (code.includes('useCurrentFrame') && !code.includes('interpolate') && !code.includes('spring')) {
    result.suggestions.push('使用了 useCurrentFrame 但没有用于动画，考虑使用 interpolate 或 spring');
  }

  // 检查样式对象
  if (code.includes('style={{') && code.includes('backgroundColor')) {
    if (!code.includes('#') && !code.includes('rgb')) {
      result.warnings.push({
        line: 0,
        column: 0,
        message: 'backgroundColor 的值可能不正确，应该是颜色值',
        severity: 'warning',
      });
    }
  }

  // 检查是否使用了中文引号
  if (code.includes('\u201c') || code.includes('\u201d') || code.includes('\u2018') || code.includes('\u2019')) {
    result.errors.push({
      line: 0,
      column: 0,
      message: '代码中包含中文引号，应该使用英文引号',
      severity: 'error',
    });
  }
}

/**
 * Props 使用检查
 */
function validatePropsUsage(
  code: string,
  props: Record<string, any>,
  result: ValidationResult
): void {
  // 检查 props 中的属性是否在代码中使用
  for (const propName of Object.keys(props)) {
    if (!code.includes(propName)) {
      result.warnings.push({
        line: 0,
        column: 0,
        message: `Props 中定义了 "${propName}" 但在代码中未使用`,
        severity: 'warning',
      });
    }
  }
}

// ============================================================================
// 工具注册
// ============================================================================

/**
 * 注册代码验证工具到 FastMCP 服务器
 */
export function registerCodeValidationTools(mcp: FastMCP): void {
  mcp.addTool({
    name: 'video_validate_composition',
    description: '验证 Remotion 组合代码的语法和正确性',
    parameters: ValidateCompositionInputSchema,
    execute: handleValidateComposition,
  });
}

// 导出处理函数供测试使用
export { handleValidateComposition };
