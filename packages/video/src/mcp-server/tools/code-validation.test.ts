/**
 * 代码验证工具单元测试
 */

import { describe, expect, test } from 'bun:test';
import { handleValidateComposition, ValidateCompositionInputSchema } from './code-validation';

describe('代码验证工具', () => {
  const testWorkspace = '/tmp/test-workspace';

  describe('handleValidateComposition', () => {
    test('应该验证有效的 Remotion 代码', async () => {
      const validCode = `
import React from 'react';
import { AbsoluteFill, useCurrentFrame, useVideoConfig } from 'remotion';

export const MyComposition: React.FC = () => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();

  return (
    <AbsoluteFill style={{ backgroundColor: '#000' }}>
      <h1>Frame: {frame}</h1>
    </AbsoluteFill>
  );
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: validCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(true);
      expect(response.data.errors.length).toBe(0);
    });

    test('应该检测缺少 React 导入', async () => {
      const invalidCode = `
import { AbsoluteFill } from 'remotion';

export const MyComposition = () => {
  return <AbsoluteFill />;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('React'))).toBe(true);
    });

    test('应该检测缺少 Remotion 导入', async () => {
      const invalidCode = `
import React from 'react';

export const MyComposition: React.FC = () => {
  return <div>Hello</div>;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('Remotion'))).toBe(true);
    });

    test('应该检测括号不匹配', async () => {
      const invalidCode = `
import React from 'react';
import { AbsoluteFill } from 'remotion';

export const MyComposition: React.FC = () => {
  return (
    <AbsoluteFill>
      <div>
        <h1>Title</h1>
      </div>
    </AbsoluteFill>
  );
  // 缺少一个闭合括号
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('括号'))).toBe(true);
    });

    test('应该检测使用了 useCurrentFrame 但未导入', async () => {
      const invalidCode = `
import React from 'react';
import { AbsoluteFill } from 'remotion';

export const MyComposition: React.FC = () => {
  const frame = useCurrentFrame();
  return <AbsoluteFill>{frame}</AbsoluteFill>;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('useCurrentFrame'))).toBe(true);
    });

    test('应该检测使用了 useVideoConfig 但未导入', async () => {
      const invalidCode = `
import React from 'react';
import { AbsoluteFill } from 'remotion';

export const MyComposition: React.FC = () => {
  const { fps } = useVideoConfig();
  return <AbsoluteFill>{fps}</AbsoluteFill>;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('useVideoConfig'))).toBe(true);
    });

    test('应该检测缺少组件导出', async () => {
      const invalidCode = `
import React from 'react';
import { AbsoluteFill } from 'remotion';

const MyComposition: React.FC = () => {
  return <AbsoluteFill />;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('导出'))).toBe(true);
    });

    test('应该检测中文引号', async () => {
      // 使用 Unicode 转义序列确保是中文引号
      const invalidCode = `
import React from 'react';
import { AbsoluteFill } from "remotion";

export const MyComposition: React.FC = () => {
  const title = \u201c测试标题\u201d;
  return <AbsoluteFill style={{ backgroundColor: "#000" }} />;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.valid).toBe(false);
      expect(response.data.errors.some((e: any) => e.message.includes('中文引号'))).toBe(true);
    });

    test('应该提供建议使用 AbsoluteFill', async () => {
      const codeWithoutAbsoluteFill = `
import React from 'react';
import { useCurrentFrame } from 'remotion';

export const MyComposition: React.FC = () => {
  const frame = useCurrentFrame();
  return <div>{frame}</div>;
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: codeWithoutAbsoluteFill,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.suggestions.some((s: string) => s.includes('AbsoluteFill'))).toBe(true);
    });

    test('应该检测 Sequence 缺少必需属性', async () => {
      const invalidCode = `
import React from 'react';
import { AbsoluteFill, Sequence } from 'remotion';

export const MyComposition: React.FC = () => {
  return (
    <AbsoluteFill>
      <Sequence>
        <div>Scene 1</div>
      </Sequence>
    </AbsoluteFill>
  );
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: invalidCode,
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.warnings.some((w: any) => w.message.includes('Sequence'))).toBe(true);
    });

    test('应该检测 Props 使用', async () => {
      const code = `
import React from 'react';
import { AbsoluteFill } from 'remotion';

interface Props {
  title: string;
}

export const MyComposition: React.FC<Props> = ({ title }) => {
  return (
    <AbsoluteFill>
      <h1>{title}</h1>
    </AbsoluteFill>
  );
};
`;

      const input = {
        workspacePath: testWorkspace,
        code: code,
        props: {
          title: 'Test Title',
          unusedProp: 'Test Value', // 这个 prop 在代码中完全没有出现
        },
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.warnings.some((w: any) => w.message.includes('unusedProp'))).toBe(true);
    });

    test('应该返回完整的验证结果结构', async () => {
      const input = {
        workspacePath: testWorkspace,
        code: 'import React from "react";',
      };

      const result = await handleValidateComposition(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data).toHaveProperty('valid');
      expect(response.data).toHaveProperty('errors');
      expect(response.data).toHaveProperty('warnings');
      expect(response.data).toHaveProperty('suggestions');
      expect(Array.isArray(response.data.errors)).toBe(true);
      expect(Array.isArray(response.data.warnings)).toBe(true);
      expect(Array.isArray(response.data.suggestions)).toBe(true);
    });
  });

  describe('ValidateCompositionInputSchema', () => {
    test('应该验证有效输入', () => {
      const input = {
        workspacePath: '/valid/path',
        code: 'const x = 1;',
      };

      const result = ValidateCompositionInputSchema.safeParse(input);
      expect(result.success).toBe(true);
    });

    test('应该接受可选的 props', () => {
      const input = {
        workspacePath: '/valid/path',
        code: 'const x = 1;',
        props: {
          title: 'Test',
          count: 42,
        },
      };

      const result = ValidateCompositionInputSchema.safeParse(input);
      expect(result.success).toBe(true);
    });

    test('应该要求 workspacePath', () => {
      const input = {
        code: 'const x = 1;',
      };

      const result = ValidateCompositionInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });

    test('应该要求 code', () => {
      const input = {
        workspacePath: '/valid/path',
      };

      const result = ValidateCompositionInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });
  });
});
