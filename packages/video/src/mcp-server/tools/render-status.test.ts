/**
 * 渲染状态工具单元测试
 */

import { describe, expect, test, beforeAll, afterAll, beforeEach } from 'bun:test';
import { mkdirSync, rmSync, writeFileSync, existsSync, unlinkSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';
import {
  handleGetRenderStatus,
  GetRenderStatusInputSchema,
  updateRenderStatus,
  clearRenderStatus,
} from './render-status';
import { getProjectPath } from '../utils/paths';

describe('渲染状态工具', () => {
  // 创建临时测试目录
  const testWorkspace = join(tmpdir(), `render-status-test-${Date.now()}`);
  const projectName = '测试项目';
  const projectPath = getProjectPath(testWorkspace, projectName);
  const statusFilePath = join(projectPath, '.render-status.json');

  beforeAll(() => {
    // 创建测试项目目录
    mkdirSync(projectPath, { recursive: true });
  });

  afterAll(() => {
    rmSync(testWorkspace, { recursive: true, force: true });
  });

  beforeEach(() => {
    // 清理状态文件
    if (existsSync(statusFilePath)) {
      unlinkSync(statusFilePath);
    }
  });

  describe('handleGetRenderStatus', () => {
    test('项目不存在应返回错误', async () => {
      const input = {
        workspacePath: testWorkspace,
        projectName: '不存在的项目',
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(false);
      expect(response.error).toBeDefined();
    });

    test('没有状态文件应返回 idle 状态', async () => {
      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.status).toBe('idle');
      expect(response.data.progress).toBe(0);
    });

    test('应该读取渲染状态文件', async () => {
      // 创建状态文件
      const statusData = {
        status: 'rendering',
        progress: 0.5,
        currentFrame: 150,
        totalFrames: 300,
        startTime: new Date().toISOString(),
      };
      writeFileSync(statusFilePath, JSON.stringify(statusData, null, 2));

      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.status).toBe('rendering');
      expect(response.data.progress).toBe(0.5);
      expect(response.data.currentFrame).toBe(150);
      expect(response.data.totalFrames).toBe(300);
    });

    test('应该计算预估剩余时间', async () => {
      // 创建状态文件，模拟已渲染 50%，耗时 10 秒
      const startTime = new Date(Date.now() - 10000).toISOString();
      const statusData = {
        status: 'rendering',
        progress: 0.5,
        currentFrame: 150,
        totalFrames: 300,
        startTime: startTime,
      };
      writeFileSync(statusFilePath, JSON.stringify(statusData, null, 2));

      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.estimatedTimeRemaining).toBeDefined();
      expect(response.data.estimatedTimeRemaining).toBeGreaterThan(0);
      // 应该大约是 10 秒（剩余 50%）
      expect(response.data.estimatedTimeRemaining).toBeGreaterThan(8);
      expect(response.data.estimatedTimeRemaining).toBeLessThan(12);
    });

    test('应该支持按 renderId 筛选', async () => {
      // 创建状态文件
      const statusData = {
        renderId: 'render-123',
        status: 'rendering',
        progress: 0.3,
      };
      writeFileSync(statusFilePath, JSON.stringify(statusData, null, 2));

      // 查询不匹配的 renderId
      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
        renderId: 'render-456',
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.status).toBe('idle');
      expect(response.data.progress).toBe(0);
    });

    test('应该返回完整的状态信息', async () => {
      const statusData = {
        status: 'completed',
        progress: 1,
        currentFrame: 300,
        totalFrames: 300,
        outputPath: '/path/to/output.mp4',
        startTime: new Date(Date.now() - 60000).toISOString(),
        endTime: new Date().toISOString(),
      };
      writeFileSync(statusFilePath, JSON.stringify(statusData, null, 2));

      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.status).toBe('completed');
      expect(response.data.progress).toBe(1);
      expect(response.data.outputPath).toBe('/path/to/output.mp4');
      expect(response.data.startTime).toBeDefined();
      expect(response.data.endTime).toBeDefined();
    });

    test('应该处理失败状态', async () => {
      const statusData = {
        status: 'failed',
        progress: 0.7,
        error: 'Rendering failed: Out of memory',
      };
      writeFileSync(statusFilePath, JSON.stringify(statusData, null, 2));

      const input = {
        workspacePath: testWorkspace,
        projectName: projectName,
      };

      const result = await handleGetRenderStatus(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.status).toBe('failed');
      expect(response.data.error).toBe('Rendering failed: Out of memory');
    });
  });

  describe('updateRenderStatus', () => {
    test('应该创建新的状态文件', () => {
      updateRenderStatus(testWorkspace, projectName, {
        status: 'rendering',
        progress: 0.1,
        currentFrame: 30,
        totalFrames: 300,
      });

      expect(existsSync(statusFilePath)).toBe(true);

      const statusData = JSON.parse(require('fs').readFileSync(statusFilePath, 'utf-8'));
      expect(statusData.status).toBe('rendering');
      expect(statusData.progress).toBe(0.1);
      expect(statusData.currentFrame).toBe(30);
    });

    test('应该更新现有状态文件', () => {
      // 创建初始状态
      const initialStatus = {
        status: 'rendering',
        progress: 0.3,
        currentFrame: 90,
        totalFrames: 300,
      };
      writeFileSync(statusFilePath, JSON.stringify(initialStatus, null, 2));

      // 更新状态
      updateRenderStatus(testWorkspace, projectName, {
        progress: 0.6,
        currentFrame: 180,
      });

      const statusData = JSON.parse(require('fs').readFileSync(statusFilePath, 'utf-8'));
      expect(statusData.status).toBe('rendering'); // 保持不变
      expect(statusData.progress).toBe(0.6); // 已更新
      expect(statusData.currentFrame).toBe(180); // 已更新
      expect(statusData.totalFrames).toBe(300); // 保持不变
    });

    test('应该合并新旧状态', () => {
      const initialStatus = {
        status: 'rendering',
        progress: 0.5,
        startTime: '2024-01-01T00:00:00Z',
      };
      writeFileSync(statusFilePath, JSON.stringify(initialStatus, null, 2));

      updateRenderStatus(testWorkspace, projectName, {
        status: 'completed',
        progress: 1,
        endTime: '2024-01-01T00:05:00Z',
      });

      const statusData = JSON.parse(require('fs').readFileSync(statusFilePath, 'utf-8'));
      expect(statusData.status).toBe('completed');
      expect(statusData.progress).toBe(1);
      expect(statusData.startTime).toBe('2024-01-01T00:00:00Z'); // 保留
      expect(statusData.endTime).toBe('2024-01-01T00:05:00Z'); // 新增
    });

    test('项目不存在应该静默失败', () => {
      // 不应该抛出错误
      expect(() => {
        updateRenderStatus(testWorkspace, '不存在的项目', {
          status: 'rendering',
          progress: 0.5,
        });
      }).not.toThrow();
    });
  });

  describe('clearRenderStatus', () => {
    test('应该删除状态文件', () => {
      // 创建状态文件
      writeFileSync(statusFilePath, JSON.stringify({ status: 'completed' }));
      expect(existsSync(statusFilePath)).toBe(true);

      // 清除状态
      clearRenderStatus(testWorkspace, projectName);

      expect(existsSync(statusFilePath)).toBe(false);
    });

    test('状态文件不存在应该静默成功', () => {
      expect(existsSync(statusFilePath)).toBe(false);

      // 不应该抛出错误
      expect(() => {
        clearRenderStatus(testWorkspace, projectName);
      }).not.toThrow();
    });

    test('项目不存在应该静默失败', () => {
      // 不应该抛出错误
      expect(() => {
        clearRenderStatus(testWorkspace, '不存在的项目');
      }).not.toThrow();
    });
  });

  describe('GetRenderStatusInputSchema', () => {
    test('应该验证有效输入', () => {
      const input = {
        workspacePath: '/valid/path',
        projectName: '测试项目',
      };

      const result = GetRenderStatusInputSchema.safeParse(input);
      expect(result.success).toBe(true);
    });

    test('应该接受可选的 renderId', () => {
      const input = {
        workspacePath: '/valid/path',
        projectName: '测试项目',
        renderId: 'render-123',
      };

      const result = GetRenderStatusInputSchema.safeParse(input);
      expect(result.success).toBe(true);
    });

    test('应该要求 workspacePath', () => {
      const input = {
        projectName: '测试项目',
      };

      const result = GetRenderStatusInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });

    test('应该要求 projectName', () => {
      const input = {
        workspacePath: '/valid/path',
      };

      const result = GetRenderStatusInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });
  });
});
