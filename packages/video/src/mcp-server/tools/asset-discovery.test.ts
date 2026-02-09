/**
 * 素材发现工具单元测试
 */

import { describe, expect, test, beforeAll, afterAll } from 'bun:test';
import { mkdirSync, rmSync, writeFileSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';
import { handleListAvailableAssets, ListAvailableAssetsInputSchema } from './asset-discovery';

describe('素材发现工具', () => {
  // 创建临时测试目录
  const testWorkspace = join(tmpdir(), `asset-discovery-test-${Date.now()}`);
  const assetsDir = join(testWorkspace, 'assets');
  const imagesDir = join(assetsDir, 'images');
  const videosDir = join(assetsDir, 'videos');
  const audioDir = join(assetsDir, 'audio');
  const fontsDir = join(assetsDir, 'fonts');

  beforeAll(() => {
    // 创建测试目录结构
    mkdirSync(imagesDir, { recursive: true });
    mkdirSync(videosDir, { recursive: true });
    mkdirSync(audioDir, { recursive: true });
    mkdirSync(fontsDir, { recursive: true });

    // 创建测试文件
    writeFileSync(join(imagesDir, 'test1.png'), 'fake png content');
    writeFileSync(join(imagesDir, 'test2.jpg'), 'fake jpg content');
    writeFileSync(join(videosDir, 'video1.mp4'), 'fake mp4 content');
    writeFileSync(join(audioDir, 'audio1.mp3'), 'fake mp3 content');
    writeFileSync(join(fontsDir, 'font1.ttf'), 'fake ttf content');

    // 创建一些隐藏文件（应该被忽略）
    writeFileSync(join(imagesDir, '.hidden.png'), 'hidden file');

    // 创建 node_modules 目录（应该被忽略）
    const nodeModulesDir = join(testWorkspace, 'node_modules');
    mkdirSync(nodeModulesDir, { recursive: true });
    writeFileSync(join(nodeModulesDir, 'test.png'), 'should be ignored');
  });

  afterAll(() => {
    rmSync(testWorkspace, { recursive: true, force: true });
  });

  describe('handleListAvailableAssets', () => {
    test('应该列出所有素材', async () => {
      const input = {
        workspacePath: testWorkspace,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.assets.length).toBeGreaterThan(0);
      expect(response.data.total).toBeGreaterThan(0);
    });

    test('应该按类型筛选图片', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'image' as const,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.length).toBe(2); // test1.png, test2.jpg
      expect(assets.every((a: any) => a.type === 'image')).toBe(true);
    });

    test('应该按类型筛选视频', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'video' as const,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.length).toBe(1); // video1.mp4
      expect(assets[0].type).toBe('video');
    });

    test('应该按类型筛选音频', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'audio' as const,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.length).toBe(1); // audio1.mp3
      expect(assets[0].type).toBe('audio');
    });

    test('应该按类型筛选字体', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'font' as const,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.length).toBe(1); // font1.ttf
      expect(assets[0].type).toBe('font');
    });

    test('应该支持搜索模式', async () => {
      const input = {
        workspacePath: testWorkspace,
        searchPattern: '*.png',
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.length).toBe(1); // test1.png
      expect(assets[0].name).toBe('test1.png');
    });

    test('应该限制返回数量', async () => {
      const input = {
        workspacePath: testWorkspace,
        maxResults: 2,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      expect(response.data.returned).toBeLessThanOrEqual(2);
    });

    test('应该忽略隐藏文件', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'image' as const,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.every((a: any) => !a.name.startsWith('.'))).toBe(true);
    });

    test('应该忽略 node_modules 目录', async () => {
      const input = {
        workspacePath: testWorkspace,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const assets = response.data.assets;
      expect(assets.every((a: any) => !a.path.includes('node_modules'))).toBe(true);
    });

    test('工作区不存在应返回错误', async () => {
      const input = {
        workspacePath: '/nonexistent/path',
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(false);
      expect(response.error).toBeDefined();
    });

    test('应该返回素材的完整信息', async () => {
      const input = {
        workspacePath: testWorkspace,
        assetType: 'image' as const,
        maxResults: 1,
      };

      const result = await handleListAvailableAssets(input);
      const response = JSON.parse(result);

      expect(response.success).toBe(true);
      const asset = response.data.assets[0];
      expect(asset).toHaveProperty('path');
      expect(asset).toHaveProperty('name');
      expect(asset).toHaveProperty('type');
      expect(asset).toHaveProperty('size');
      expect(asset.size).toBeGreaterThan(0);
    });
  });

  describe('ListAvailableAssetsInputSchema', () => {
    test('应该验证有效输入', () => {
      const input = {
        workspacePath: '/valid/path',
        assetType: 'image',
        searchPattern: '*.png',
        maxResults: 50,
      };

      const result = ListAvailableAssetsInputSchema.safeParse(input);
      expect(result.success).toBe(true);
    });

    test('应该使用默认值', () => {
      const input = {
        workspacePath: '/valid/path',
      };

      const result = ListAvailableAssetsInputSchema.safeParse(input);
      expect(result.success).toBe(true);
      if (result.success) {
        expect(result.data.maxResults).toBe(100);
      }
    });

    test('应该拒绝无效的素材类型', () => {
      const input = {
        workspacePath: '/valid/path',
        assetType: 'invalid',
      };

      const result = ListAvailableAssetsInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });

    test('应该要求 workspacePath', () => {
      const input = {};

      const result = ListAvailableAssetsInputSchema.safeParse(input);
      expect(result.success).toBe(false);
    });
  });
});
