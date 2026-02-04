/**
 * 验证器单元测试
 */

import { describe, expect, test, beforeAll, afterAll } from 'bun:test';
import { mkdirSync, rmSync, writeFileSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';
import {
  validateProjectName,
  validateWorkspacePath,
  validateAssetType,
  validateAssetExtension,
  inferAssetType,
  validateAssetFile,
  isValidUUID,
  validateProjectId,
  validateCompositionId,
  validateAssetId,
} from './validators';

describe('验证器', () => {
  // 创建临时测试目录
  const testDir = join(tmpdir(), `mcp-video-test-${Date.now()}`);
  const testFile = join(testDir, 'test-image.png');

  beforeAll(() => {
    mkdirSync(testDir, { recursive: true });
    writeFileSync(testFile, 'fake image content');
  });

  afterAll(() => {
    rmSync(testDir, { recursive: true, force: true });
  });

  describe('validateProjectName', () => {
    test('有效的项目名称应通过验证', () => {
      const result = validateProjectName('我的视频项目');
      expect(result.valid).toBe(true);
    });

    test('英文项目名称应通过验证', () => {
      const result = validateProjectName('My Video Project');
      expect(result.valid).toBe(true);
    });

    test('空字符串应验证失败', () => {
      const result = validateProjectName('');
      expect(result.valid).toBe(false);
    });

    test('仅空白字符应验证失败', () => {
      const result = validateProjectName('   ');
      expect(result.valid).toBe(false);
    });

    test('包含无效字符应验证失败', () => {
      const result = validateProjectName('project<name>');
      expect(result.valid).toBe(false);
    });

    test('以点号开头应验证失败', () => {
      const result = validateProjectName('.hidden');
      expect(result.valid).toBe(false);
    });

    test('以点号结尾应验证失败', () => {
      const result = validateProjectName('project.');
      expect(result.valid).toBe(false);
    });

    test('Windows 保留名称应验证失败', () => {
      const result = validateProjectName('CON');
      expect(result.valid).toBe(false);
    });

    test('超长名称应验证失败', () => {
      const longName = 'a'.repeat(101);
      const result = validateProjectName(longName);
      expect(result.valid).toBe(false);
    });
  });

  describe('validateWorkspacePath', () => {
    test('有效的绝对路径应通过验证', () => {
      const result = validateWorkspacePath(testDir);
      expect(result.valid).toBe(true);
    });

    test('空字符串应验证失败', () => {
      const result = validateWorkspacePath('');
      expect(result.valid).toBe(false);
    });

    test('相对路径应验证失败', () => {
      const result = validateWorkspacePath('./relative/path');
      expect(result.valid).toBe(false);
    });

    test('不存在的路径应验证失败', () => {
      const result = validateWorkspacePath('/nonexistent/path/12345');
      expect(result.valid).toBe(false);
    });
  });

  describe('validateAssetType', () => {
    test('image 类型应通过验证', () => {
      const result = validateAssetType('image');
      expect(result.valid).toBe(true);
      if (result.valid) {
        expect(result.type).toBe('image');
      }
    });

    test('video 类型应通过验证', () => {
      const result = validateAssetType('video');
      expect(result.valid).toBe(true);
    });

    test('audio 类型应通过验证', () => {
      const result = validateAssetType('audio');
      expect(result.valid).toBe(true);
    });

    test('font 类型应通过验证', () => {
      const result = validateAssetType('font');
      expect(result.valid).toBe(true);
    });

    test('无效类型应验证失败', () => {
      const result = validateAssetType('unknown');
      expect(result.valid).toBe(false);
    });
  });

  describe('validateAssetExtension', () => {
    test('PNG 图片扩展名应通过验证', () => {
      const result = validateAssetExtension('/path/to/image.png', 'image');
      expect(result.valid).toBe(true);
      if (result.valid) {
        expect(result.extension).toBe('.png');
      }
    });

    test('JPG 图片扩展名应通过验证', () => {
      const result = validateAssetExtension('/path/to/image.jpg', 'image');
      expect(result.valid).toBe(true);
    });

    test('MP4 视频扩展名应通过验证', () => {
      const result = validateAssetExtension('/path/to/video.mp4', 'video');
      expect(result.valid).toBe(true);
    });

    test('MP3 音频扩展名应通过验证', () => {
      const result = validateAssetExtension('/path/to/audio.mp3', 'audio');
      expect(result.valid).toBe(true);
    });

    test('TTF 字体扩展名应通过验证', () => {
      const result = validateAssetExtension('/path/to/font.ttf', 'font');
      expect(result.valid).toBe(true);
    });

    test('不支持的扩展名应验证失败', () => {
      const result = validateAssetExtension('/path/to/file.txt', 'image');
      expect(result.valid).toBe(false);
    });

    test('类型不匹配应验证失败', () => {
      const result = validateAssetExtension('/path/to/image.png', 'video');
      expect(result.valid).toBe(false);
    });
  });

  describe('inferAssetType', () => {
    test('应正确推断图片类型', () => {
      expect(inferAssetType('/path/to/image.png')).toBe('image');
      expect(inferAssetType('/path/to/image.jpg')).toBe('image');
      expect(inferAssetType('/path/to/image.jpeg')).toBe('image');
      expect(inferAssetType('/path/to/image.gif')).toBe('image');
      expect(inferAssetType('/path/to/image.webp')).toBe('image');
      expect(inferAssetType('/path/to/image.svg')).toBe('image');
    });

    test('应正确推断视频类型', () => {
      expect(inferAssetType('/path/to/video.mp4')).toBe('video');
      expect(inferAssetType('/path/to/video.webm')).toBe('video');
      expect(inferAssetType('/path/to/video.mov')).toBe('video');
    });

    test('应正确推断音频类型', () => {
      expect(inferAssetType('/path/to/audio.mp3')).toBe('audio');
      expect(inferAssetType('/path/to/audio.wav')).toBe('audio');
      expect(inferAssetType('/path/to/audio.ogg')).toBe('audio');
      expect(inferAssetType('/path/to/audio.m4a')).toBe('audio');
    });

    test('应正确推断字体类型', () => {
      expect(inferAssetType('/path/to/font.ttf')).toBe('font');
      expect(inferAssetType('/path/to/font.otf')).toBe('font');
      expect(inferAssetType('/path/to/font.woff')).toBe('font');
      expect(inferAssetType('/path/to/font.woff2')).toBe('font');
    });

    test('未知扩展名应返回 null', () => {
      expect(inferAssetType('/path/to/file.txt')).toBeNull();
      expect(inferAssetType('/path/to/file.xyz')).toBeNull();
    });
  });

  describe('validateAssetFile', () => {
    test('存在的有效文件应通过验证', () => {
      const result = validateAssetFile(testFile, 'image');
      expect(result.valid).toBe(true);
    });

    test('不存在的文件应验证失败', () => {
      const result = validateAssetFile('/nonexistent/file.png', 'image');
      expect(result.valid).toBe(false);
    });

    test('空路径应验证失败', () => {
      const result = validateAssetFile('', 'image');
      expect(result.valid).toBe(false);
    });
  });

  describe('isValidUUID', () => {
    test('有效的 UUID v4 应返回 true', () => {
      expect(isValidUUID('550e8400-e29b-41d4-a716-446655440000')).toBe(true);
      expect(isValidUUID('6ba7b810-9dad-41d4-80b4-00c04fd430c8')).toBe(true);
    });

    test('无效的 UUID 应返回 false', () => {
      expect(isValidUUID('not-a-uuid')).toBe(false);
      expect(isValidUUID('550e8400-e29b-31d4-a716-446655440000')).toBe(false); // v3 not v4
      expect(isValidUUID('')).toBe(false);
    });
  });

  describe('validateProjectId', () => {
    test('有效的项目 ID 应通过验证', () => {
      const result = validateProjectId('550e8400-e29b-41d4-a716-446655440000');
      expect(result.valid).toBe(true);
    });

    test('空字符串应验证失败', () => {
      const result = validateProjectId('');
      expect(result.valid).toBe(false);
    });

    test('无效格式应验证失败', () => {
      const result = validateProjectId('invalid-id');
      expect(result.valid).toBe(false);
    });
  });

  describe('validateCompositionId', () => {
    test('有效的组合 ID 应通过验证', () => {
      expect(validateCompositionId('main').valid).toBe(true);
      expect(validateCompositionId('intro-scene').valid).toBe(true);
      expect(validateCompositionId('scene_01').valid).toBe(true);
    });

    test('空字符串应验证失败', () => {
      const result = validateCompositionId('');
      expect(result.valid).toBe(false);
    });

    test('包含特殊字符应验证失败', () => {
      const result = validateCompositionId('scene@01');
      expect(result.valid).toBe(false);
    });
  });

  describe('validateAssetId', () => {
    test('有效的素材 ID 应通过验证', () => {
      const result = validateAssetId('550e8400-e29b-41d4-a716-446655440000');
      expect(result.valid).toBe(true);
    });

    test('空字符串应验证失败', () => {
      const result = validateAssetId('');
      expect(result.valid).toBe(false);
    });

    test('无效格式应验证失败', () => {
      const result = validateAssetId('invalid-id');
      expect(result.valid).toBe(false);
    });
  });
});
