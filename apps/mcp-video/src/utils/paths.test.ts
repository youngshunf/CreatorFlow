/**
 * 路径工具单元测试
 */

import { describe, expect, test } from 'bun:test';
import {
  getVideoProjectsDir,
  getProjectPath,
  getProjectConfigPath,
  getAssetsPath,
  getAssetTypePath,
  getCompositionsPath,
  getCompositionFilePath,
  getOutputPath,
  getOutputFilePath,
  extractProjectName,
  getAssetRelativePath,
  getCompositionRelativePath,
  getOutputRelativePath,
  VIDEO_PROJECTS_DIR_NAME,
  ASSETS_DIR_NAME,
  COMPOSITIONS_DIR_NAME,
  OUTPUT_DIR_NAME,
} from './paths';

describe('路径工具', () => {
  const workspacePath = '/Users/user/workspace';
  const projectName = '我的视频项目';

  describe('getVideoProjectsDir', () => {
    test('应返回视频项目根目录路径', () => {
      const result = getVideoProjectsDir(workspacePath);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}`);
    });
  });

  describe('getProjectPath', () => {
    test('应返回项目目录路径', () => {
      const result = getProjectPath(workspacePath, projectName);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}`);
    });
  });

  describe('getProjectConfigPath', () => {
    test('应返回项目配置文件路径', () => {
      const result = getProjectConfigPath(workspacePath, projectName);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}/project.json`);
    });
  });

  describe('getAssetsPath', () => {
    test('应返回素材目录路径', () => {
      const result = getAssetsPath(workspacePath, projectName);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}/${ASSETS_DIR_NAME}`);
    });
  });

  describe('getAssetTypePath', () => {
    test('应返回图片素材子目录路径', () => {
      const result = getAssetTypePath(workspacePath, projectName, 'image');
      expect(result).toContain('images');
    });

    test('应返回视频素材子目录路径', () => {
      const result = getAssetTypePath(workspacePath, projectName, 'video');
      expect(result).toContain('videos');
    });

    test('应返回音频素材子目录路径', () => {
      const result = getAssetTypePath(workspacePath, projectName, 'audio');
      expect(result).toContain('audio');
    });

    test('应返回字体素材子目录路径', () => {
      const result = getAssetTypePath(workspacePath, projectName, 'font');
      expect(result).toContain('fonts');
    });
  });

  describe('getCompositionsPath', () => {
    test('应返回组合目录路径', () => {
      const result = getCompositionsPath(workspacePath, projectName);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}/${COMPOSITIONS_DIR_NAME}`);
    });
  });

  describe('getCompositionFilePath', () => {
    test('应返回组合文件路径', () => {
      const result = getCompositionFilePath(workspacePath, projectName, 'main');
      expect(result).toContain('main.tsx');
    });
  });

  describe('getOutputPath', () => {
    test('应返回输出目录路径', () => {
      const result = getOutputPath(workspacePath, projectName);
      expect(result).toBe(`${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}/${OUTPUT_DIR_NAME}`);
    });
  });

  describe('getOutputFilePath', () => {
    test('应返回输出文件路径', () => {
      const result = getOutputFilePath(workspacePath, projectName, 'render-001.mp4');
      expect(result).toContain('render-001.mp4');
    });
  });

  describe('extractProjectName', () => {
    test('应从项目路径中提取项目名称', () => {
      const projectPath = `${workspacePath}/${VIDEO_PROJECTS_DIR_NAME}/${projectName}`;
      const result = extractProjectName(projectPath);
      expect(result).toBe(projectName);
    });

    test('无效路径应返回 null', () => {
      const result = extractProjectName('/invalid/path');
      expect(result).toBeNull();
    });
  });

  describe('getAssetRelativePath', () => {
    test('应返回图片素材的相对路径', () => {
      const result = getAssetRelativePath('image', 'background.png');
      expect(result).toBe(`${ASSETS_DIR_NAME}/images/background.png`);
    });

    test('应返回视频素材的相对路径', () => {
      const result = getAssetRelativePath('video', 'intro.mp4');
      expect(result).toBe(`${ASSETS_DIR_NAME}/videos/intro.mp4`);
    });
  });

  describe('getCompositionRelativePath', () => {
    test('应返回组合的相对路径', () => {
      const result = getCompositionRelativePath('main');
      expect(result).toBe(`${COMPOSITIONS_DIR_NAME}/main.tsx`);
    });
  });

  describe('getOutputRelativePath', () => {
    test('应返回输出文件的相对路径', () => {
      const result = getOutputRelativePath('render-001.mp4');
      expect(result).toBe(`${OUTPUT_DIR_NAME}/render-001.mp4`);
    });
  });
});
