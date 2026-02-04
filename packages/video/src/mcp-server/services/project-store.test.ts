/**
 * ProjectStore 服务单元测试
 */

import { describe, expect, test, beforeEach, afterEach } from 'bun:test';
import { mkdirSync, rmSync, existsSync, readFileSync, writeFileSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';
import { ProjectStore, serializeProject, deserializeProject, isValidProject } from './project-store';
import { getProjectConfigPath, getAssetsPath, getCompositionsPath, getOutputPath, getAssetTypePath } from '../utils/paths';
import type { VideoProject } from '../types';

describe('ProjectStore', () => {
  // 创建临时测试目录
  let testDir: string;
  let store: ProjectStore;

  beforeEach(() => {
    testDir = join(tmpdir(), `mcp-video-test-${Date.now()}-${Math.random().toString(36).slice(2)}`);
    mkdirSync(testDir, { recursive: true });
    store = new ProjectStore({ workspacePath: testDir });
  });

  afterEach(() => {
    rmSync(testDir, { recursive: true, force: true });
  });

  describe('createProject', () => {
    test('应创建新项目并返回项目对象', async () => {
      const project = await store.createProject({
        name: '测试项目',
        durationInSeconds: 10,
      });

      expect(project).toBeDefined();
      expect(project.id).toBeDefined();
      expect(project.name).toBe('测试项目');
      expect(project.config.width).toBe(1920);
      expect(project.config.height).toBe(1080);
      expect(project.config.fps).toBe(30);
      expect(project.config.durationInFrames).toBe(300); // 10 * 30
      expect(project.compositions).toEqual([]);
      expect(project.assets).toEqual([]);
      expect(project.renders).toEqual([]);
    });

    test('应使用自定义配置创建项目', async () => {
      const project = await store.createProject({
        name: '自定义项目',
        width: 1280,
        height: 720,
        fps: 60,
        durationInSeconds: 5,
        description: '这是一个测试项目',
      });

      expect(project.config.width).toBe(1280);
      expect(project.config.height).toBe(720);
      expect(project.config.fps).toBe(60);
      expect(project.config.durationInFrames).toBe(300); // 5 * 60
      expect(project.description).toBe('这是一个测试项目');
    });

    test('应创建项目目录结构', async () => {
      const project = await store.createProject({
        name: '目录测试',
        durationInSeconds: 10,
      });

      // 检查项目配置文件
      const configPath = getProjectConfigPath(testDir, '目录测试');
      expect(existsSync(configPath)).toBe(true);

      // 检查素材目录
      const assetsPath = getAssetsPath(testDir, '目录测试');
      expect(existsSync(assetsPath)).toBe(true);
      expect(existsSync(join(assetsPath, 'images'))).toBe(true);
      expect(existsSync(join(assetsPath, 'videos'))).toBe(true);
      expect(existsSync(join(assetsPath, 'audio'))).toBe(true);
      expect(existsSync(join(assetsPath, 'fonts'))).toBe(true);

      // 检查组合目录
      const compositionsPath = getCompositionsPath(testDir, '目录测试');
      expect(existsSync(compositionsPath)).toBe(true);

      // 检查输出目录
      const outputPath = getOutputPath(testDir, '目录测试');
      expect(existsSync(outputPath)).toBe(true);
    });

    test('应持久化项目数据到 JSON 文件', async () => {
      const project = await store.createProject({
        name: '持久化测试',
        durationInSeconds: 10,
      });

      const configPath = getProjectConfigPath(testDir, '持久化测试');
      const content = readFileSync(configPath, 'utf-8');
      const savedProject = JSON.parse(content);

      expect(savedProject.id).toBe(project.id);
      expect(savedProject.name).toBe(project.name);
      expect(savedProject.config).toEqual(project.config);
    });

    test('项目名称为空时应抛出错误', async () => {
      await expect(store.createProject({
        name: '',
        durationInSeconds: 10,
      })).rejects.toThrow();
    });

    test('项目已存在时应抛出错误', async () => {
      await store.createProject({
        name: '重复项目',
        durationInSeconds: 10,
      });

      await expect(store.createProject({
        name: '重复项目',
        durationInSeconds: 10,
      })).rejects.toThrow();
    });
  });

  describe('getProject', () => {
    test('应通过 ID 获取项目', async () => {
      const created = await store.createProject({
        name: '获取测试',
        durationInSeconds: 10,
      });

      const project = await store.getProject(created.id);

      expect(project).toBeDefined();
      expect(project?.id).toBe(created.id);
      expect(project?.name).toBe('获取测试');
    });

    test('项目不存在时应返回 null', async () => {
      const project = await store.getProject('nonexistent-id');
      expect(project).toBeNull();
    });
  });

  describe('getProjectByName', () => {
    test('应通过名称获取项目', async () => {
      const created = await store.createProject({
        name: '名称获取测试',
        durationInSeconds: 10,
      });

      const project = await store.getProjectByName('名称获取测试');

      expect(project).toBeDefined();
      expect(project?.id).toBe(created.id);
    });

    test('项目不存在时应返回 null', async () => {
      const project = await store.getProjectByName('不存在的项目');
      expect(project).toBeNull();
    });
  });

  describe('listProjects', () => {
    test('应列出所有项目', async () => {
      await store.createProject({ name: '项目1', durationInSeconds: 10 });
      await store.createProject({ name: '项目2', durationInSeconds: 10 });
      await store.createProject({ name: '项目3', durationInSeconds: 10 });

      const projects = await store.listProjects();

      expect(projects.length).toBe(3);
    });

    test('应按更新时间降序排列', async () => {
      const p1 = await store.createProject({ name: '项目A', durationInSeconds: 10 });
      await new Promise(resolve => setTimeout(resolve, 10)); // 等待一小段时间
      const p2 = await store.createProject({ name: '项目B', durationInSeconds: 10 });
      await new Promise(resolve => setTimeout(resolve, 10));
      const p3 = await store.createProject({ name: '项目C', durationInSeconds: 10 });

      const projects = await store.listProjects();

      // 最新创建的应该在前面
      expect(projects[0]!.name).toBe('项目C');
      expect(projects[1]!.name).toBe('项目B');
      expect(projects[2]!.name).toBe('项目A');
    });

    test('没有项目时应返回空数组', async () => {
      const projects = await store.listProjects();
      expect(projects).toEqual([]);
    });

    test('应返回正确的项目摘要信息', async () => {
      await store.createProject({
        name: '摘要测试',
        durationInSeconds: 10,
        description: '测试描述',
      });

      const projects = await store.listProjects();

      expect(projects[0]!.name).toBe('摘要测试');
      expect(projects[0]!.description).toBe('测试描述');
      expect(projects[0]!.compositionCount).toBe(0);
      expect(projects[0]!.assetCount).toBe(0);
    });
  });

  describe('updateProject', () => {
    test('应更新项目名称', async () => {
      const created = await store.createProject({
        name: '原始名称',
        durationInSeconds: 10,
      });

      const updated = await store.updateProject(created.id, {
        name: '新名称',
      });

      expect(updated.name).toBe('新名称');

      // 验证文件系统中的目录也被重命名
      const newConfigPath = getProjectConfigPath(testDir, '新名称');
      expect(existsSync(newConfigPath)).toBe(true);
    });

    test('应更新项目描述', async () => {
      const created = await store.createProject({
        name: '描述测试',
        durationInSeconds: 10,
      });

      const updated = await store.updateProject(created.id, {
        description: '新描述',
      });

      expect(updated.description).toBe('新描述');
    });

    test('应更新项目配置', async () => {
      const created = await store.createProject({
        name: '配置测试',
        durationInSeconds: 10,
      });

      const updated = await store.updateProject(created.id, {
        config: {
          width: 1280,
          height: 720,
        },
      });

      expect(updated.config.width).toBe(1280);
      expect(updated.config.height).toBe(720);
      expect(updated.config.fps).toBe(30); // 未更新的字段保持不变
    });

    test('应自动更新 updatedAt 时间戳', async () => {
      const created = await store.createProject({
        name: '时间戳测试',
        durationInSeconds: 10,
      });

      await new Promise(resolve => setTimeout(resolve, 10));

      const updated = await store.updateProject(created.id, {
        description: '更新描述',
      });

      expect(new Date(updated.updatedAt).getTime()).toBeGreaterThan(
        new Date(created.updatedAt).getTime()
      );
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.updateProject('nonexistent-id', {
        name: '新名称',
      })).rejects.toThrow();
    });

    test('新名称已存在时应抛出错误', async () => {
      await store.createProject({ name: '项目A', durationInSeconds: 10 });
      const projectB = await store.createProject({ name: '项目B', durationInSeconds: 10 });

      await expect(store.updateProject(projectB.id, {
        name: '项目A',
      })).rejects.toThrow();
    });
  });

  describe('deleteProject', () => {
    test('应删除项目', async () => {
      const created = await store.createProject({
        name: '删除测试',
        durationInSeconds: 10,
      });

      const result = await store.deleteProject(created.id);

      expect(result).toBe(true);

      // 验证项目已被删除
      const project = await store.getProject(created.id);
      expect(project).toBeNull();

      // 验证目录已被删除
      const configPath = getProjectConfigPath(testDir, '删除测试');
      expect(existsSync(configPath)).toBe(false);
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.deleteProject('nonexistent-id')).rejects.toThrow();
    });
  });

  // ==========================================================================
  // 素材管理测试
  // ==========================================================================

  describe('addAsset', () => {
    let testAssetDir: string;

    beforeEach(() => {
      // 创建测试素材目录
      testAssetDir = join(testDir, 'test-assets');
      mkdirSync(testAssetDir, { recursive: true });
    });

    test('应添加图片素材到项目', async () => {
      // 创建测试图片文件
      const testImagePath = join(testAssetDir, 'test-image.png');
      writeFileSync(testImagePath, 'fake png content');

      const project = await store.createProject({
        name: '素材测试',
        durationInSeconds: 10,
      });

      const asset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      expect(asset).toBeDefined();
      expect(asset.id).toBeDefined();
      expect(asset.type).toBe('image');
      expect(asset.name).toBe('test-image.png');
      expect(asset.path).toContain('素材/images/test-image.png');
    });

    test('应复制文件到项目素材目录', async () => {
      const testImagePath = join(testAssetDir, 'copy-test.png');
      writeFileSync(testImagePath, 'test content for copy');

      const project = await store.createProject({
        name: '复制测试',
        durationInSeconds: 10,
      });

      await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 验证文件已复制到项目目录
      const targetPath = join(getAssetTypePath(testDir, '复制测试', 'image'), 'copy-test.png');
      expect(existsSync(targetPath)).toBe(true);
      expect(readFileSync(targetPath, 'utf-8')).toBe('test content for copy');
    });

    test('应在项目配置中注册素材', async () => {
      const testImagePath = join(testAssetDir, 'register-test.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '注册测试',
        durationInSeconds: 10,
      });

      await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 重新获取项目验证素材已注册
      const updatedProject = await store.getProject(project.id);
      expect(updatedProject?.assets.length).toBe(1);
      expect(updatedProject?.assets[0]?.type).toBe('image');
    });

    test('应使用自定义名称', async () => {
      const testImagePath = join(testAssetDir, 'original-name.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '自定义名称测试',
        durationInSeconds: 10,
      });

      const asset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
        name: '我的自定义图片',
      });

      expect(asset.name).toBe('我的自定义图片');
    });

    test('应处理重复文件名', async () => {
      const testImagePath = join(testAssetDir, 'duplicate.png');
      writeFileSync(testImagePath, 'first content');

      const project = await store.createProject({
        name: '重复文件名测试',
        durationInSeconds: 10,
      });

      // 添加第一个文件
      const asset1 = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 修改源文件内容
      writeFileSync(testImagePath, 'second content');

      // 添加同名文件
      const asset2 = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 两个素材应该有不同的路径
      expect(asset1.path).not.toBe(asset2.path);
      expect(asset2.path).toContain('duplicate_1.png');
    });

    test('应支持不同类型的素材', async () => {
      const project = await store.createProject({
        name: '多类型测试',
        durationInSeconds: 10,
      });

      // 创建不同类型的测试文件
      const videoPath = join(testAssetDir, 'test.mp4');
      const audioPath = join(testAssetDir, 'test.mp3');
      const fontPath = join(testAssetDir, 'test.ttf');

      writeFileSync(videoPath, 'fake video');
      writeFileSync(audioPath, 'fake audio');
      writeFileSync(fontPath, 'fake font');

      const videoAsset = await store.addAsset(project.id, {
        sourcePath: videoPath,
        assetType: 'video',
      });
      const audioAsset = await store.addAsset(project.id, {
        sourcePath: audioPath,
        assetType: 'audio',
      });
      const fontAsset = await store.addAsset(project.id, {
        sourcePath: fontPath,
        assetType: 'font',
      });

      expect(videoAsset.type).toBe('video');
      expect(videoAsset.path).toContain('videos');
      expect(audioAsset.type).toBe('audio');
      expect(audioAsset.path).toContain('audio');
      expect(fontAsset.type).toBe('font');
      expect(fontAsset.path).toContain('fonts');
    });

    test('项目不存在时应抛出错误', async () => {
      const testImagePath = join(testAssetDir, 'test.png');
      writeFileSync(testImagePath, 'test content');

      await expect(store.addAsset('nonexistent-id', {
        sourcePath: testImagePath,
        assetType: 'image',
      })).rejects.toThrow();
    });

    test('文件不存在时应抛出错误', async () => {
      const project = await store.createProject({
        name: '文件不存在测试',
        durationInSeconds: 10,
      });

      await expect(store.addAsset(project.id, {
        sourcePath: '/nonexistent/path/image.png',
        assetType: 'image',
      })).rejects.toThrow();
    });

    test('不支持的格式应抛出错误', async () => {
      const testFilePath = join(testAssetDir, 'test.xyz');
      writeFileSync(testFilePath, 'test content');

      const project = await store.createProject({
        name: '格式测试',
        durationInSeconds: 10,
      });

      await expect(store.addAsset(project.id, {
        sourcePath: testFilePath,
        assetType: 'image',
      })).rejects.toThrow();
    });
  });

  describe('removeAsset', () => {
    let testAssetDir: string;

    beforeEach(() => {
      testAssetDir = join(testDir, 'test-assets');
      mkdirSync(testAssetDir, { recursive: true });
    });

    test('应从项目中移除素材注册', async () => {
      const testImagePath = join(testAssetDir, 'remove-test.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '移除测试',
        durationInSeconds: 10,
      });

      const asset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      const result = await store.removeAsset(project.id, asset.id);

      expect(result).toBe(true);

      // 验证素材已从项目中移除
      const updatedProject = await store.getProject(project.id);
      expect(updatedProject?.assets.length).toBe(0);
    });

    test('默认不删除实际文件', async () => {
      const testImagePath = join(testAssetDir, 'keep-file.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '保留文件测试',
        durationInSeconds: 10,
      });

      const asset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 获取复制后的文件路径
      const copiedFilePath = join(store.getProjectPath('保留文件测试'), asset.path);

      await store.removeAsset(project.id, asset.id, false);

      // 文件应该仍然存在
      expect(existsSync(copiedFilePath)).toBe(true);
    });

    test('deleteFile=true 时应删除实际文件', async () => {
      const testImagePath = join(testAssetDir, 'delete-file.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '删除文件测试',
        durationInSeconds: 10,
      });

      const asset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
      });

      // 获取复制后的文件路径
      const copiedFilePath = join(store.getProjectPath('删除文件测试'), asset.path);
      expect(existsSync(copiedFilePath)).toBe(true);

      await store.removeAsset(project.id, asset.id, true);

      // 文件应该被删除
      expect(existsSync(copiedFilePath)).toBe(false);
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.removeAsset('nonexistent-id', 'asset-id')).rejects.toThrow();
    });

    test('素材不存在时应抛出错误', async () => {
      const project = await store.createProject({
        name: '素材不存在测试',
        durationInSeconds: 10,
      });

      await expect(store.removeAsset(project.id, 'nonexistent-asset-id')).rejects.toThrow();
    });
  });

  describe('listAssets', () => {
    let testAssetDir: string;

    beforeEach(() => {
      testAssetDir = join(testDir, 'test-assets');
      mkdirSync(testAssetDir, { recursive: true });
    });

    test('应按类型分组返回素材', async () => {
      const project = await store.createProject({
        name: '分组测试',
        durationInSeconds: 10,
      });

      // 创建不同类型的测试文件
      const imagePath = join(testAssetDir, 'test.png');
      const videoPath = join(testAssetDir, 'test.mp4');
      const audioPath = join(testAssetDir, 'test.mp3');

      writeFileSync(imagePath, 'fake image');
      writeFileSync(videoPath, 'fake video');
      writeFileSync(audioPath, 'fake audio');

      await store.addAsset(project.id, { sourcePath: imagePath, assetType: 'image' });
      await store.addAsset(project.id, { sourcePath: videoPath, assetType: 'video' });
      await store.addAsset(project.id, { sourcePath: audioPath, assetType: 'audio' });

      const assetsByType = await store.listAssets(project.id);

      expect(assetsByType.image.length).toBe(1);
      expect(assetsByType.video.length).toBe(1);
      expect(assetsByType.audio.length).toBe(1);
      expect(assetsByType.font.length).toBe(0);
    });

    test('没有素材时应返回空分组', async () => {
      const project = await store.createProject({
        name: '空素材测试',
        durationInSeconds: 10,
      });

      const assetsByType = await store.listAssets(project.id);

      expect(assetsByType.image).toEqual([]);
      expect(assetsByType.video).toEqual([]);
      expect(assetsByType.audio).toEqual([]);
      expect(assetsByType.font).toEqual([]);
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.listAssets('nonexistent-id')).rejects.toThrow();
    });
  });

  describe('getAsset', () => {
    let testAssetDir: string;

    beforeEach(() => {
      testAssetDir = join(testDir, 'test-assets');
      mkdirSync(testAssetDir, { recursive: true });
    });

    test('应返回素材详情', async () => {
      const testImagePath = join(testAssetDir, 'get-test.png');
      writeFileSync(testImagePath, 'test content');

      const project = await store.createProject({
        name: '获取素材测试',
        durationInSeconds: 10,
      });

      const addedAsset = await store.addAsset(project.id, {
        sourcePath: testImagePath,
        assetType: 'image',
        name: '测试图片',
      });

      const asset = await store.getAsset(project.id, addedAsset.id);

      expect(asset).toBeDefined();
      expect(asset?.id).toBe(addedAsset.id);
      expect(asset?.name).toBe('测试图片');
      expect(asset?.type).toBe('image');
    });

    test('素材不存在时应返回 null', async () => {
      const project = await store.createProject({
        name: '素材不存在测试',
        durationInSeconds: 10,
      });

      const asset = await store.getAsset(project.id, 'nonexistent-asset-id');
      expect(asset).toBeNull();
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.getAsset('nonexistent-id', 'asset-id')).rejects.toThrow();
    });
  });

  // ==========================================================================
  // 组合管理测试
  // ==========================================================================

  describe('addComposition', () => {
    test('应添加组合到项目', async () => {
      const project = await store.createProject({
        name: '添加组合测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '主场景',
        code: 'export default function Main() { return <div>Hello</div>; }',
        props: { title: '测试标题' },
      });

      expect(composition).toBeDefined();
      expect(composition.id).toBeDefined();
      expect(composition.name).toBe('主场景');
      expect(composition.code).toBe('export default function Main() { return <div>Hello</div>; }');
      expect(composition.props).toEqual({ title: '测试标题' });
    });

    test('应在项目配置中注册组合', async () => {
      const project = await store.createProject({
        name: '注册组合测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'export default function Test() { return null; }',
      });

      // 重新获取项目验证组合已注册
      const updatedProject = await store.getProject(project.id);
      expect(updatedProject?.compositions).toHaveLength(1);
      expect(updatedProject?.compositions[0]?.id).toBe(composition.id);
      expect(updatedProject?.compositions[0]?.name).toBe('测试组合');
    });

    test('应支持不带 props 的组合', async () => {
      const project = await store.createProject({
        name: '无props组合测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '简单组合',
        code: 'export default function Simple() { return <div />; }',
      });

      expect(composition.props).toEqual({});
    });

    test('应更新项目的 updatedAt 时间戳', async () => {
      const project = await store.createProject({
        name: '时间戳测试',
        durationInSeconds: 10,
      });

      const originalUpdatedAt = project.updatedAt;

      // 等待一小段时间确保时间戳不同
      await new Promise(resolve => setTimeout(resolve, 10));

      await store.addComposition(project.id, {
        name: '新组合',
        code: 'export default function New() { return null; }',
      });

      const updatedProject = await store.getProject(project.id);
      expect(new Date(updatedProject!.updatedAt).getTime()).toBeGreaterThan(
        new Date(originalUpdatedAt).getTime()
      );
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(
        store.addComposition('nonexistent-id', {
          name: '测试',
          code: 'code',
        })
      ).rejects.toThrow();
    });

    test('组合名称为空时应抛出错误', async () => {
      const project = await store.createProject({
        name: '空名称测试',
        durationInSeconds: 10,
      });

      await expect(
        store.addComposition(project.id, {
          name: '',
          code: 'export default function Test() { return null; }',
        })
      ).rejects.toThrow();
    });

    test('组合代码为空时应抛出错误', async () => {
      const project = await store.createProject({
        name: '空代码测试',
        durationInSeconds: 10,
      });

      await expect(
        store.addComposition(project.id, {
          name: '测试组合',
          code: '',
        })
      ).rejects.toThrow();
    });
  });

  describe('updateComposition', () => {
    test('应更新组合名称', async () => {
      const project = await store.createProject({
        name: '更新名称测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '原始名称',
        code: 'export default function Test() { return null; }',
      });

      const updated = await store.updateComposition(project.id, composition.id, {
        name: '新名称',
      });

      expect(updated.name).toBe('新名称');
      expect(updated.code).toBe(composition.code); // 代码不变
    });

    test('应更新组合代码', async () => {
      const project = await store.createProject({
        name: '更新代码测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'export default function Old() { return null; }',
      });

      const updated = await store.updateComposition(project.id, composition.id, {
        code: 'export default function New() { return <div>New</div>; }',
      });

      expect(updated.code).toBe('export default function New() { return <div>New</div>; }');
      expect(updated.name).toBe(composition.name); // 名称不变
    });

    test('应更新组合属性', async () => {
      const project = await store.createProject({
        name: '更新属性测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'export default function Test() { return null; }',
        props: { oldProp: 'old' },
      });

      const updated = await store.updateComposition(project.id, composition.id, {
        props: { newProp: 'new', count: 42 },
      });

      expect(updated.props).toEqual({ newProp: 'new', count: 42 });
    });

    test('应同时更新多个字段', async () => {
      const project = await store.createProject({
        name: '多字段更新测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '原始组合',
        code: 'old code',
        props: { old: true },
      });

      const updated = await store.updateComposition(project.id, composition.id, {
        name: '新组合',
        code: 'new code',
        props: { new: true },
      });

      expect(updated.name).toBe('新组合');
      expect(updated.code).toBe('new code');
      expect(updated.props).toEqual({ new: true });
    });

    test('应持久化更新到文件', async () => {
      const project = await store.createProject({
        name: '持久化更新测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'original code',
      });

      await store.updateComposition(project.id, composition.id, {
        name: '更新后的组合',
      });

      // 创建新的 store 实例验证持久化
      const newStore = new ProjectStore({ workspacePath: testDir });
      const loadedProject = await newStore.getProject(project.id);

      expect(loadedProject?.compositions[0]?.name).toBe('更新后的组合');
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(
        store.updateComposition('nonexistent-id', 'comp-id', { name: 'new' })
      ).rejects.toThrow();
    });

    test('组合不存在时应抛出错误', async () => {
      const project = await store.createProject({
        name: '组合不存在测试',
        durationInSeconds: 10,
      });

      await expect(
        store.updateComposition(project.id, 'nonexistent-comp-id', { name: 'new' })
      ).rejects.toThrow();
    });

    test('更新名称为空时应抛出错误', async () => {
      const project = await store.createProject({
        name: '空名称更新测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'code',
      });

      await expect(
        store.updateComposition(project.id, composition.id, { name: '' })
      ).rejects.toThrow();
    });

    test('更新代码为空时应抛出错误', async () => {
      const project = await store.createProject({
        name: '空代码更新测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'code',
      });

      await expect(
        store.updateComposition(project.id, composition.id, { code: '' })
      ).rejects.toThrow();
    });
  });

  describe('removeComposition', () => {
    test('应从项目中移除组合', async () => {
      const project = await store.createProject({
        name: '移除组合测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '待移除组合',
        code: 'export default function Test() { return null; }',
      });

      const result = await store.removeComposition(project.id, composition.id);
      expect(result).toBe(true);

      // 验证组合已被移除
      const updatedProject = await store.getProject(project.id);
      expect(updatedProject?.compositions).toHaveLength(0);
    });

    test('应只移除指定的组合', async () => {
      const project = await store.createProject({
        name: '选择性移除测试',
        durationInSeconds: 10,
      });

      const comp1 = await store.addComposition(project.id, {
        name: '组合1',
        code: 'code1',
      });

      const comp2 = await store.addComposition(project.id, {
        name: '组合2',
        code: 'code2',
      });

      await store.removeComposition(project.id, comp1.id);

      const updatedProject = await store.getProject(project.id);
      expect(updatedProject?.compositions).toHaveLength(1);
      expect(updatedProject?.compositions[0]?.id).toBe(comp2.id);
    });

    test('应更新项目的 updatedAt 时间戳', async () => {
      const project = await store.createProject({
        name: '移除时间戳测试',
        durationInSeconds: 10,
      });

      const composition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'code',
      });

      const projectAfterAdd = await store.getProject(project.id);
      const addedAt = projectAfterAdd!.updatedAt;

      // 等待一小段时间确保时间戳不同
      await new Promise(resolve => setTimeout(resolve, 10));

      await store.removeComposition(project.id, composition.id);

      const updatedProject = await store.getProject(project.id);
      expect(new Date(updatedProject!.updatedAt).getTime()).toBeGreaterThan(
        new Date(addedAt).getTime()
      );
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(
        store.removeComposition('nonexistent-id', 'comp-id')
      ).rejects.toThrow();
    });

    test('组合不存在时应抛出错误', async () => {
      const project = await store.createProject({
        name: '组合不存在移除测试',
        durationInSeconds: 10,
      });

      await expect(
        store.removeComposition(project.id, 'nonexistent-comp-id')
      ).rejects.toThrow();
    });
  });

  describe('listCompositions', () => {
    test('应返回项目中的所有组合', async () => {
      const project = await store.createProject({
        name: '列出组合测试',
        durationInSeconds: 10,
      });

      await store.addComposition(project.id, {
        name: '组合1',
        code: 'code1',
      });

      await store.addComposition(project.id, {
        name: '组合2',
        code: 'code2',
      });

      const compositions = await store.listCompositions(project.id);

      expect(compositions).toHaveLength(2);
      expect(compositions.map(c => c.name)).toContain('组合1');
      expect(compositions.map(c => c.name)).toContain('组合2');
    });

    test('没有组合时应返回空数组', async () => {
      const project = await store.createProject({
        name: '空组合列表测试',
        durationInSeconds: 10,
      });

      const compositions = await store.listCompositions(project.id);
      expect(compositions).toEqual([]);
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(store.listCompositions('nonexistent-id')).rejects.toThrow();
    });
  });

  describe('getComposition', () => {
    test('应返回组合详情', async () => {
      const project = await store.createProject({
        name: '获取组合测试',
        durationInSeconds: 10,
      });

      const addedComposition = await store.addComposition(project.id, {
        name: '测试组合',
        code: 'test code',
        props: { key: 'value' },
      });

      const composition = await store.getComposition(project.id, addedComposition.id);

      expect(composition).toBeDefined();
      expect(composition?.id).toBe(addedComposition.id);
      expect(composition?.name).toBe('测试组合');
      expect(composition?.code).toBe('test code');
      expect(composition?.props).toEqual({ key: 'value' });
    });

    test('组合不存在时应返回 null', async () => {
      const project = await store.createProject({
        name: '组合不存在获取测试',
        durationInSeconds: 10,
      });

      const composition = await store.getComposition(project.id, 'nonexistent-comp-id');
      expect(composition).toBeNull();
    });

    test('项目不存在时应抛出错误', async () => {
      await expect(
        store.getComposition('nonexistent-id', 'comp-id')
      ).rejects.toThrow();
    });
  });

  describe('路径工具方法', () => {
    test('getProjectPath 应返回正确路径', () => {
      const path = store.getProjectPath('测试项目');
      expect(path).toContain('测试项目');
    });

    test('getAssetsPath 应返回正确路径', () => {
      const path = store.getAssetsPath('测试项目');
      expect(path).toContain('素材');
    });

    test('getOutputPath 应返回正确路径', () => {
      const path = store.getOutputPath('测试项目');
      expect(path).toContain('输出');
    });

    test('getCompositionsPath 应返回正确路径', () => {
      const path = store.getCompositionsPath('测试项目');
      expect(path).toContain('组合');
    });
  });

  describe('静态工厂方法', () => {
    test('create 应返回 ProjectStore 实例', () => {
      const newStore = ProjectStore.create(testDir);
      expect(newStore).toBeInstanceOf(ProjectStore);
    });
  });
});

describe('序列化/反序列化', () => {
  const sampleProject: VideoProject = {
    id: '550e8400-e29b-41d4-a716-446655440000',
    name: '测试项目',
    description: '测试描述',
    createdAt: '2024-01-01T00:00:00.000Z',
    updatedAt: '2024-01-01T00:00:00.000Z',
    config: {
      width: 1920,
      height: 1080,
      fps: 30,
      durationInFrames: 300,
    },
    compositions: [],
    assets: [],
    renders: [],
  };

  describe('serializeProject', () => {
    test('应将项目序列化为 JSON 字符串', () => {
      const json = serializeProject(sampleProject);
      expect(typeof json).toBe('string');
      expect(JSON.parse(json)).toEqual(sampleProject);
    });
  });

  describe('deserializeProject', () => {
    test('应将 JSON 字符串反序列化为项目对象', () => {
      const json = JSON.stringify(sampleProject);
      const project = deserializeProject(json);
      expect(project).toEqual(sampleProject);
    });

    test('无效 JSON 应返回 null', () => {
      const project = deserializeProject('invalid json');
      expect(project).toBeNull();
    });

    test('缺少必要字段应返回 null', () => {
      const project = deserializeProject('{"name": "test"}');
      expect(project).toBeNull();
    });
  });

  describe('isValidProject', () => {
    test('有效项目应返回 true', () => {
      expect(isValidProject(sampleProject)).toBe(true);
    });

    test('无效项目应返回 false', () => {
      expect(isValidProject({})).toBe(false);
      expect(isValidProject({ name: 'test' })).toBe(false);
      expect(isValidProject(null)).toBe(false);
      expect(isValidProject(undefined)).toBe(false);
    });
  });

  describe('序列化往返测试', () => {
    test('序列化后反序列化应得到等价对象', () => {
      const json = serializeProject(sampleProject);
      const deserialized = deserializeProject(json);
      expect(deserialized).toEqual(sampleProject);
    });
  });
});
