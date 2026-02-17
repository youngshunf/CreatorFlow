================================================================================
视频创作功能设计方案 v3 — 场景编排器（Scene Composer）
================================================================================

状态：设计评审中
日期：2026-02-14
前置文档：remotion-integration-plan-v2.md


一、问题分析
================================================================================

1.1 当前架构的三个断裂点
--------------------------------------------------------------------------------

【预览断裂】
VideoPreview.tsx 的 COMPOSITION_MAP 只映射了 4 个内置组合（TitleAnimation、
Slideshow、DataVisualization、ProductShowcase）。Agent 通过 MCP 创建的组合存储
为代码字符串（compositionCode），无法被 @remotion/player 渲染。

相关文件: apps/electron/src/renderer/components/video/VideoPreview.tsx

【渲染断裂】
RenderEngine.findEntryPoint() 在用户项目目录查找 Root.tsx，但项目目录只有
project.json 和空文件夹。@remotion/bundler.bundle() 找不到入口文件，渲染必然失败。

相关文件: packages/video/src/mcp-server/services/render-engine.ts

【用户体验断裂】
当前架构要求用户理解代码、安装依赖、使用终端。目标用户是非程序员的自媒体创作者，
大部分使用 Windows 系统。


1.2 Agent 对话记录中的实际失败
--------------------------------------------------------------------------------

来自 session 260213-brave-island：

- video_preview_start → PREVIEW_FAILED: 服务器进程意外退出 (code: 1)
  原因: PreviewServer 在项目目录执行 npx remotion studio，但目录无 Remotion 环境

- video_render → RENDER_FAILED: 视频项目打包失败
  原因: RenderEngine.bundle() 在项目目录找不到 Root.tsx 入口文件


1.3 模板 compositionCode 是死代码
--------------------------------------------------------------------------------

模板文件（social-media.ts、marketing.ts、tutorial.ts）中存储了几百行 TSX 代码字符串，
但这些代码：

- 不会被 @remotion/player 执行（Player 需要 React 组件引用）
- 不会被写入文件系统（MCP 只存到 project.json）
- 不会被 @remotion/renderer 编译（项目目录无 Remotion 环境）


1.4 Remotion 最佳实践违规
--------------------------------------------------------------------------------

- 模板代码中使用 CSS transition 属性 → 违反：所有动画必须用 useCurrentFrame()+interpolate()
- 组合缺少 Zod Schema → 违反：Remotion 推荐用 Zod schema 让组合可参数化
- 未使用 calculateMetadata → 违反：无法动态计算视频时长和尺寸
- 未使用 TransitionSeries → 违反：缺少场景间过渡效果能力
- 未使用 premountFor → 违反：Sequence 应始终设置 premount


二、设计目标
================================================================================

2.1 用户体验目标
--------------------------------------------------------------------------------

::

  用户在 APP 中说："帮我做一个产品推广视频"
      ↓
  Agent 自动编排多个场景片段，配置参数
      ↓
  用户在 APP 内实时看到视频预览
      ↓
  用户说："标题改大一点，颜色换成蓝色"
      ↓
  Agent 更新参数，预览实时刷新
      ↓
  用户满意 → 点击导出 → 获得 MP4 文件

关键约束：

- 用户不接触代码、终端、依赖安装
- 预览和导出都在 APP 内完成
- Windows/Mac/Linux 全平台兼容
- Agent 有足够的创作自由度


2.2 技术目标
--------------------------------------------------------------------------------

- 预览：@remotion/player 在 Electron renderer 进程直接渲染
- 导出：Electron 主进程用 @remotion/renderer + packages/video/src/Root.tsx 统一入口
- 数据：项目存储为 JSON（场景列表 + 参数），不存储代码
- 扩展：新增内置组合 = 新增"积木块"，Agent 创作能力随之增长


三、核心架构：场景编排器（Scene Composer）
================================================================================

3.1 设计理念
--------------------------------------------------------------------------------

不是让 agent 写代码，也不是限制 agent 只能选一个模板。
而是给 agent 一个「场景编排系统」：

- Agent 可以自由组合多个内置场景片段
- 配置每个片段的参数（文字、颜色、图片、动画、布局）
- 设置片段间的过渡效果（fade、slide、wipe、flip）
- 控制每个片段的时长
- 最终由 SceneComposer 组件用 TransitionSeries 拼接渲染

这就像视频剪辑软件的"轨道"概念 — 每个轨道放一个组合片段，agent 自由编排。


3.2 数据流
--------------------------------------------------------------------------------

::

  Agent (MCP) ──选择场景+配置参数──→ project.json (scenes[] + transitions[])
                                           │
                      ┌────────────────────┤
                      ▼                    ▼
              Electron 预览            Electron 导出
           (@remotion/player)      (主进程 @remotion/renderer)
           SceneComposer 组件       使用 packages/video/Root.tsx
           + TransitionSeries       作为统一入口打包渲染
           直接在 renderer 进程     输出 mp4/webm/gif 文件


3.3 数据模型
--------------------------------------------------------------------------------

VideoProject（修改后）::

  interface VideoProject {
    id: string;
    name: string;
    description?: string;
    createdAt: string;
    updatedAt: string;
    config: VideoConfig;        // { width, height, fps }
    scenes: Scene[];            // 场景片段列表（有序）
    transitions: Transition[];  // 片段间过渡效果（长度 = scenes.length - 1）
    assets: Asset[];
    renders: RenderHistory[];
  }

Scene（新增）::

  interface Scene {
    id: string;
    name: string;
    compositionId: string;      // 引用内置组合 ID，如 "TitleAnimation"
    durationInFrames: number;
    props: Record<string, any>; // 传给组合组件的参数
  }

Transition（新增）::

  interface Transition {
    type: 'fade' | 'slide' | 'wipe' | 'flip' | 'clock-wipe' | 'none';
    durationInFrames: number;
    direction?: 'from-left' | 'from-right' | 'from-top' | 'from-bottom';
  }

对比旧的 Composition::

  // 旧（删除）
  interface Composition {
    id: string;
    name: string;
    code: string;           // ← 代码字符串，无法执行
    props: Record<string, any>;
  }

  // 新（Scene 替代）
  interface Scene {
    id: string;
    name: string;
    compositionId: string;  // ← 引用内置组合 ID
    durationInFrames: number;
    props: Record<string, any>;
  }


四、内置组合库
================================================================================

4.1 现有组合（保留）
--------------------------------------------------------------------------------

- TitleAnimation — 标题动画（fade/slide/scale/spring 四种风格）
- Slideshow — 图片幻灯片（fade/slide/zoom/crossfade 过渡）
- DataVisualization — 数据图表（bar/line/pie/donut 四种图表）
- ProductShowcase — 产品展示（centered/split/features-grid 三种布局）


4.2 新增组合（从模板代码提取）
--------------------------------------------------------------------------------

从现有模板的 compositionCode 字符串提取为真正的 React 组件：

- SocialMediaVertical — 竖版社交媒体（9:16，来自 social-media.ts）
- SocialMediaSquare — 方形社交媒体（1:1，来自 social-media.ts）
- StepByStepTutorial — 分步教程（来自 tutorial.ts）
- Explainer — 概念讲解（来自 tutorial.ts）
- Tips — 技巧清单（来自 tutorial.ts）
- ProductMarketing — 产品营销（来自 marketing.ts）
- PromoAd — 促销广告（来自 marketing.ts）

提取时需修复 Remotion 违规：
- 删除所有 CSS transition 属性
- 用 useCurrentFrame() + interpolate() 替代


4.3 每个组合必须提供
--------------------------------------------------------------------------------

- React 组件（.tsx 文件）
- Zod Schema（用于参数验证和 Remotion Studio 编辑）
- 默认 props
- calculateMetadata 函数（动态计算时长）


4.4 组合注册
--------------------------------------------------------------------------------

所有组合在 packages/video/src/Root.tsx 中注册::

  // Root.tsx
  import { Composition, Folder } from 'remotion';
  import { SceneComposer, SceneComposerSchema } from './compositions/SceneComposer';
  import { TitleAnimation, TitleAnimationSchema } from './compositions/TitleAnimation';
  // ... 其他组合

  export const RemotionRoot = () => (
    <>
      {/* 核心：场景编排器组合 */}
      <Composition
        id="SceneComposer"
        component={SceneComposer}
        schema={SceneComposerSchema}
        calculateMetadata={calculateSceneComposerMetadata}
        // ...
      />

      {/* 单独组合（也可独立使用） */}
      <Folder name="基础">
        <Composition id="TitleAnimation" component={TitleAnimation} schema={TitleAnimationSchema} ... />
        <Composition id="Slideshow" component={Slideshow} schema={SlideshowSchema} ... />
      </Folder>
      <Folder name="社交媒体">
        <Composition id="SocialMediaVertical" ... />
        <Composition id="SocialMediaSquare" ... />
      </Folder>
      {/* ... */}
    </>
  );


五、SceneComposer 核心组件
================================================================================

5.1 组件设计
--------------------------------------------------------------------------------

SceneComposer 是一个特殊的 Remotion 组合，用 TransitionSeries 编排多个场景::

  // packages/video/src/compositions/SceneComposer.tsx

  import { TransitionSeries, linearTiming } from '@remotion/transitions';
  import { fade } from '@remotion/transitions/fade';
  import { slide } from '@remotion/transitions/slide';
  import { wipe } from '@remotion/transitions/wipe';
  import { flip } from '@remotion/transitions/flip';
  import { clockWipe } from '@remotion/transitions/clock-wipe';

  const COMPOSITION_MAP: Record<string, React.FC<any>> = {
    TitleAnimation,
    Slideshow,
    DataVisualization,
    ProductShowcase,
    SocialMediaVertical,
    SocialMediaSquare,
    StepByStepTutorial,
    Explainer,
    Tips,
    ProductMarketing,
    PromoAd,
  };

  const getPresentation = (transition: TransitionConfig) => {
    switch (transition.type) {
      case 'fade': return fade();
      case 'slide': return slide({ direction: transition.direction || 'from-right' });
      case 'wipe': return wipe({ direction: transition.direction || 'from-left' });
      case 'flip': return flip({ direction: transition.direction || 'from-left' });
      case 'clock-wipe': return clockWipe();
      default: return fade();
    }
  };

  export const SceneComposer: React.FC<SceneComposerProps> = ({ scenes, transitions }) => {
    return (
      <TransitionSeries>
        {scenes.map((scene, i) => (
          <React.Fragment key={scene.id}>
            <TransitionSeries.Sequence
              durationInFrames={scene.durationInFrames}
              premountFor={30}
            >
              <SceneRenderer
                compositionId={scene.compositionId}
                props={scene.props}
              />
            </TransitionSeries.Sequence>
            {transitions[i] && transitions[i].type !== 'none' && (
              <TransitionSeries.Transition
                presentation={getPresentation(transitions[i])}
                timing={linearTiming({
                  durationInFrames: transitions[i].durationInFrames
                })}
              />
            )}
          </React.Fragment>
        ))}
      </TransitionSeries>
    );
  };


5.2 calculateMetadata
--------------------------------------------------------------------------------

动态计算总时长::

  const calculateSceneComposerMetadata: CalculateMetadataFunction<SceneComposerProps> =
    async ({ props }) => {
      const scenesTotal = props.scenes.reduce((sum, s) => sum + s.durationInFrames, 0);
      const transitionsTotal = props.transitions
        .filter(t => t.type !== 'none')
        .reduce((sum, t) => sum + t.durationInFrames, 0);

      return {
        durationInFrames: scenesTotal - transitionsTotal,
      };
    };


六、MCP 工具重构
================================================================================

6.1 新增/修改的工具
--------------------------------------------------------------------------------

video_list_available_compositions（新增）::

  描述: 列出所有可用的内置组合及其 props schema
  输入: 无
  输出: [{ id, name, description, category, propsSchema, defaultProps }]

video_add_scene（新增，替代 video_add_composition）::

  描述: 向项目添加一个场景片段
  输入: {
    workspacePath, projectId,
    compositionId,        // 内置组合 ID
    durationInFrames,     // 时长（帧数）
    props,                // 组合参数
    insertAt?             // 插入位置（默认末尾）
  }
  输出: { sceneId, projectState }

video_update_scene（新增，替代 video_update_composition）::

  描述: 更新场景片段的参数
  输入: {
    workspacePath, projectId, sceneId,
    props?,               // 更新参数
    durationInFrames?,    // 更新时长
    compositionId?        // 更换组合类型
  }

video_remove_scene（新增，替代 video_remove_composition）::

  描述: 移除场景片段
  输入: { workspacePath, projectId, sceneId }

video_reorder_scenes（新增）::

  描述: 重新排列场景顺序
  输入: { workspacePath, projectId, sceneIds[] }

video_set_transitions（新增）::

  描述: 设置场景间的过渡效果
  输入: {
    workspacePath, projectId,
    transitions: [{ type, durationInFrames, direction? }]
  }

video_get_project_preview（新增）::

  描述: 获取项目预览数据（供前端 Player 渲染）
  输入: { workspacePath, projectId }
  输出: { scenes[], transitions[], config }


6.2 保留的工具
--------------------------------------------------------------------------------

- video_create_project — 创建项目（修改：不再创建 compositions，改为 scenes）
- video_list_projects — 列出项目
- video_get_project — 获取项目详情
- video_delete_project — 删除项目
- video_add_asset — 添加素材
- video_remove_asset — 移除素材
- video_list_templates — 列出模板（修改：返回 compositionId 而非 compositionCode）
- video_create_from_template — 从模板创建（修改：生成 scenes 而非 compositions）
- video_render — 渲染导出（修改：使用统一入口）
- video_get_render_status — 渲染状态


6.3 删除的工具
--------------------------------------------------------------------------------

- video_add_composition — 被 video_add_scene 替代
- video_update_composition — 被 video_update_scene 替代
- video_remove_composition — 被 video_remove_scene 替代
- video_preview_start — 预览改为前端 Player 直接渲染
- video_preview_stop — 同上


七、渲染引擎修复
================================================================================

7.1 RenderEngine 修改
--------------------------------------------------------------------------------

核心变化：始终使用 packages/video/src/Root.tsx 作为入口::

  // packages/video/src/mcp-server/services/render-engine.ts

  class RenderEngine {
    // 旧: 在项目目录找入口
    // findEntryPoint(projectPath) → projectPath/Root.tsx (不存在!)

    // 新: 始终使用 packages/video 作为入口
    private getEntryPoint(): string {
      return path.resolve(__dirname, '../../Root.tsx');
    }

    async render(options: {
      compositionId: string;  // 固定为 "SceneComposer"
      inputProps: SceneComposerProps;  // scenes + transitions
      outputPath: string;
      format: OutputFormat;
      quality: QualityPreset;
    }) {
      const bundled = await bundle({
        entryPoint: this.getEntryPoint(),
        // ...
      });

      await renderMedia({
        composition: await selectComposition({
          serveUrl: bundled,
          id: 'SceneComposer',
          inputProps: options.inputProps,
        }),
        outputLocation: options.outputPath,
        codec: options.format === 'mp4' ? 'h264' : options.format,
        // ...
      });
    }
  }


7.2 RenderWorker 修改
--------------------------------------------------------------------------------

修改 apps/electron/src/main/video/render-worker.ts::

  - 渲染脚本路径改为 packages/video/ 目录
  - 传入 compositionId="SceneComposer" + inputProps={scenes, transitions}
  - 不再依赖用户项目目录的 Remotion 环境


八、前端预览修复
================================================================================

8.1 VideoPreview 修改
--------------------------------------------------------------------------------

修改 apps/electron/src/renderer/components/video/VideoPreview.tsx::

  // 旧: 只支持 4 个组合的直接渲染
  const COMPOSITION_MAP = { TitleAnimation, Slideshow, ... };

  // 新: 使用 SceneComposer 渲染整个项目
  import { SceneComposer } from '@sprouty-ai/video/compositions';

  const VideoPreview = ({ project }) => {
    const composerProps = {
      scenes: project.scenes,
      transitions: project.transitions,
    };

    return (
      <Player
        component={SceneComposer}
        inputProps={composerProps}
        durationInFrames={calculateTotalDuration(project)}
        compositionWidth={project.config.width}
        compositionHeight={project.config.height}
        fps={project.config.fps}
        // ...
      />
    );
  };


九、工作区路径集成
================================================================================

9.1 路径规则
--------------------------------------------------------------------------------

有 creator.db 时::

  工作区/{项目名}/{序号-内容标题}/视频/
  例: ~/saas/我的自媒体/同遇APP官方推广/01-产品介绍视频/视频/
      ├── project.json
      ├── 素材/
      └── 输出/

无 creator.db 时::

  工作区/视频创作/{项目名}/
  例: ~/saas/我的自媒体/视频创作/智小芽推广视频/
      ├── project.json
      ├── 素材/
      └── 输出/


9.2 路径解析
--------------------------------------------------------------------------------

修改 packages/video/src/mcp-server/utils/paths.ts::

  async function resolveVideoProjectPath(
    workspacePath: string,
    projectName: string,
    options?: {
      contentTitle?: string;
      contentIndex?: number;
    }
  ): Promise<string> {
    const dbPath = path.join(workspacePath, '.sprouty-ai', 'db', 'creator.db');

    if (await fileExists(dbPath)) {
      // 从 DB 获取当前活跃项目
      const project = await getActiveProject(dbPath);
      const index = options?.contentIndex ?? await getNextContentIndex(dbPath, project.id);
      const title = options?.contentTitle ?? projectName;
      return path.join(workspacePath, project.name, `${index}-${title}`, '视频');
    }

    return path.join(workspacePath, '视频创作', projectName);
  }


十、Agent 创作示例
================================================================================

10.1 完整创作流程
--------------------------------------------------------------------------------

::

  用户："帮我做一个智小芽 APP 的推广视频"

  Agent 调用 MCP：

  1. video_list_available_compositions
     → 获取所有可用组合及其参数说明

  2. video_create_project("智小芽推广视频", {
       width: 1080, height: 1920, fps: 30
     })
     → 创建项目

  3. video_add_scene(projectId, {
       compositionId: "TitleAnimation",
       durationInFrames: 90,  // 3秒
       props: {
         title: "智小芽",
         subtitle: "你的 AI 学习伙伴",
         animationStyle: "spring",
         colors: { primary: "#6366f1", background: "#0f0f23", text: "#ffffff" }
       }
     })

  4. video_add_scene(projectId, {
       compositionId: "ProductShowcase",
       durationInFrames: 150,  // 5秒
       props: {
         productName: "智小芽 AI 助手",
         layout: "features-grid",
         features: [
           { title: "智能对话", description: "自然语言交互", icon: "💬" },
           { title: "学习辅导", description: "个性化学习方案", icon: "📚" },
           { title: "知识管理", description: "高效整理笔记", icon: "🗂️" }
         ]
       }
     })

  5. video_add_scene(projectId, {
       compositionId: "DataVisualization",
       durationInFrames: 120,  // 4秒
       props: {
         chartType: "bar",
         title: "用户增长趋势",
         data: [
           { label: "Q1", value: 1200 },
           { label: "Q2", value: 3500 },
           { label: "Q3", value: 8900 },
           { label: "Q4", value: 15000 }
         ]
       }
     })

  6. video_add_scene(projectId, {
       compositionId: "SocialMediaVertical",
       durationInFrames: 90,  // 3秒
       props: {
         title: "立即下载",
         subtitle: "开启智能学习之旅",
         cta: { text: "扫码下载" }
       }
     })

  7. video_set_transitions(projectId, [
       { type: "fade", durationInFrames: 15 },
       { type: "slide", durationInFrames: 20, direction: "from-right" },
       { type: "fade", durationInFrames: 15 }
     ])

  → 前端收到项目更新，SceneComposer 用 TransitionSeries 渲染预览
  → 总时长: (90+150+120+90) - (15+20+15) = 400 帧 ≈ 13.3秒


10.2 迭代修改
--------------------------------------------------------------------------------

::

  用户："标题改大一点，第二个场景的颜色换成蓝色"

  Agent 调用 MCP：

  1. video_update_scene(projectId, scene1Id, {
       props: { ...existingProps, fontSize: 80 }
     })

  2. video_update_scene(projectId, scene2Id, {
       props: { ...existingProps, colors: { primary: "#3b82f6", ... } }
     })

  → 前端预览实时刷新


十一、实施阶段
================================================================================

阶段 1：扩展内置组合库
--------------------------------------------------------------------------------

- 新建 7 个组合文件（从模板代码提取）
- 为每个组合添加 Zod Schema
- 修复 CSS transition 违规
- 更新 compositions/index.ts 导出
- 更新 Root.tsx 注册

阶段 2：新增 SceneComposer 组件
--------------------------------------------------------------------------------

- 新建 SceneComposer.tsx（TransitionSeries 编排）
- 新建 SceneRenderer.tsx（组合 ID → 组件映射）
- 添加 calculateMetadata 动态时长计算
- 安装 @remotion/transitions 依赖

阶段 3：重构数据模型
--------------------------------------------------------------------------------

- 修改 types.ts（Scene 替代 Composition，新增 Transition）
- 修改 project-store.ts（存储 scenes + transitions）
- 修改模板定义（compositionId 替代 compositionCode）

阶段 4：重构 MCP 工具
--------------------------------------------------------------------------------

- 新增 scene 相关工具（add/update/remove/reorder）
- 新增 transition 工具
- 新增 list_available_compositions 工具
- 删除旧的 composition 工具
- 删除 preview start/stop 工具

阶段 5：修复渲染引擎
--------------------------------------------------------------------------------

- 修改 RenderEngine（统一入口 Root.tsx）
- 修改 RenderWorker（传入 SceneComposer props）
- 修改渲染脚本路径

阶段 6：修复前端预览
--------------------------------------------------------------------------------

- 修改 VideoPreview（使用 SceneComposer）
- 修改 VideoEditor（场景列表 UI）
- 修改 VideoExport（传入正确参数）

阶段 7：工作区路径集成
--------------------------------------------------------------------------------

- 修改 paths.ts（支持 creator.db 路径模式）
- 修改 project-store.ts（路径解析）

阶段 8：清理
--------------------------------------------------------------------------------

- 删除模板中的 compositionCode 字符串
- 删除 preview-server.ts
- 删除 renderer/index.ts 空占位
- 清理无用的导入和类型


十二、风险评估
================================================================================

| 风险                                | 级别 | 缓解措施                                          |
|-------------------------------------|------|---------------------------------------------------|
| 内置组合不够灵活                    | 中   | 每个组合提供丰富 props，后续持续新增组合           |
| @remotion/renderer Windows 兼容性   | 中   | 渲染在 Electron 主进程，依赖随 APP 打包分发        |
| 现有 project.json 数据迁移          | 低   | 添加迁移逻辑，旧 code 字段映射到 compositionId    |
| @remotion/transitions 包体积        | 低   | 按需导入，tree-shaking 有效                        |
| TransitionSeries 渲染性能           | 低   | premountFor 预加载 + 场景数量通常 < 10             |


十三、后续增强（V2）
================================================================================

- 运行时编译：esbuild 编译 agent 写的 TSX → 高级用户的完全自由创作
- 音频轨道：背景音乐、配音、音效
- 字幕系统：@remotion/captions 集成
- 视频嵌入：在场景中嵌入用户上传的视频片段
- AI 配音：ElevenLabs TTS 集成
- 更多过渡效果：自定义过渡动画
