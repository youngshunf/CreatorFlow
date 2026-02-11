/**
 * 选题推荐服务 — 基于热榜和账号画像生成个性化选题
 */

import type { CreatorMediaDB } from '../db/connection.ts';
import type {
  HotTopic,
  RecommendedTopic,
  CreateRecommendedTopic,
  AccountProfile,
  Content,
  PipelineMode,
} from '../db/types.ts';
import * as recommendedTopicsRepo from '../db/repositories/recommended-topics.ts';
import * as contentsRepo from '../db/repositories/contents.ts';
import type { RecommendedTopicFilters } from '../db/repositories/recommended-topics.ts';

/** 生成短 ID */
function generateId(): string {
  return Math.random().toString(36).substring(2, 15) + Math.random().toString(36).substring(2, 15);
}

/** 简单的字符串哈希（取前8位） */
function simpleHash(str: string): string {
  let hash = 0;
  for (let i = 0; i < str.length; i++) {
    const char = str.charCodeAt(i);
    hash = ((hash << 5) - hash) + char;
    hash |= 0;
  }
  return Math.abs(hash).toString(16).padStart(8, '0').substring(0, 8);
}

/** 获取当天日期 YYYY-MM-DD */
function getTodayDate(): string {
  return new Date().toISOString().split('T')[0]!;
}

/** Agent 提示词 */
export interface TopicPrompt {
  system: string;
  user: string;
}

/** AI 返回的单条选题原始数据 */
export interface RawTopicFromAI {
  title: string;
  potential_score: number;
  heat_index: number;
  reason?: string;
  keywords?: unknown;
  platform_heat?: unknown;
  heat_sources?: unknown;
  trend?: unknown;
  target_audience?: unknown;
  creative_angles?: unknown;
  content_outline?: unknown;
  format_suggestions?: unknown;
  material_clues?: unknown;
  risk_notes?: unknown;
}

/** 采纳结果 */
export interface AdoptResult {
  topic: RecommendedTopic;
  content: Content;
}

/** 选题推荐服务 */
export class TopicRecommendService {
  private db: CreatorMediaDB;

  constructor(db: CreatorMediaDB) {
    this.db = db;
  }

  /**
   * 构建 SproutyAgent 选题生成提示词
   */
  buildPrompt(profile: AccountProfile, hotTopics: HotTopic[]): TopicPrompt {
    // 解析 JSON 字段
    const contentPillars = profile.content_pillars ? JSON.parse(profile.content_pillars) : [];
    const pillarWeights = profile.pillar_weights ? JSON.parse(profile.pillar_weights) : {};
    const keywords = profile.keywords ? JSON.parse(profile.keywords) : [];

    const system = `你是一个资深自媒体选题分析专家。根据用户的账号画像和当前热点数据，生成个性化的选题推荐。

输出要求：严格输出 JSON 数组，每个选题包含以下字段：
- title: string — 选题标题（具体、有吸引力、包含关键词）
- potential_score: number — 潜力分数（0-100）
- heat_index: number — 热度指数（0-100）
- reason: string — 推荐理由（50-100字，说明为什么适合该账号）
- keywords: Array<{tag: string, weight: number, type: "core"|"long_tail"|"trend"}> — 关键词标签
- platform_heat: Record<string, {score: number, rank: number, mentions: number}> — 平台热度分布
- heat_sources: Array<{platform_id: string, title: string, url: string, rank: number}> — 热度来源
- trend: {direction: "up"|"down"|"flat", level: "weak"|"medium"|"strong", confidence: number} — 趋势信息
- target_audience: {age_range: string, gender_ratio: string, interests: string[], description: string} — 目标受众
- creative_angles: Array<{angle: string, hook: string, description: string}> — 创作角度（1-2个）
- content_outline: Array<{title: string, points: string[]}> — 内容结构要点（2-3段）
- format_suggestions: Array<{type: string, duration: string, platforms: string[]}> — 形式建议
- material_clues: Array<{type: string, description: string, source: string}> — 素材线索
- risk_notes: Array<{level: "low"|"medium"|"high", risk: string, advice: string}> — 风险提示

约束：
- 选题必须匹配用户的内容支柱和目标受众
- format_suggestions 优先推荐用户发布平台对应的内容形式
- 生成 3-5 个选题
- 只输出 JSON 数组，不要输出其他内容`;

    // 构建热榜摘要（限制长度）
    const hotTopicsSummary = hotTopics.slice(0, 50).map((t, i) => ({
      rank: t.rank ?? i + 1,
      title: t.title,
      platform: t.platform_name,
      heat_score: t.heat_score,
      url: t.url,
    }));

    const user = `## 账号画像
- 领域: ${profile.niche}${profile.sub_niche ? ' / ' + profile.sub_niche : ''}
- 人设: ${profile.persona || '未设置'}
- 内容支柱: ${JSON.stringify(contentPillars)}（权重: ${JSON.stringify(pillarWeights)}）
- 目标受众: ${profile.target_audience || '未设置'}
- 语气风格: ${profile.tone || '未设置'}
- 关键词: ${JSON.stringify(keywords)}
- 发布频率: ${profile.posting_frequency || '未设置'}

## 当前热点数据（最新批次）
${JSON.stringify(hotTopicsSummary, null, 2)}

请生成 3-5 个与该账号高度匹配的选题推荐。只输出 JSON 数组。`;

    return { system, user };
  }

  /**
   * 将 AI 返回的原始数据转换为数据库记录
   */
  parseAIResponse(projectId: string, rawTopics: RawTopicFromAI[]): CreateRecommendedTopic[] {
    const batchDate = getTodayDate();

    return rawTopics.map(raw => {
      const titleHash = simpleHash(raw.title.substring(0, 20));
      return {
        id: generateId(),
        project_id: projectId,
        title: raw.title,
        industry_id: null,
        potential_score: raw.potential_score ?? 0,
        heat_index: raw.heat_index ?? 0,
        reason: raw.reason ?? null,
        keywords: raw.keywords ? JSON.stringify(raw.keywords) : null,
        platform_heat: raw.platform_heat ? JSON.stringify(raw.platform_heat) : null,
        heat_sources: raw.heat_sources ? JSON.stringify(raw.heat_sources) : null,
        trend: raw.trend ? JSON.stringify(raw.trend) : null,
        industry_tags: null,
        target_audience: raw.target_audience ? JSON.stringify(raw.target_audience) : null,
        creative_angles: raw.creative_angles ? JSON.stringify(raw.creative_angles) : null,
        content_outline: raw.content_outline ? JSON.stringify(raw.content_outline) : null,
        format_suggestions: raw.format_suggestions ? JSON.stringify(raw.format_suggestions) : null,
        material_clues: raw.material_clues ? JSON.stringify(raw.material_clues) : null,
        risk_notes: raw.risk_notes ? JSON.stringify(raw.risk_notes) : null,
        source_info: JSON.stringify({ batch_date: batchDate, generator: 'sprouty-agent' }),
        batch_date: batchDate,
        source_uid: `${projectId}::${batchDate}::${titleHash}`,
        status: 0 as const,
        content_id: null,
      };
    });
  }

  /**
   * 保存生成的选题到数据库
   */
  saveTopics(topics: CreateRecommendedTopic[]): void {
    if (topics.length > 0) {
      recommendedTopicsRepo.batchCreate(this.db, topics);
    }
  }

  /**
   * 采纳选题 → 创建内容 → 关联
   */
  adoptTopic(topicId: string, projectId: string, pipelineMode: PipelineMode = 'semi-auto'): AdoptResult | null {
    const topic = recommendedTopicsRepo.get(this.db, topicId);
    if (!topic) return null;

    // 在 contents 表创建记录
    const contentId = generateId();
    const content = contentsRepo.createContent(this.db, {
      id: contentId,
      project_id: projectId,
      title: topic.title,
      topic: topic.title,
      topic_source: 'hot_topic',
      source_topic_id: topicId,
      script_path: null,
      status: 'idea',
      content_type: null,
      target_platforms: null,
      pipeline_mode: pipelineMode,
      pipeline_state: null,
      viral_pattern_id: null,
      tags: topic.keywords,
      scheduled_at: null,
      files: null,
      metadata: null,
      review_summary: null,
    });

    // 关联选题和内容
    recommendedTopicsRepo.linkContent(this.db, topicId, contentId);

    // 重新获取更新后的 topic
    const updatedTopic = recommendedTopicsRepo.get(this.db, topicId)!;

    return { topic: updatedTopic, content };
  }

  /** 忽略选题 */
  ignoreTopic(topicId: string): boolean {
    return recommendedTopicsRepo.updateStatus(this.db, topicId, 2);
  }

  /** 批量忽略 */
  batchIgnore(topicIds: string[]): number {
    return recommendedTopicsRepo.batchUpdateStatus(this.db, topicIds, 2);
  }

  /** 查询选题列表 */
  getTopics(projectId: string, filters?: RecommendedTopicFilters): RecommendedTopic[] {
    return recommendedTopicsRepo.listByProject(this.db, projectId, filters);
  }

  /** 获取单条选题详情 */
  getTopic(topicId: string): RecommendedTopic | undefined {
    return recommendedTopicsRepo.get(this.db, topicId);
  }
}
