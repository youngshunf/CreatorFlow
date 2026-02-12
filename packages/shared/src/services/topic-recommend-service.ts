/**
 * 选题推荐服务 — 精简版（详情存 Markdown 文件）
 */

import type { CreatorMediaDB } from '../db/connection.ts';
import type {
  RecommendedTopic,
  CreateRecommendedTopic,
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
      tags: null,
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
