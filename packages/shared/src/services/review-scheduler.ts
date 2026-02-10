/**
 * 采集调度服务 — 定时轮询 review_tasks 并执行数据采集
 *
 * 发布后按 T+1h, T+6h, T+24h, T+72h, T+168h 五个时间点采集平台数据，
 * 最后一次采集(T+168h)触发反馈闭环。
 */

import type { CreatorMediaDB } from '../db/connection.ts';
import type { ReviewTask } from '../db/types.ts';
import {
  listPendingTasks,
  markExecuting,
  markCompleted,
  markFailed,
  createReviewTask,
  getReviewTask,
  updateReviewTask,
} from '../db/repositories/review-tasks.ts';
import {
  getPublishRecord,
  updatePublishMetrics,
  updatePublishRecord,
} from '../db/repositories/publish-records.ts';

/** 采集时间点配置（相对发布时间的毫秒偏移） */
const REVIEW_SCHEDULE: Array<{ type: '1h' | '6h' | '24h' | '72h' | '168h'; offsetMs: number }> = [
  { type: '1h', offsetMs: 1 * 60 * 60 * 1000 },
  { type: '6h', offsetMs: 6 * 60 * 60 * 1000 },
  { type: '24h', offsetMs: 24 * 60 * 60 * 1000 },
  { type: '72h', offsetMs: 72 * 60 * 60 * 1000 },
  { type: '168h', offsetMs: 168 * 60 * 60 * 1000 },
];

/** 最大重试次数 */
const MAX_RETRY_COUNT = 3;

/** 重试延迟（30 分钟） */
const RETRY_DELAY_MS = 30 * 60 * 1000;

/** 默认轮询间隔（5 分钟） */
const DEFAULT_POLL_INTERVAL_MS = 5 * 60 * 1000;

/** 数据采集回调 — 由外部注入具体的平台采集逻辑 */
export interface MetricsCollector {
  /** 采集平台数据，返回指标快照 */
  collect(platform: string, publishUrl: string): Promise<{
    views: number;
    likes: number;
    comments: number;
    shares: number;
    favorites: number;
    metrics_json?: string;
  }>;
}

/** 反馈闭环回调 — 最后一次采集后触发 */
export interface FeedbackHandler {
  /** 处理反馈闭环（分析数据趋势、更新支柱权重等） */
  processFeedback(publishRecordId: string): Promise<void>;
}

/** ReviewScheduler 配置 */
export interface ReviewSchedulerOptions {
  db: CreatorMediaDB;
  collector: MetricsCollector;
  feedbackHandler?: FeedbackHandler;
  pollIntervalMs?: number;
  batchSize?: number;
}

/** 采集调度服务 */
export class ReviewScheduler {
  private db: CreatorMediaDB;
  private collector: MetricsCollector;
  private feedbackHandler?: FeedbackHandler;
  private pollIntervalMs: number;
  private batchSize: number;
  private timer: ReturnType<typeof setInterval> | null = null;
  private running = false;

  constructor(options: ReviewSchedulerOptions) {
    this.db = options.db;
    this.collector = options.collector;
    this.feedbackHandler = options.feedbackHandler;
    this.pollIntervalMs = options.pollIntervalMs ?? DEFAULT_POLL_INTERVAL_MS;
    this.batchSize = options.batchSize ?? 10;
  }

  /** 启动定时轮询 */
  start(): void {
    if (this.timer) return;
    this.timer = setInterval(() => this.poll(), this.pollIntervalMs);
    // 立即执行一次
    this.poll();
  }

  /** 停止轮询 */
  stop(): void {
    if (this.timer) {
      clearInterval(this.timer);
      this.timer = null;
    }
  }

  /** 查询待执行的 review_tasks，逐个执行 */
  async poll(): Promise<void> {
    if (this.running) return;
    this.running = true;

    try {
      const tasks = listPendingTasks(this.db, this.batchSize);
      for (const task of tasks) {
        await this.executeReviewTask(task);
      }
    } finally {
      this.running = false;
    }
  }

  /** 执行单个采集任务 */
  async executeReviewTask(task: ReviewTask): Promise<void> {
    // 标记为执行中
    markExecuting(this.db, task.id);

    try {
      // 获取发布记录
      const record = getPublishRecord(this.db, task.publish_record_id);
      if (!record || !record.publish_url) {
        markFailed(this.db, task.id, '发布记录不存在或缺少发布链接');
        return;
      }

      // 调用平台数据采集
      const metrics = await this.collector.collect(record.platform, record.publish_url);

      // 保存采集结果快照
      const resultSnapshot = JSON.stringify({
        ...metrics,
        collected_at: new Date().toISOString(),
        review_type: task.review_type,
      });
      markCompleted(this.db, task.id, resultSnapshot);

      // 更新 publish_records 的指标字段
      updatePublishMetrics(this.db, record.id, {
        views: metrics.views,
        likes: metrics.likes,
        comments: metrics.comments,
        shares: metrics.shares,
        favorites: metrics.favorites,
        metrics_json: metrics.metrics_json ?? null,
      });

      // 更新采集计数
      updatePublishRecord(this.db, record.id, {
        review_count: (record.review_count ?? 0) + 1,
      });

      // 如果是最后一次采集(T+168h)，触发反馈闭环
      if (task.review_type === '168h') {
        updatePublishRecord(this.db, record.id, {
          feedback_processed: 1,
        });
        if (this.feedbackHandler) {
          await this.feedbackHandler.processFeedback(record.id);
        }
      }
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : String(error);
      this.handleFailure(task, errorMessage);
    }
  }

  /**
   * 为发布记录创建 5 个采集任务（T+1h, T+6h, T+24h, T+72h, T+168h）
   * @param publishRecordId - 发布记录 ID
   * @param publishedAt - 发布时间（ISO 字符串），默认当前时间
   */
  createReviewTasks(publishRecordId: string, publishedAt?: string): ReviewTask[] {
    const baseTime = publishedAt ? new Date(publishedAt).getTime() : Date.now();
    const tasks: ReviewTask[] = [];

    for (const schedule of REVIEW_SCHEDULE) {
      const scheduledAt = new Date(baseTime + schedule.offsetMs).toISOString();
      const id = `rt_${publishRecordId}_${schedule.type}`;
      const task = createReviewTask(this.db, {
        id,
        publish_record_id: publishRecordId,
        scheduled_at: scheduledAt,
        executed_at: null,
        status: 'pending',
        review_type: schedule.type,
        result_snapshot: null,
        error_message: null,
        retry_count: 0,
      });
      tasks.push(task);
    }

    // 更新发布记录的采集计划
    const scheduleJson = JSON.stringify(
      REVIEW_SCHEDULE.map(s => ({
        type: s.type,
        scheduled_at: new Date(baseTime + s.offsetMs).toISOString(),
      }))
    );
    updatePublishRecord(this.db, publishRecordId, {
      review_schedule: scheduleJson,
      next_review_at: tasks[0]?.scheduled_at,
    });

    return tasks;
  }

  /** 失败处理：retry_count < 3 则延迟 30 分钟重试 */
  private handleFailure(task: ReviewTask, errorMessage: string): void {
    const updated = markFailed(this.db, task.id, errorMessage);
    if (!updated) return;

    if (updated.retry_count < MAX_RETRY_COUNT) {
      // 重新设为 pending，延迟 30 分钟
      const retryAt = new Date(Date.now() + RETRY_DELAY_MS).toISOString();
      updateReviewTask(this.db, task.id, {
        status: 'pending',
        scheduled_at: retryAt,
      });
    }
  }
}

/** 工厂函数 */
export function createReviewScheduler(options: ReviewSchedulerOptions): ReviewScheduler {
  return new ReviewScheduler(options);
}
