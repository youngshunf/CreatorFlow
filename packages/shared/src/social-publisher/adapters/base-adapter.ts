/**
 * 平台适配器基类
 *
 * 提供通用的发布/抓取流程，子类只需实现平台特定操作
 */

import type { Page, BrowserContext } from 'playwright';
import type {
  SocialPlatform,
  PlatformAdapter,
  PublishContent,
  PublishResult,
  MetricsResult,
  CommentResult,
} from '../types.ts';
import { HumanBehavior } from '../browser/human-behavior.ts';

export abstract class BaseAdapter implements PlatformAdapter {
  abstract platform: SocialPlatform;
  protected humanBehavior = new HumanBehavior();

  abstract getLoginUrl(): string;
  abstract getPublishUrl(): string;
  abstract detectLoginSuccess(page: Page): Promise<boolean>;
  abstract executePublish(page: Page, content: PublishContent): Promise<PublishResult>;
  abstract scrapeMetrics(page: Page): Promise<MetricsResult>;
  abstract scrapeComments(page: Page, limit: number): Promise<CommentResult[]>;

  /** 回复评论 — 默认未实现，子类可覆盖 */
  async replyComment(
    _page: Page,
    _commentId: string,
    _replyText: string,
  ): Promise<{ success: boolean; error?: string }> {
    return { success: false, error: '该平台暂不支持自动回复评论' };
  }

  /**
   * 通用发布流程
   * 打开新页面 → 随机等待 → 导航到发布页 → 执行平台特定发布
   */
  async publish(context: BrowserContext, content: PublishContent): Promise<PublishResult> {
    const page = await context.newPage();
    try {
      await this.humanBehavior.randomWait(1000, 3000);
      await page.goto(this.getPublishUrl(), { waitUntil: 'domcontentloaded' });
      await this.humanBehavior.randomWait(2000, 4000);
      return await this.executePublish(page, content);
    } catch (error) {
      return {
        success: false,
        error: `发布失败: ${error instanceof Error ? error.message : String(error)}`,
      };
    } finally {
      await page.close().catch(() => {});
    }
  }

  /**
   * 通用指标抓取流程
   */
  async scrape(context: BrowserContext, publishUrl: string): Promise<MetricsResult> {
    const page = await context.newPage();
    try {
      await this.humanBehavior.randomWait(1000, 2000);
      await page.goto(publishUrl, { waitUntil: 'domcontentloaded' });
      await this.humanBehavior.randomWait(2000, 4000);
      return await this.scrapeMetrics(page);
    } catch (error) {
      return {
        success: false,
        error: `抓取指标失败: ${error instanceof Error ? error.message : String(error)}`,
      };
    } finally {
      await page.close().catch(() => {});
    }
  }

  /**
   * 通用评论抓取流程
   */
  async fetchComments(context: BrowserContext, publishUrl: string, limit: number): Promise<CommentResult[]> {
    const page = await context.newPage();
    try {
      await this.humanBehavior.randomWait(1000, 2000);
      await page.goto(publishUrl, { waitUntil: 'domcontentloaded' });
      await this.humanBehavior.randomWait(2000, 4000);
      return await this.scrapeComments(page, limit);
    } catch (error) {
      return [];
    } finally {
      await page.close().catch(() => {});
    }
  }
}
