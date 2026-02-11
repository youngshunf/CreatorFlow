/**
 * 小红书平台适配器 (P0)
 *
 * 实现小红书的登录检测、内容发布、数据抓取
 */

import type { Page } from 'playwright';
import type {
  SocialPlatform,
  PublishContent,
  PublishResult,
  MetricsResult,
  CommentResult,
} from '../types.ts';
import { BaseAdapter } from './base-adapter.ts';

export class XiaohongshuAdapter extends BaseAdapter {
  platform: SocialPlatform = 'xiaohongshu';

  getLoginUrl(): string {
    return 'https://www.xiaohongshu.com';
  }

  getPublishUrl(): string {
    return 'https://creator.xiaohongshu.com/publish/publish';
  }

  /**
   * 检测小红书登录状态
   * 通过检查用户头像元素或特定 Cookie 判断
   */
  async detectLoginSuccess(page: Page): Promise<boolean> {
    try {
      // 导航到创作者中心
      await page.goto('https://creator.xiaohongshu.com', {
        waitUntil: 'domcontentloaded',
        timeout: 15000,
      });
      await this.humanBehavior.randomWait(2000, 3000);

      // 检查是否存在用户头像或登录后的元素
      const loggedIn = await page.evaluate(() => {
        // 检查是否被重定向到登录页
        if (window.location.href.includes('login')) return false;
        // 检查是否存在用户信息元素
        const avatar = document.querySelector('.user-avatar, .avatar, [class*="avatar"]');
        const userName = document.querySelector('.user-name, [class*="username"], [class*="nick"]');
        return !!(avatar || userName);
      });

      return loggedIn;
    } catch {
      return false;
    }
  }

  /**
   * 执行小红书发布
   */
  async executePublish(page: Page, content: PublishContent): Promise<PublishResult> {
    try {
      // 1. 上传素材
      if (content.media_files.length > 0) {
        const fileInput = await page.waitForSelector('input[type="file"]', { timeout: 10000 });
        if (fileInput) {
          await fileInput.setInputFiles(content.media_files);
          // 等待上传完成
          await this.humanBehavior.randomWait(3000, 6000);
        }
      }

      // 2. 填写标题
      if (content.title) {
        const titleInput = await page.waitForSelector(
          '[placeholder*="标题"], [class*="title"] input, [class*="title"] textarea',
          { timeout: 5000 },
        ).catch(() => null);
        if (titleInput) {
          await titleInput.click();
          await this.humanBehavior.typeText(page, '', content.title);
        }
      }

      // 3. 填写正文
      const bodyInput = await page.waitForSelector(
        '[placeholder*="正文"], [class*="content"] [contenteditable], [class*="editor"] [contenteditable]',
        { timeout: 5000 },
      ).catch(() => null);
      if (bodyInput) {
        await bodyInput.click();
        await this.humanBehavior.randomWait(500, 1000);
        await this.humanBehavior.typeText(page, '', content.text);
      }

      // 4. 添加话题标签
      for (const tag of content.tags) {
        await this.humanBehavior.randomWait(500, 1500);
        // 尝试找到标签输入区域
        const tagInput = await page.waitForSelector(
          '[placeholder*="话题"], [placeholder*="标签"], [class*="tag"] input',
          { timeout: 3000 },
        ).catch(() => null);
        if (tagInput) {
          await tagInput.click();
          await this.humanBehavior.typeText(page, '', `#${tag}`);
          await page.keyboard.press('Enter');
        }
      }

      // 5. 点击发布按钮
      await this.humanBehavior.randomWait(1000, 2000);
      const publishBtn = await page.waitForSelector(
        'button:has-text("发布"), [class*="publish"] button, [class*="submit"] button',
        { timeout: 5000 },
      ).catch(() => null);

      if (!publishBtn) {
        return { success: false, error: '未找到发布按钮' };
      }

      await this.humanBehavior.click(page, 'button:has-text("发布")');

      // 6. 等待发布确认
      await this.humanBehavior.randomWait(3000, 6000);

      // 7. 尝试提取发布后的 URL
      const currentUrl = page.url();
      const postIdMatch = currentUrl.match(/\/(\w{24})/);

      return {
        success: true,
        platform_url: currentUrl,
        platform_post_id: postIdMatch?.[1],
      };
    } catch (error) {
      return {
        success: false,
        error: `小红书发布失败: ${error instanceof Error ? error.message : String(error)}`,
      };
    }
  }

  /**
   * 抓取小红书笔记数据指标
   */
  async scrapeMetrics(page: Page): Promise<MetricsResult> {
    try {
      await this.humanBehavior.randomWait(2000, 3000);

      const metrics = await page.evaluate(() => {
        /** 从文本中提取数字 */
        const extractNumber = (text: string | null): number => {
          if (!text) return 0;
          const cleaned = text.replace(/[,，]/g, '');
          // 处理 "1.2万" 格式
          const wanMatch = cleaned.match(/([\d.]+)\s*万/);
          if (wanMatch) return Math.round(parseFloat(wanMatch[1]!) * 10000);
          const num = parseInt(cleaned.replace(/\D/g, ''), 10);
          return isNaN(num) ? 0 : num;
        };

        // 尝试多种选择器匹配互动数据
        const selectors = {
          likes: ['[class*="like"] span', '[class*="zan"] span', '.like-count'],
          comments: ['[class*="comment"] span', '[class*="chat"] span', '.comment-count'],
          favorites: ['[class*="collect"] span', '[class*="star"] span', '.collect-count'],
          shares: ['[class*="share"] span', '.share-count'],
        };

        const result: Record<string, number> = { views: 0, likes: 0, comments: 0, shares: 0, favorites: 0 };

        for (const [key, sels] of Object.entries(selectors)) {
          for (const sel of sels) {
            const el = document.querySelector(sel);
            if (el?.textContent) {
              result[key] = extractNumber(el.textContent);
              break;
            }
          }
        }

        return result;
      });

      return {
        success: true,
        metrics: {
          views: metrics.views ?? 0,
          likes: metrics.likes ?? 0,
          comments: metrics.comments ?? 0,
          shares: metrics.shares ?? 0,
          favorites: metrics.favorites ?? 0,
          raw_json: JSON.stringify(metrics),
        },
      };
    } catch (error) {
      return {
        success: false,
        error: `抓取小红书指标失败: ${error instanceof Error ? error.message : String(error)}`,
      };
    }
  }

  /**
   * 抓取小红书评论列表
   */
  async scrapeComments(page: Page, limit: number): Promise<CommentResult[]> {
    try {
      // 滚动加载评论
      await this.humanBehavior.scroll(page, 'down', 500);
      await this.humanBehavior.randomWait(2000, 3000);

      const comments = await page.evaluate((maxCount: number) => {
        const results: Array<{
          id: string;
          author: string;
          author_avatar?: string;
          content: string;
          likes: number;
          created_at: string;
        }> = [];

        // 尝试多种评论容器选择器
        const commentEls = document.querySelectorAll(
          '[class*="comment-item"], [class*="commentItem"], .comment-list > div'
        );

        for (let i = 0; i < Math.min(commentEls.length, maxCount); i++) {
          const el = commentEls[i];
          const authorEl = el.querySelector('[class*="author"], [class*="name"], .nickname');
          const avatarEl = el.querySelector('img[class*="avatar"]');
          const contentEl = el.querySelector('[class*="content"], .comment-text');
          const likeEl = el.querySelector('[class*="like"] span');
          const timeEl = el.querySelector('[class*="time"], .date');

          results.push({
            id: el.getAttribute('data-id') || `comment-${i}`,
            author: authorEl?.textContent?.trim() || '未知用户',
            author_avatar: avatarEl?.getAttribute('src') || undefined,
            content: contentEl?.textContent?.trim() || '',
            likes: parseInt(likeEl?.textContent?.replace(/\D/g, '') || '0', 10),
            created_at: timeEl?.textContent?.trim() || new Date().toISOString(),
          });
        }

        return results;
      }, limit);

      return comments;
    } catch {
      return [];
    }
  }

  /**
   * 回复小红书评论
   */
  async replyComment(
    page: Page,
    commentId: string,
    replyText: string,
  ): Promise<{ success: boolean; error?: string }> {
    try {
      // 找到评论的回复按钮
      const replyBtn = await page.waitForSelector(
        `[data-id="${commentId}"] [class*="reply"], [data-id="${commentId}"] button:has-text("回复")`,
        { timeout: 5000 },
      ).catch(() => null);

      if (!replyBtn) {
        return { success: false, error: '未找到回复按钮' };
      }

      await this.humanBehavior.click(page, `[data-id="${commentId}"] [class*="reply"]`);
      await this.humanBehavior.randomWait(500, 1000);

      // 输入回复内容
      const replyInput = await page.waitForSelector(
        '[class*="reply"] textarea, [class*="reply"] [contenteditable]',
        { timeout: 5000 },
      ).catch(() => null);

      if (!replyInput) {
        return { success: false, error: '未找到回复输入框' };
      }

      await replyInput.click();
      await this.humanBehavior.typeText(page, '', replyText);
      await this.humanBehavior.randomWait(500, 1000);

      // 点击发送
      await page.keyboard.press('Enter');
      await this.humanBehavior.randomWait(1000, 2000);

      return { success: true };
    } catch (error) {
      return {
        success: false,
        error: `回复评论失败: ${error instanceof Error ? error.message : String(error)}`,
      };
    }
  }
}
