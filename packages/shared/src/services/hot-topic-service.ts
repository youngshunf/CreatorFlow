/**
 * 热榜服务 — 从 NewsNow API 获取热榜数据并存入本地数据库
 */

import type { CreatorMediaDB } from '../db/connection.ts';
import type { HotTopic, CreateHotTopic, HotTopicFetchSource } from '../db/types.ts';
import * as hotTopicsRepo from '../db/repositories/hot-topics.ts';
import type { HotTopicFilters } from '../db/repositories/hot-topics.ts';

/** NewsNow API 默认地址 */
const DEFAULT_API_BASE = 'https://newsnow.busiyi.world/api/s';

/** 浏览器 UA，避免被 Cloudflare 拦截 */
const BROWSER_UA = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36';

/** NewsNow API 返回的单条数据 */
interface NewsNowItem {
  id?: string;
  title: string;
  url?: string;
  mobileUrl?: string;
  extra?: {
    hover?: number;
    heat?: number;
    icon?: unknown;
  };
}

/** NewsNow 单平台响应（GET /api/s?id=xxx） */
interface NewsNowSingleResponse {
  status?: string;
  id?: string;
  updatedTime?: number;
  items?: NewsNowItem[];
}

/** NewsNow 批量响应（POST /api/s/entire） */
interface NewsNowEntireResponse {
  [sourceId: string]: NewsNowSingleResponse;
}

/** 热榜平台配置 */
interface HotTopicPlatformConfig {
  id: string;
  name: string;
}

/** 支持的热榜平台列表（对齐 NewsNow 实际可用 sources） */
export const HOT_TOPIC_PLATFORMS: HotTopicPlatformConfig[] = [
  { id: 'weibo', name: '微博热搜' },
  { id: 'zhihu', name: '知乎热榜' },
  { id: 'bilibili-hot-search', name: 'B站热搜' },
  { id: 'douyin', name: '抖音热搜' },
  { id: 'baidu', name: '百度热搜' },
  { id: 'toutiao', name: '今日头条' },
  { id: 'juejin', name: '掘金' },
  { id: '36kr-renqi', name: '36氪' },
  { id: 'thepaper', name: '澎湃新闻' },
  { id: 'wallstreetcn-hot', name: '华尔街见闻' },
  { id: 'cls-hot', name: '财联社' },
  { id: 'douban', name: '豆瓣' },
  { id: 'tieba', name: '贴吧' },
  { id: 'sspai', name: '少数派' },
  { id: 'iqiyi-hot-ranklist', name: '爱奇艺' },
  { id: 'tencent-hot', name: '腾讯新闻' },
  { id: 'hupu', name: '虎扑' },
  { id: 'coolapk', name: '酷安' },
  { id: 'github-trending-today', name: 'GitHub Trending' },
  { id: 'hackernews', name: 'Hacker News' },
  { id: 'producthunt', name: 'Product Hunt' },
];

/** 服务配置 */
export interface HotTopicServiceConfig {
  apiBase?: string;
}

/** 生成短 ID */
function generateId(): string {
  return Math.random().toString(36).substring(2, 15) + Math.random().toString(36).substring(2, 15);
}

/** 获取当天日期 YYYY-MM-DD */
function getTodayDate(): string {
  return new Date().toISOString().split('T')[0]!;
}

/** 热榜服务 */
export class HotTopicService {
  private db: CreatorMediaDB;
  private apiBase: string;

  constructor(db: CreatorMediaDB, config?: HotTopicServiceConfig) {
    this.db = db;
    this.apiBase = config?.apiBase ?? DEFAULT_API_BASE;
  }

  /** 判断数据来源 */
  private getFetchSource(): HotTopicFetchSource {
    return this.apiBase === DEFAULT_API_BASE ? 'newsnow' : 'self-hosted';
  }

  /** 公共请求头 */
  private get headers(): Record<string, string> {
    return {
      'User-Agent': BROWSER_UA,
      'Accept': 'application/json',
      'Referer': this.apiBase.replace('/api/s', '/'),
    };
  }

  /**
   * 批量拉取 — POST /api/s/entire
   * 一次请求获取所有平台数据，效率更高
   */
  private async fetchEntire(platformIds: string[]): Promise<Map<string, NewsNowItem[]>> {
    const url = `${this.apiBase}/entire`;
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 30000);
    const result = new Map<string, NewsNowItem[]>();

    try {
      const response = await fetch(url, {
        method: 'POST',
        signal: controller.signal,
        headers: { ...this.headers, 'Content-Type': 'application/json' },
        body: JSON.stringify({ sources: platformIds }),
      });
      if (!response.ok) return result;
      const json = await response.json() as NewsNowEntireResponse;
      for (const [sourceId, data] of Object.entries(json)) {
        if (data?.items?.length) {
          result.set(sourceId, data.items);
        }
      }
    } catch {
      // 批量接口失败，调用方会降级到逐个拉取
    } finally {
      clearTimeout(timeout);
    }
    return result;
  }

  /** 从 NewsNow API 获取单个平台的热榜（降级方案） */
  private async fetchPlatform(platformId: string): Promise<NewsNowItem[]> {
    const url = `${this.apiBase}?id=${encodeURIComponent(platformId)}`;
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 10000);

    try {
      const response = await fetch(url, {
        signal: controller.signal,
        headers: this.headers,
      });
      if (!response.ok) return [];
      const json = await response.json() as NewsNowSingleResponse;
      return json.items ?? [];
    } catch {
      return [];
    } finally {
      clearTimeout(timeout);
    }
  }

  /** 将 NewsNowItem 列表转为 CreateHotTopic 列表 */
  private toCreateTopics(
    platformId: string,
    platformName: string,
    items: NewsNowItem[],
    fetchSource: HotTopicFetchSource,
    fetchedAt: string,
    batchDate: string,
  ): CreateHotTopic[] {
    const topics: CreateHotTopic[] = [];
    for (let i = 0; i < items.length; i++) {
      const item = items[i]!;
      if (!item.title) continue;
      topics.push({
        id: generateId(),
        platform_id: platformId,
        platform_name: platformName,
        title: item.title,
        url: item.url || item.mobileUrl || null,
        rank: i + 1,
        heat_score: item.extra?.hover ?? item.extra?.heat ?? null,
        fetch_source: fetchSource,
        fetched_at: fetchedAt,
        batch_date: batchDate,
      });
    }
    return topics;
  }

  /**
   * 拉取热榜数据并存入本地数据库
   * 优先使用批量接口，失败则降级逐个拉取
   * @param platforms 要拉取的平台 ID 列表，不传则拉取所有平台
   */
  async fetchHotTopics(platforms?: string[]): Promise<{ count: number; source: HotTopicFetchSource }> {
    const targetPlatforms = platforms
      ? HOT_TOPIC_PLATFORMS.filter(p => platforms.includes(p.id))
      : HOT_TOPIC_PLATFORMS;

    const batchDate = getTodayDate();
    const fetchedAt = new Date().toISOString();
    const fetchSource = this.getFetchSource();
    const allTopics: CreateHotTopic[] = [];
    const platformIds = targetPlatforms.map(p => p.id);

    // 优先批量拉取
    const batchResult = await this.fetchEntire(platformIds);

    // 处理批量结果 + 降级逐个拉取失败的平台
    for (const platform of targetPlatforms) {
      let items = batchResult.get(platform.id);
      if (!items?.length) {
        // 批量接口没拿到该平台数据，降级单独拉取
        items = await this.fetchPlatform(platform.id);
      }
      if (items.length > 0) {
        allTopics.push(...this.toCreateTopics(
          platform.id, platform.name, items, fetchSource, fetchedAt, batchDate,
        ));
      }
    }

    // 批量写入数据库
    if (allTopics.length > 0) {
      hotTopicsRepo.batchInsert(this.db, allTopics);
    }

    return { count: allTopics.length, source: fetchSource };
  }

  /** 查询本地热榜数据 */
  getHotTopics(filters?: HotTopicFilters): HotTopic[] {
    return hotTopicsRepo.listByFilters(this.db, filters);
  }

  /** 获取最新批次信息 */
  getLatestBatch() {
    return hotTopicsRepo.getLatestBatch(this.db);
  }

  /** 清理旧数据 */
  cleanOldData(keepDays: number = 7): number {
    return hotTopicsRepo.cleanOldData(this.db, keepDays);
  }
}
