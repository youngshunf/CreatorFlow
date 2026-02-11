/**
 * social-publisher — 人类行为模拟
 *
 * 模拟真实用户的输入、鼠标移动、点击和滚动行为，
 * 通过随机延迟和贝塞尔曲线轨迹降低自动化检测风险。
 */

import type { Page } from 'playwright';

// ============================================================
// HumanBehavior
// ============================================================

export class HumanBehavior {
  /**
   * 模拟人类键盘输入
   *
   * 逐字符输入，每个字符之间添加 50-150ms 的随机延迟，
   * 模拟真实打字节奏。
   */
  async typeText(page: Page, selector: string, text: string): Promise<void> {
    await page.click(selector);
    await this.randomWait(200, 500);

    for (const char of text) {
      await page.keyboard.type(char, { delay: 0 });
      await sleep(randomInt(50, 150));
    }
  }

  /**
   * 贝塞尔曲线鼠标移动
   *
   * 从当前位置移动到目标坐标，使用三次贝塞尔曲线生成平滑轨迹，
   * 每步之间添加微小延迟模拟真实鼠标运动。
   */
  async moveTo(page: Page, targetX: number, targetY: number): Promise<void> {
    // 获取当前鼠标位置（默认从页面中心开始）
    const startX = randomInt(100, 300);
    const startY = randomInt(100, 300);

    // 生成贝塞尔曲线控制点
    const cp1x = startX + (targetX - startX) * 0.25 + randomInt(-30, 30);
    const cp1y = startY + (targetY - startY) * 0.1 + randomInt(-30, 30);
    const cp2x = startX + (targetX - startX) * 0.75 + randomInt(-20, 20);
    const cp2y = startY + (targetY - startY) * 0.9 + randomInt(-20, 20);

    // 沿曲线移动（分 20-30 步）
    const steps = randomInt(20, 30);
    for (let i = 0; i <= steps; i++) {
      const t = i / steps;
      const x = cubicBezier(t, startX, cp1x, cp2x, targetX);
      const y = cubicBezier(t, startY, cp1y, cp2y, targetY);
      await page.mouse.move(x, y);
      await sleep(randomInt(5, 15));
    }
  }

  /**
   * 模拟人类点击
   *
   * 先将鼠标移动到元素附近（添加 ±3px 随机偏移），
   * 然后执行点击操作。
   */
  async click(page: Page, selector: string): Promise<void> {
    const element = await page.waitForSelector(selector, { timeout: 10000 });
    if (!element) throw new Error(`元素未找到: ${selector}`);

    const box = await element.boundingBox();
    if (!box) throw new Error(`无法获取元素边界: ${selector}`);

    // 计算点击位置：元素中心 + 随机偏移
    const x = box.x + box.width / 2 + randomInt(-3, 3);
    const y = box.y + box.height / 2 + randomInt(-3, 3);

    // 先移动鼠标到目标位置
    await this.moveTo(page, x, y);
    await sleep(randomInt(50, 150));

    // 执行点击
    await page.mouse.click(x, y);
  }

  /**
   * 随机等待
   *
   * 在指定范围内随机等待一段时间，模拟用户思考或浏览的间隔。
   */
  async randomWait(minMs: number = 2000, maxMs: number = 5000): Promise<void> {
    await sleep(randomInt(minMs, maxMs));
  }

  /**
   * 模拟页面滚动
   *
   * 使用鼠标滚轮事件分多步滚动，每步之间添加随机延迟，
   * 模拟真实用户的滚动行为。
   *
   * @param direction - 滚动方向：'up' 向上，'down' 向下
   * @param distance - 滚动总距离（像素），默认随机 300-800
   */
  async scroll(
    page: Page,
    direction: 'up' | 'down',
    distance?: number,
  ): Promise<void> {
    const totalDistance = distance ?? randomInt(300, 800);
    const steps = randomInt(5, 12);
    const stepDistance = totalDistance / steps;
    const sign = direction === 'down' ? 1 : -1;

    for (let i = 0; i < steps; i++) {
      // 每步滚动距离添加 ±20% 的随机波动
      const jitter = stepDistance * (0.8 + Math.random() * 0.4);
      await page.mouse.wheel(0, sign * jitter);
      await sleep(randomInt(50, 200));
    }
  }
}

// ============================================================
// 工具函数
// ============================================================

/** 三次贝塞尔曲线插值 */
function cubicBezier(t: number, p0: number, p1: number, p2: number, p3: number): number {
  const u = 1 - t;
  return u * u * u * p0 + 3 * u * u * t * p1 + 3 * u * t * t * p2 + t * t * t * p3;
}

/** 生成 [min, max] 范围内的随机整数 */
function randomInt(min: number, max: number): number {
  return Math.floor(Math.random() * (max - min + 1)) + min;
}

/** 异步等待指定毫秒数 */
function sleep(ms: number): Promise<void> {
  return new Promise(resolve => setTimeout(resolve, ms));
}
