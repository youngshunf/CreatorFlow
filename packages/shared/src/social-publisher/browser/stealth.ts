/**
 * social-publisher — 反检测配置
 *
 * 提供 stealth 插件初始化和额外的反检测脚本注入，
 * 用于绕过平台的自动化检测机制。
 */

// ============================================================
// Stealth 插件
// ============================================================

/**
 * 获取配置好的 stealth 插件实例
 *
 * 使用 puppeteer-extra-plugin-stealth 提供基础反检测能力，
 * 通过动态导入避免硬依赖缺失时的运行时错误。
 */
export async function getStealthPlugin(): Promise<unknown | null> {
  try {
    // puppeteer-extra-plugin-stealth 需要动态导入
    const { default: StealthPlugin } = await import('puppeteer-extra-plugin-stealth');
    return StealthPlugin();
  } catch {
    console.warn('[social-publisher] puppeteer-extra-plugin-stealth 未安装，跳过 stealth 插件');
    return null;
  }
}

// ============================================================
// 额外反检测脚本
// ============================================================

/**
 * 获取需要注入页面的反检测脚本
 *
 * 这些脚本在页面加载前注入，修复常见的自动化检测点：
 * - navigator.webdriver 属性
 * - chrome.runtime 对象
 * - Permissions API 行为
 * - WebGL 调试信息
 * - 插件和语言列表
 */
export function getAntiDetectScripts(): string {
  return `
    // 1. 移除 navigator.webdriver 标记
    Object.defineProperty(navigator, 'webdriver', {
      get: () => undefined,
    });

    // 2. 模拟 chrome.runtime 对象（非自动化环境下存在）
    if (!window.chrome) {
      window.chrome = {};
    }
    if (!window.chrome.runtime) {
      window.chrome.runtime = {
        connect: function() {},
        sendMessage: function() {},
      };
    }

    // 3. 修复 Permissions API — 使 'notifications' 查询返回 'prompt'
    const originalQuery = window.navigator.permissions.query;
    window.navigator.permissions.query = function(parameters) {
      if (parameters.name === 'notifications') {
        return Promise.resolve({ state: Notification.permission });
      }
      return originalQuery.call(this, parameters);
    };

    // 4. 修复 plugins 数组长度（正常浏览器至少有几个插件）
    Object.defineProperty(navigator, 'plugins', {
      get: () => {
        const plugins = [
          { name: 'Chrome PDF Plugin', filename: 'internal-pdf-viewer', description: 'Portable Document Format' },
          { name: 'Chrome PDF Viewer', filename: 'mhjfbmdgcfjbbpaeojofohoefgiehjai', description: '' },
          { name: 'Native Client', filename: 'internal-nacl-plugin', description: '' },
        ];
        plugins.length = 3;
        return plugins;
      },
    });

    // 5. 修复 languages 属性
    Object.defineProperty(navigator, 'languages', {
      get: () => ['zh-CN', 'zh', 'en-US', 'en'],
    });

    // 6. 防止通过 Error.stack 检测自动化框架
    const originalError = Error;
    Error = class extends originalError {
      constructor(message) {
        super(message);
        if (this.stack) {
          this.stack = this.stack.replace(/\\n.*playwright.*|\\n.*puppeteer.*/gi, '');
        }
      }
    };
  `;
}
