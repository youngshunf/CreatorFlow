#!/usr/bin/env bun
/**
 * i18n 自动翻译脚本
 * 
 * 功能：
 * 1. 扫描代码中的 t() 和 $t() 调用
 * 2. 提取中文文本
 * 3. 使用 AI API 进行翻译
 * 4. 生成各语言的 JSON 文件
 * 
 * 使用方法：
 *   bun run scripts/i18n-translate.ts [options]
 * 
 * 选项：
 *   --scan       仅扫描并提取文本，不翻译
 *   --translate  执行翻译
 *   --lang=xx    指定目标语言 (默认: en)
 *   --all        翻译所有配置的语言
 *   --api=xxx    指定 API (anthropic/openai, 默认: anthropic)
 * 
 * 环境变量：
 *   ANTHROPIC_API_KEY - Anthropic API 密钥
 *   OPENAI_API_KEY - OpenAI API 密钥
 */

import * as fs from 'fs';
import * as path from 'path';
import { glob } from 'glob';

// 配置
const CONFIG = {
  // 扫描的目录
  scanDirs: [
    'apps/electron/src/renderer',
    'packages/shared/src',
    'packages/ui/src',
  ],
  
  // 文件扩展名
  extensions: ['.ts', '.tsx', '.js', '.jsx'],
  
  // 忽略的目录
  ignoreDirs: ['node_modules', 'dist', 'build', '.git'],
  
  // 语言包目录
  localesDir: 'apps/electron/src/renderer/locales',
  
  // 支持的目标语言
  targetLanguages: ['en', 'ja', 'ko', 'zh-tw', 'es', 'fr'],
  
  // 源语言
  sourceLanguage: 'zh-cn',
};

// 语言名称映射
const LANGUAGE_NAMES: Record<string, string> = {
  'zh-cn': '简体中文',
  'zh-tw': '繁體中文',
  'en': 'English',
  'ja': '日本語',
  'ko': '한국어',
  'es': 'Español',
  'fr': 'Français',
  'de': 'Deutsch',
  'it': 'Italiano',
  'pt': 'Português',
  'ru': 'Русский',
};

/**
 * 从代码中提取 t() 和 $t() 调用
 */
function extractTranslations(code: string): string[] {
  const translations: string[] = [];
  
  // 匹配 t('xxx') 和 $t('xxx', ...)
  // 支持单引号、双引号和反引号
  const patterns = [
    /\bt\s*\(\s*['"`]([^'"`]+)['"`]\s*\)/g,
    /\$t\s*\(\s*['"`]([^'"`]+)['"`]\s*[,)]/g,
  ];
  
  for (const pattern of patterns) {
    let match;
    while ((match = pattern.exec(code)) !== null) {
      const text = match[1];
      if (text && !translations.includes(text)) {
        translations.push(text);
      }
    }
  }
  
  return translations;
}

/**
 * 扫描项目文件
 */
async function scanProject(): Promise<string[]> {
  const projectRoot = process.cwd();
  const allTranslations: Set<string> = new Set();
  
  for (const scanDir of CONFIG.scanDirs) {
    const fullDir = path.join(projectRoot, scanDir);
    
    if (!fs.existsSync(fullDir)) {
      console.log(`目录不存在，跳过: ${scanDir}`);
      continue;
    }
    
    // 构建 glob 模式
    const patterns = CONFIG.extensions.map(ext => 
      path.join(fullDir, '**', `*${ext}`)
    );
    
    for (const pattern of patterns) {
      const files = await glob(pattern, {
        ignore: CONFIG.ignoreDirs.map(dir => `**/${dir}/**`),
      });
      
      for (const file of files) {
        try {
          const content = fs.readFileSync(file, 'utf-8');
          const translations = extractTranslations(content);
          
          for (const t of translations) {
            allTranslations.add(t);
          }
        } catch (error) {
          console.error(`读取文件失败: ${file}`, error);
        }
      }
    }
  }
  
  return Array.from(allTranslations).sort();
}

/**
 * 读取现有的语言包
 */
function loadLocaleFile(lang: string): Map<string, string> {
  const filePath = path.join(process.cwd(), CONFIG.localesDir, `${lang}.json`);
  const map = new Map<string, string>();
  
  if (fs.existsSync(filePath)) {
    try {
      const content = fs.readFileSync(filePath, 'utf-8');
      const data: [string, string][] = JSON.parse(content);
      
      for (const [key, value] of data) {
        map.set(key, value);
      }
    } catch (error) {
      console.error(`读取语言包失败: ${filePath}`, error);
    }
  }
  
  return map;
}

/**
 * 保存语言包
 */
function saveLocaleFile(lang: string, data: Map<string, string>): void {
  const filePath = path.join(process.cwd(), CONFIG.localesDir, `${lang}.json`);
  
  // 转换为数组格式
  const array: [string, string][] = Array.from(data.entries()).sort((a, b) => 
    a[0].localeCompare(b[0], 'zh-CN')
  );
  
  // 格式化 JSON
  const json = JSON.stringify(array, null, 2);
  
  // 确保目录存在
  const dir = path.dirname(filePath);
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir, { recursive: true });
  }
  
  fs.writeFileSync(filePath, json, 'utf-8');
  console.log(`已保存: ${filePath}`);
}

/**
 * 使用 Anthropic API 翻译
 */
async function translateWithAnthropic(
  texts: string[],
  targetLang: string
): Promise<Map<string, string>> {
  const apiKey = process.env.ANTHROPIC_API_KEY;
  
  if (!apiKey) {
    throw new Error('请设置 ANTHROPIC_API_KEY 环境变量');
  }
  
  const langName = LANGUAGE_NAMES[targetLang] || targetLang;
  
  const prompt = `请将以下中文文本翻译成${langName}。

要求：
1. 保持 {xxx} 格式的占位符不变
2. 翻译要准确、自然
3. 对于 UI 文本，要简洁明了
4. 返回格式：每行一个翻译对，使用 ||| 分隔原文和译文

待翻译文本：
${texts.map((t, i) => `${i + 1}. ${t}`).join('\n')}

输出格式示例：
原文1|||译文1
原文2|||译文2

请按此格式输出，每行一对，不要其他内容。`;

  const response = await fetch('http://claude.api.dcfuture.cn/v1/messages', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'x-api-key': apiKey,
      'anthropic-version': '2025-09-29',
    },
    body: JSON.stringify({
      model: 'claude-sonnet-4-5-20250929',
      max_tokens: 4096,
      messages: [
        { role: 'user', content: prompt }
      ],
    }),
  });
  
  if (!response.ok) {
    const error = await response.text();
    throw new Error(`Anthropic API 错误: ${response.status} - ${error}`);
  }
  
  const result = await response.json() as {
    content: Array<{ type: string; text: string }>;
  };
  
  const text = result.content[0]?.text || '';
  console.log('AI 响应:', text.substring(0, 200) + '...');
  
  // 解析文本格式 - 每行一个翻译对，使用 ||| 分隔
  const lines = text
    .split('\n')
    .map(line => line.trim())
    .filter(line => line && line.includes('|||'));
  
  if (lines.length === 0) {
    console.error('原始响应:', text);
    throw new Error('无法解析翻译结果，未找到有效的翻译对（格式：原文|||译文）');
  }
  
  const map = new Map<string, string>();
  
  for (const line of lines) {
    const parts = line.split('|||');
    if (parts.length !== 2) {
      console.warn('跳过无效行:', line);
      continue;
    }
    
    const source = parts[0].trim().replace(/^\d+\.\s*/, ''); // 移除行号
    const target = parts[1].trim();
    
    if (source && target) {
      map.set(source, target);
    }
  }
  
  console.log(`成功解析 ${map.size} 个翻译对`);
  return map;
}

/**
 * 使用 OpenAI API 翻译
 */
async function translateWithOpenAI(
  texts: string[],
  targetLang: string
): Promise<Map<string, string>> {
  const apiKey = process.env.OPENAI_API_KEY;
  
  if (!apiKey) {
    throw new Error('请设置 OPENAI_API_KEY 环境变量');
  }
  
  const langName = LANGUAGE_NAMES[targetLang] || targetLang;
  
  const prompt = `请将以下中文文本翻译成${langName}。

要求：
1. 保持 {xxx} 格式的占位符不变
2. 翻译要准确、自然
3. 对于 UI 文本，要简洁明了
4. 返回格式：每行一个翻译对，使用 ||| 分隔原文和译文

待翻译文本：
${texts.map((t, i) => `${i + 1}. ${t}`).join('\n')}

输出格式示例：
原文1|||译文1
原文2|||译文2

请按此格式输出，每行一对，不要其他内容。`;

  const response = await fetch('https://api.openai.com/v1/chat/completions', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${apiKey}`,
    },
    body: JSON.stringify({
      model: 'gpt-4o',
      messages: [
        { role: 'user', content: prompt }
      ],
      temperature: 0.3,
    }),
  });
  
  if (!response.ok) {
    const error = await response.text();
    throw new Error(`OpenAI API 错误: ${response.status} - ${error}`);
  }
  
  const result = await response.json() as {
    choices: Array<{ message: { content: string } }>;
  };
  
  const text = result.choices[0]?.message?.content || '';
  console.log('AI 响应:', text.substring(0, 200) + '...');
  
  // 解析文本格式 - 每行一个翻译对，使用 ||| 分隔
  const lines = text
    .split('\n')
    .map(line => line.trim())
    .filter(line => line && line.includes('|||'));
  
  if (lines.length === 0) {
    console.error('原始响应:', text);
    throw new Error('无法解析翻译结果，未找到有效的翻译对（格式：原文|||译文）');
  }
  
  const map = new Map<string, string>();
  
  for (const line of lines) {
    const parts = line.split('|||');
    if (parts.length !== 2) {
      console.warn('跳过无效行:', line);
      continue;
    }
    
    const source = parts[0].trim().replace(/^\d+\.\s*/, ''); // 移除行号
    const target = parts[1].trim();
    
    if (source && target) {
      map.set(source, target);
    }
  }
  
  console.log(`成功解析 ${map.size} 个翻译对`);
  return map;
}

/**
 * 翻译函数
 */
async function translate(
  texts: string[],
  targetLang: string,
  api: 'anthropic' | 'openai' = 'anthropic'
): Promise<Map<string, string>> {
  if (api === 'openai') {
    return translateWithOpenAI(texts, targetLang);
  }
  return translateWithAnthropic(texts, targetLang);
}

/**
 * 主函数
 */
async function main() {
  const args = process.argv.slice(2);
  
  const scanOnly = args.includes('--scan');
  const doTranslate = args.includes('--translate');
  const translateAll = args.includes('--all');
  
  // 解析目标语言
  const langArg = args.find(arg => arg.startsWith('--lang='));
  const targetLang = langArg ? langArg.split('=')[1] : 'en';
  
  // 解析 API
  const apiArg = args.find(arg => arg.startsWith('--api='));
  const api = (apiArg?.split('=')[1] || 'anthropic') as 'anthropic' | 'openai';
  
  console.log('🔍 扫描项目文件...\n');
  const allTexts = await scanProject();
  
  console.log(`找到 ${allTexts.length} 个待翻译文本:\n`);
  allTexts.forEach((text, i) => {
    console.log(`  ${i + 1}. ${text}`);
  });
  
  // 更新源语言文件
  const sourceData = loadLocaleFile(CONFIG.sourceLanguage);
  for (const text of allTexts) {
    if (!sourceData.has(text)) {
      sourceData.set(text, '');
    }
  }
  saveLocaleFile(CONFIG.sourceLanguage, sourceData);
  
  if (scanOnly) {
    console.log('\n✅ 扫描完成（仅扫描模式）');
    return;
  }
  
  if (!doTranslate) {
    console.log('\n提示: 使用 --translate 执行翻译');
    return;
  }
  
  // 执行翻译
  const languages = translateAll ? CONFIG.targetLanguages : [targetLang];
  
  for (const lang of languages) {
    console.log(`\n🌍 翻译到 ${LANGUAGE_NAMES[lang] || lang}...`);
    
    // 加载现有翻译
    const existingData = loadLocaleFile(lang);
    
    // 找出需要翻译的文本
    const toTranslate = allTexts.filter(text => {
      const existing = existingData.get(text);
      return !existing || existing === '';
    });
    
    if (toTranslate.length === 0) {
      console.log('  所有文本已翻译');
      continue;
    }
    
    console.log(`  需要翻译 ${toTranslate.length} 条文本`);
    
    try {
      // 分批翻译（每批 50 条）
      const batchSize = 50;
      for (let i = 0; i < toTranslate.length; i += batchSize) {
        const batch = toTranslate.slice(i, i + batchSize);
        console.log(`  翻译第 ${Math.floor(i / batchSize) + 1} 批...`);
        
        const translations = await translate(batch, lang, api);
        
        // 合并翻译结果
        for (const [source, target] of translations) {
          existingData.set(source, target);
        }
      }
      
      // 保存结果
      saveLocaleFile(lang, existingData);
      console.log(`  ✅ ${lang} 翻译完成`);
      
    } catch (error) {
      console.error(`  ❌ 翻译失败:`, error);
    }
  }
  
  console.log('\n✅ 全部完成');
}

// 运行
main().catch(console.error);
