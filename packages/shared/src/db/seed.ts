/**
 * 自媒体创作 APP v2.0 — 默认种子数据
 *
 * 全局爆款模式库（project_id = NULL，所有项目可用）
 */

import { randomUUID } from 'crypto';

/** 生成种子数据 SQL */
export function getSeedSQL(): string {
  const patterns = [
    // hook — 开头钩子
    {
      id: randomUUID(),
      category: 'hook',
      name: '悬念开头法',
      description: '用一个出人意料的结论或问题开头，制造悬念吸引继续阅读',
      template: '你知道吗？{震撼事实}。但真相远比你想的更{形容词}...',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'hook',
      name: '数字冲击法',
      description: '用具体数字开头，制造冲击力和可信度',
      template: '{数字}个{对象}中，只有{小数字}个能{成就}。原因竟然是...',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'hook',
      name: '反常识开头',
      description: '提出与常识相反的观点，引发好奇心',
      template: '别再{常见做法}了！{权威来源}证明，{反常识结论}',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'hook',
      name: '痛点共鸣法',
      description: '直击目标受众的痛点，引发情感共鸣',
      template: '是不是每次{痛点场景}都感觉{负面情绪}？今天教你{解决方案}',
      source: 'manual',
    },

    // structure — 内容结构
    {
      id: randomUUID(),
      category: 'structure',
      name: '3段式测评',
      description: '外观→性能→总结的经典测评结构',
      template: '## 外观设计\n{外观描述}\n\n## 核心性能\n{性能测试}\n\n## 总结推荐\n{总结}',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'structure',
      name: '问题-方案-结果',
      description: '先抛出问题，再给出方案，最后展示结果',
      template: '## 问题\n{痛点描述}\n\n## 方案\n{解决步骤}\n\n## 结果\n{效果展示}',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'structure',
      name: '清单体',
      description: '用编号列表组织内容，适合干货分享',
      template: '# {主题}必知的{N}件事\n\n1. {要点1}\n2. {要点2}\n3. {要点3}\n...',
      source: 'manual',
    },

    // title — 标题公式
    {
      id: randomUUID(),
      category: 'title',
      name: '数字+好处',
      description: '用数字量化好处，提高点击率',
      template: '{数字}个让你{好处}的{方法/技巧/秘诀}',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'title',
      name: '对比反差',
      description: '用前后对比制造反差感',
      template: '从{差状态}到{好状态}，我只做了{一件事}',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'title',
      name: '权威背书',
      description: '借助权威增加可信度',
      template: '{权威人物/机构}推荐的{N}个{对象}，第{X}个太绝了',
      source: 'manual',
    },

    // cta — 行动号召
    {
      id: randomUUID(),
      category: 'cta',
      name: '投票互动',
      description: '用投票引导评论互动',
      template: '你觉得{选项A}还是{选项B}更好？评论区告诉我！',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'cta',
      name: '收藏引导',
      description: '强调实用性引导收藏',
      template: '建议先收藏⭐，{使用场景}的时候一定用得上！',
      source: 'manual',
    },

    // visual — 视觉风格
    {
      id: randomUUID(),
      category: 'visual',
      name: '对比图封面',
      description: '左右对比或前后对比的封面设计',
      template: '封面分为左右两栏，左侧{对比A}，右侧{对比B}，中间用分割线',
      source: 'manual',
    },
    {
      id: randomUUID(),
      category: 'visual',
      name: '大字报封面',
      description: '用大号加粗文字作为封面主体',
      template: '纯色背景 + 居中大字标题（不超过8个字）+ 小字副标题',
      source: 'manual',
    },

    // rhythm — 节奏编排
    {
      id: randomUUID(),
      category: 'rhythm',
      name: '黄金3秒',
      description: '短视频前3秒必须出现核心信息',
      template: '0-3秒: {核心悬念/结论}\n3-15秒: {展开论述}\n15-30秒: {总结+CTA}',
      source: 'manual',
    },
  ];

  const values = patterns.map(p => {
    const examples = JSON.stringify([]);
    const tags = JSON.stringify([p.category]);
    return `('${p.id}', NULL, NULL, '${p.category}', '${p.name}', '${p.description}', '${p.template.replace(/'/g, "''")}', '${examples}', '${p.source}', 0, NULL, '${tags}', CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)`;
  });

  return `
INSERT OR IGNORE INTO viral_patterns (id, project_id, platform, category, name, description, template, examples, source, usage_count, success_rate, tags, created_at, updated_at)
VALUES
${values.join(',\n')};
`;
}
