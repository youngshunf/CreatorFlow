/**
 * 项目上下文生成器
 *
 * 从数据库读取项目信息，生成 Markdown 格式的 Agent 上下文
 */

import type { CreatorMediaDB } from './connection.ts';
import type { Project, AccountProfile, PlatformAccount, Competitor } from './types.ts';

/**
 * 生成项目上下文 Markdown
 * 用于注入到 Agent 的系统提示中
 */
export function generateProjectContext(db: CreatorMediaDB, projectId: string): string {
  // 读取项目
  const project = db.prepare<Project>('SELECT * FROM projects WHERE id = ?').get(projectId);
  if (!project) return '';

  // 读取画像
  const profile = db.prepare<AccountProfile>(
    'SELECT * FROM account_profiles WHERE project_id = ?'
  ).get(projectId);

  // 读取平台账号
  const accounts = db.prepare<PlatformAccount>(
    'SELECT * FROM platform_accounts WHERE project_id = ? ORDER BY is_primary DESC'
  ).all(projectId);

  // 读取竞品
  const competitors = db.prepare<Competitor>(
    'SELECT * FROM competitors WHERE project_id = ? ORDER BY created_at ASC'
  ).all(projectId);

  // 构建 Markdown
  const lines: string[] = [];
  lines.push('## 当前项目上下文');
  lines.push('');
  lines.push(`**项目**: ${project.name}`);

  if (profile) {
    lines.push(`**领域**: ${profile.niche}${profile.sub_niche ? ` > ${profile.sub_niche}` : ''}`);
    if (profile.persona) lines.push(`**人设**: ${profile.persona}`);
    if (profile.target_audience) lines.push(`**目标受众**: ${profile.target_audience}`);
    if (profile.tone) lines.push(`**调性**: ${profile.tone}`);

    if (profile.keywords) {
      try {
        const keywords = JSON.parse(profile.keywords) as string[];
        lines.push(`**核心关键词**: ${keywords.join('、')}`);
      } catch { /* 忽略解析错误 */ }
    }

    if (profile.content_pillars) {
      try {
        const pillars = JSON.parse(profile.content_pillars) as string[];
        lines.push(`**内容支柱**: ${pillars.join('、')}`);
      } catch { /* 忽略解析错误 */ }
    }

    if (profile.posting_frequency) lines.push(`**发布频率**: ${profile.posting_frequency}`);

    if (profile.taboo_topics) {
      try {
        const taboos = JSON.parse(profile.taboo_topics) as string[];
        if (taboos.length > 0) lines.push(`**禁忌话题**: ${taboos.join('、')}`);
      } catch { /* 忽略解析错误 */ }
    }
  }

  // 平台账号
  if (accounts.length > 0) {
    lines.push('');
    lines.push('### 平台账号');
    lines.push('| 平台 | 昵称 | 粉丝 | 登录状态 |');
    lines.push('|------|------|------|---------|');
    for (const acc of accounts) {
      const primary = acc.is_primary ? '(主)' : '';
      const followers = formatFollowers(acc.followers);
      const authIcon = acc.auth_status === 'logged_in' ? '✅ 已登录' : '❌ 未登录';
      lines.push(`| ${acc.platform}${primary} | ${acc.nickname ?? '-'} | ${followers} | ${authIcon} |`);
    }
  }

  // 竞品
  if (competitors.length > 0) {
    lines.push('');
    lines.push('### 竞品参考');
    for (const comp of competitors) {
      const followers = comp.follower_count ? formatFollowers(comp.follower_count) + '粉' : '';
      const style = comp.content_style ? `: ${comp.content_style}` : '';
      lines.push(`- **${comp.name}** (${comp.platform}, ${followers})${style}`);
    }
  }

  // 创作要求
  lines.push('');
  lines.push('### 创作要求');
  lines.push('- 所有内容必须符合上述人设和调性');
  lines.push('- 热点筛选优先匹配核心关键词');
  lines.push('- 标题和封面参考竞品的成功模式');
  if (accounts.length > 0) {
    lines.push('- 自动发布仅限已登录的平台账号');
  }

  return lines.join('\n');
}

/** 格式化粉丝数 */
function formatFollowers(count: number): string {
  if (count >= 10000) {
    return `${(count / 10000).toFixed(1)}w`;
  }
  return String(count);
}
