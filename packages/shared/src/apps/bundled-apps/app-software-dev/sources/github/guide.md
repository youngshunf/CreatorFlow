# GitHub

访问 GitHub 代码仓库、Issues、Pull Requests 和项目管理功能。

## Scope

- 浏览和搜索代码仓库
- 管理 Issues 和 Pull Requests
- 查看提交历史和代码变更
- 项目看板和里程碑管理

## Guidelines

### API 认证
1. 访问 [GitHub Settings > Developer settings > Personal access tokens](https://github.com/settings/tokens)
2. 生成 Personal Access Token (classic 或 fine-grained)
3. 选择需要的权限范围 (repo, read:user 等)
4. 将 token 作为 Bearer 认证密钥

### 常用端点
- `GET /user` - 获取当前用户信息
- `GET /user/repos` - 列出用户仓库
- `GET /repos/{owner}/{repo}` - 获取仓库详情
- `GET /repos/{owner}/{repo}/issues` - 列出 Issues
- `GET /repos/{owner}/{repo}/pulls` - 列出 Pull Requests
- `GET /search/repositories?q={query}` - 搜索仓库

### 使用建议
1. 使用 fine-grained token 限制访问范围
2. 定期轮换 token 保证安全
3. 注意 API 速率限制 (每小时 5000 次)

## Context

GitHub API 提供对代码仓库和项目管理的完整访问。

速率限制：
- 认证用户: 5000 次/小时
- 未认证: 60 次/小时

## API Notes

响应格式为 JSON，分页使用 `page` 和 `per_page` 参数。
Link header 包含分页导航信息。
