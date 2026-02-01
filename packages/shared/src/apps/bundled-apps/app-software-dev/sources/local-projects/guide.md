# 本地项目

管理本地代码项目、开发文件和技术文档。

## Scope

- 浏览和搜索项目文件
- 读取代码和配置文件
- 查看项目结构和依赖
- 分析代码和生成文档

## Guidelines

### 推荐目录结构
```
Projects/
├── {项目名}/
│   ├── src/
│   ├── docs/
│   ├── tests/
│   ├── README.md
│   └── package.json
├── snippets/         # 代码片段
├── templates/        # 项目模板
└── learning/         # 学习项目
```

### 项目类型支持
- Node.js / TypeScript 项目
- Python 项目
- Go 项目
- Rust 项目
- 前端项目 (React, Vue, etc.)

### 使用建议
1. 每个项目保持独立目录
2. 使用 Git 进行版本控制
3. 保持 README 文档更新
4. 统一代码风格和结构

## Context

本地项目源用于访问和管理你的开发项目。
启用后可修改路径指向你的项目存储位置。

支持自动识别项目类型：
- package.json → Node.js
- requirements.txt / pyproject.toml → Python
- go.mod → Go
- Cargo.toml → Rust
