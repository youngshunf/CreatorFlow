# Sprouty AI 内置应用与技能包规划

## 概述

本文档定义 Sprouty AI 内置应用的行业场景、应用定义和技能包规划。

---

## 一、内容创作类

### 1. 自媒体创作 (app.creator-media) ✅ 已完成
**目标用户**：小红书、抖音、公众号、B站等平台的内容创作者

**核心技能**：
- `material-organize` - 素材整理
- `script-create` - 脚本创作
- `topic-research` - 选题研究
- `content-review` - 内容审核
- `platform-adaptation` - 平台适配
- `viral-content-analysis` - 爆款分析
- `topic-to-publish-flow` - 选题到发布流程
- `data-review` - 数据复盘

---

### 2. 技术写作 (app.tech-writing)
**目标用户**：技术博主、文档工程师、开发者教程作者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `doc-outline` | 文档大纲生成 - 根据主题生成结构化大纲 |
| `code-explain` | 代码解释 - 将代码转化为易懂的文字说明 |
| `api-doc-gen` | API 文档生成 - 从代码生成 API 文档 |
| `tutorial-create` | 教程创作 - 分步骤教程编写 |
| `tech-review` | 技术审核 - 检查技术准确性和表达清晰度 |
| `diagram-suggest` | 图表建议 - 建议合适的图表和流程图 |

**目录结构**：
```
文档/
├── 草稿/
├── 已发布/
└── 归档/
代码示例/
图片资源/
```

---

### 3. 视频脚本 (app.video-script)
**目标用户**：YouTube、B站长视频创作者、短剧编剧

**核心技能**：
| 技能 | 说明 |
|------|------|
| `video-outline` | 视频大纲 - 生成视频结构和分镜 |
| `hook-generator` | 开场钩子 - 生成吸引人的开场 |
| `script-polish` | 脚本润色 - 优化口语化表达 |
| `subtitle-gen` | 字幕生成 - 生成时间轴字幕 |
| `thumbnail-idea` | 封面创意 - 建议封面设计方案 |
| `retention-analysis` | 留存分析 - 分析脚本的观众留存点 |

---

### 4. 播客制作 (app.podcast)
**目标用户**：播客主播、音频内容创作者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `episode-plan` | 节目策划 - 规划单集内容结构 |
| `guest-research` | 嘉宾研究 - 整理嘉宾背景和问题 |
| `show-notes` | 节目笔记 - 生成发布用的节目说明 |
| `transcript-edit` | 转录编辑 - 整理和编辑音频转录 |
| `highlight-extract` | 精彩提取 - 提取可传播的精彩片段 |

---

## 二、开发者工具类

### 5. 软件开发 (app.software-dev)
**目标用户**：全栈开发者、独立开发者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `project-scaffold` | 项目脚手架 - 快速生成项目结构 |
| `code-review` | 代码审查 - 自动化代码审查建议 |
| `test-gen` | 测试生成 - 为代码生成单元测试 |
| `refactor-suggest` | 重构建议 - 识别和建议代码重构 |
| `debug-assist` | 调试助手 - 帮助分析和定位 bug |
| `commit-message` | 提交信息 - 生成规范的 commit message |
| `changelog-gen` | 变更日志 - 生成版本变更日志 |

**目录结构**：
```
src/
tests/
docs/
scripts/
```

---

### 6. 前端开发 (app.frontend-dev)
**目标用户**：前端工程师、UI 开发者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `component-gen` | 组件生成 - 生成 React/Vue 组件 |
| `style-convert` | 样式转换 - CSS/Tailwind/styled 互转 |
| `responsive-check` | 响应式检查 - 检查响应式布局问题 |
| `a11y-audit` | 无障碍审计 - 检查可访问性问题 |
| `bundle-analyze` | 包分析 - 分析和优化打包体积 |
| `design-to-code` | 设计稿转代码 - 从设计稿生成代码 |

---

### 7. DevOps 工程 (app.devops)
**目标用户**：运维工程师、SRE、DevOps 工程师

**核心技能**：
| 技能 | 说明 |
|------|------|
| `dockerfile-gen` | Dockerfile 生成 - 生成优化的 Dockerfile |
| `ci-pipeline` | CI 流水线 - 生成 CI/CD 配置 |
| `k8s-manifest` | K8s 配置 - 生成 Kubernetes 配置 |
| `log-analyze` | 日志分析 - 分析日志定位问题 |
| `monitor-setup` | 监控配置 - 配置监控和告警 |
| `incident-report` | 事故报告 - 生成事故分析报告 |

---

## 三、商业/营销类

### 8. 市场营销 (app.marketing)
**目标用户**：市场经理、品牌经理、营销专员

**核心技能**：
| 技能 | 说明 |
|------|------|
| `campaign-plan` | 活动策划 - 营销活动方案设计 |
| `copy-write` | 文案撰写 - 营销文案创作 |
| `audience-analyze` | 受众分析 - 目标用户画像分析 |
| `competitor-track` | 竞品追踪 - 竞争对手动态监控 |
| `ab-test-design` | AB 测试设计 - 设计营销实验 |
| `roi-calculate` | ROI 计算 - 营销效果评估 |

**目录结构**：
```
活动方案/
素材库/
├── 图片/
├── 视频/
└── 文案/
数据报告/
竞品资料/
```

---

### 9. 电商运营 (app.ecommerce)
**目标用户**：电商店主、运营经理、产品经理

**核心技能**：
| 技能 | 说明 |
|------|------|
| `product-listing` | 商品上架 - 生成商品标题和描述 |
| `price-strategy` | 定价策略 - 价格分析和建议 |
| `review-respond` | 评价回复 - 智能回复客户评价 |
| `inventory-alert` | 库存预警 - 库存分析和补货建议 |
| `promo-plan` | 促销策划 - 促销活动设计 |
| `customer-insight` | 客户洞察 - 分析客户行为数据 |

---

### 10. 客户服务 (app.customer-service)
**目标用户**：客服主管、客服专员、社群运营

**核心技能**：
| 技能 | 说明 |
|------|------|
| `reply-template` | 回复模板 - 生成常用回复模板 |
| `sentiment-detect` | 情绪识别 - 识别客户情绪和紧急程度 |
| `faq-build` | FAQ 构建 - 整理和更新常见问题 |
| `escalation-decide` | 升级判断 - 判断是否需要升级处理 |
| `satisfaction-analyze` | 满意度分析 - 分析客户满意度趋势 |

---

## 四、教育培训类

### 11. 在线教育 (app.online-education)
**目标用户**：在线讲师、培训师、知识付费创作者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `course-outline` | 课程大纲 - 设计课程结构和目标 |
| `lesson-plan` | 教案设计 - 单节课教学设计 |
| `quiz-gen` | 测验生成 - 生成课后测验题 |
| `assignment-create` | 作业设计 - 设计实践作业 |
| `student-feedback` | 学员反馈 - 分析和整理学员反馈 |
| `certificate-gen` | 证书生成 - 生成结业证书 |

**目录结构**：
```
课程/
├── 大纲/
├── 讲义/
├── 视频/
└── 素材/
作业/
测验/
学员数据/
```

---

### 12. 学习助手 (app.learning-assistant)
**目标用户**：学生、终身学习者、备考人员

**核心技能**：
| 技能 | 说明 |
|------|------|
| `note-organize` | 笔记整理 - 整理和结构化笔记 |
| `flashcard-gen` | 闪卡生成 - 生成记忆卡片 |
| `concept-explain` | 概念解释 - 用简单语言解释复杂概念 |
| `study-plan` | 学习计划 - 制定学习时间表 |
| `practice-problem` | 练习题 - 生成针对性练习题 |
| `knowledge-map` | 知识图谱 - 构建知识关联图 |

---

## 五、研究分析类

### 13. 数据分析 (app.data-analysis)
**目标用户**：数据分析师、产品经理、业务分析师

**核心技能**：
| 技能 | 说明 |
|------|------|
| `data-clean` | 数据清洗 - 处理和清洗原始数据 |
| `sql-gen` | SQL 生成 - 根据需求生成 SQL 查询 |
| `chart-suggest` | 图表建议 - 建议合适的可视化方式 |
| `insight-extract` | 洞察提取 - 从数据中提取关键洞察 |
| `report-gen` | 报告生成 - 生成数据分析报告 |
| `anomaly-detect` | 异常检测 - 识别数据异常 |

**目录结构**：
```
数据源/
├── 原始数据/
└── 清洗后/
分析脚本/
报告/
图表/
```

---

### 14. 市场研究 (app.market-research)
**目标用户**：市场研究员、战略分析师、投资分析师

**核心技能**：
| 技能 | 说明 |
|------|------|
| `industry-scan` | 行业扫描 - 收集行业动态信息 |
| `trend-analyze` | 趋势分析 - 识别市场趋势 |
| `competitor-profile` | 竞品画像 - 构建竞争对手档案 |
| `swot-analyze` | SWOT 分析 - 生成 SWOT 分析 |
| `report-synthesize` | 报告综合 - 整合多源信息成报告 |
| `data-visualize` | 数据可视化 - 生成研究图表 |

---

## 六、设计创意类

### 15. UI/UX 设计 (app.ui-ux-design)
**目标用户**：UI 设计师、产品设计师、交互设计师

**核心技能**：
| 技能 | 说明 |
|------|------|
| `user-flow` | 用户流程 - 设计用户操作流程 |
| `wireframe-gen` | 线框图生成 - 生成页面线框图描述 |
| `copy-ux` | UX 文案 - 编写界面文案 |
| `design-review` | 设计评审 - 提供设计改进建议 |
| `style-guide` | 风格指南 - 生成设计规范文档 |
| `usability-check` | 可用性检查 - 检查交互可用性问题 |

---

### 16. 品牌设计 (app.brand-design)
**目标用户**：品牌设计师、创意总监、品牌经理

**核心技能**：
| 技能 | 说明 |
|------|------|
| `brand-story` | 品牌故事 - 撰写品牌叙事 |
| `naming-brainstorm` | 命名头脑风暴 - 生成品牌/产品名称 |
| `tagline-gen` | 标语生成 - 生成品牌标语 |
| `color-suggest` | 配色建议 - 建议品牌配色方案 |
| `voice-tone` | 语调定义 - 定义品牌语调风格 |
| `brand-audit` | 品牌审计 - 评估品牌一致性 |

---

## 七、办公管理类

### 17. 智能办公 (app.smart-office) ⭐ P0
**目标用户**：职场人士、行政人员、文档工作者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `word-create` | Word 文档创建 - 创建报告、方案、合同等文档 |
| `word-edit` | Word 文档编辑 - 格式调整、内容修改、批注处理 |
| `ppt-create` | PPT 创建 - 根据主题生成演示文稿大纲和内容 |
| `ppt-polish` | PPT 优化 - 优化排版、建议配图、精简文字 |
| `excel-formula` | Excel 公式 - 生成和解释复杂公式 |
| `excel-analysis` | Excel 数据分析 - 数据透视、图表生成、趋势分析 |
| `excel-clean` | Excel 数据清洗 - 去重、格式统一、异常值处理 |
| `pdf-extract` | PDF 提取 - 提取文字、表格、图片 |
| `pdf-convert` | PDF 转换 - PDF 与其他格式互转 |
| `doc-template` | 文档模板 - 生成常用文档模板 |
| `batch-process` | 批量处理 - 批量重命名、格式转换、合并拆分 |

**目录结构**：
```
文档/
├── Word/
├── PPT/
├── Excel/
└── PDF/
模板/
输出/
```

**常见使用场景**：
- 周报/月报自动生成
- 会议 PPT 快速制作
- 销售数据 Excel 分析
- 合同 PDF 信息提取
- 批量文档格式转换

---

### 18. 项目管理 (app.project-management)
**目标用户**：项目经理、Scrum Master、团队负责人

**核心技能**：
| 技能 | 说明 |
|------|------|
| `task-breakdown` | 任务拆解 - 将大任务拆分为子任务 |
| `sprint-plan` | Sprint 规划 - 规划迭代内容 |
| `standup-summary` | 站会总结 - 总结和记录站会内容 |
| `risk-identify` | 风险识别 - 识别项目风险 |
| `progress-report` | 进度报告 - 生成项目进度报告 |
| `retrospective` | 复盘总结 - 引导团队复盘 |

**目录结构**：
```
项目/
├── 文档/
├── 会议记录/
└── 报告/
资源/
归档/
```

---

### 19. 会议管理 (app.meeting)
**目标用户**：会议组织者、行政助理、团队管理者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `agenda-create` | 议程创建 - 制定会议议程 |
| `minutes-gen` | 会议纪要 - 生成会议纪要 |
| `action-extract` | 行动项提取 - 提取待办事项 |
| `follow-up` | 跟进提醒 - 生成跟进邮件 |
| `schedule-optimize` | 时间优化 - 优化会议时间安排 |

---

### 20. 知识管理 (app.knowledge-management)
**目标用户**：知识管理者、文档管理员、团队 Wiki 维护者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `doc-classify` | 文档分类 - 自动分类和打标签 |
| `summary-gen` | 摘要生成 - 生成文档摘要 |
| `link-suggest` | 关联建议 - 建议相关文档链接 |
| `outdated-detect` | 过期检测 - 识别需要更新的文档 |
| `search-optimize` | 搜索优化 - 优化文档可搜索性 |
| `onboard-guide` | 入职指南 - 生成新人入职文档 |

---

## 八、专业领域类

### 21. 法律助手 (app.legal-assistant)
**目标用户**：法务人员、律师助理、合同管理人员

**核心技能**：
| 技能 | 说明 |
|------|------|
| `contract-review` | 合同审查 - 审查合同条款风险 |
| `clause-suggest` | 条款建议 - 建议合同条款 |
| `legal-research` | 法律研究 - 检索相关法律法规 |
| `compliance-check` | 合规检查 - 检查合规性问题 |
| `document-template` | 文书模板 - 生成法律文书模板 |

---

### 22. 财务分析 (app.finance-analysis)
**目标用户**：财务分析师、CFO、投资经理

**核心技能**：
| 技能 | 说明 |
|------|------|
| `financial-model` | 财务建模 - 构建财务预测模型 |
| `report-analyze` | 财报分析 - 分析财务报表 |
| `budget-plan` | 预算规划 - 制定预算方案 |
| `cashflow-forecast` | 现金流预测 - 预测现金流状况 |
| `ratio-calculate` | 比率计算 - 计算和解释财务比率 |

---

### 23. 人力资源 (app.human-resources)
**目标用户**：HR 经理、招聘专员、HRBP

**核心技能**：
| 技能 | 说明 |
|------|------|
| `job-desc-gen` | 职位描述 - 生成职位描述 |
| `resume-screen` | 简历筛选 - 筛选和评估简历 |
| `interview-question` | 面试问题 - 生成面试问题 |
| `offer-letter` | Offer 信 - 生成录用通知 |
| `onboard-checklist` | 入职清单 - 生成入职流程清单 |
| `performance-review` | 绩效评估 - 协助绩效评估 |

**目录结构**：
```
招聘/
├── 职位/
├── 简历/
└── 面试记录/
员工档案/
培训资料/
制度文档/
```

---

## 九、个人效率类

### 24. 个人助理 (app.personal-assistant)
**目标用户**：职场人士、自由职业者、创业者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `email-draft` | 邮件起草 - 起草各类邮件 |
| `schedule-manage` | 日程管理 - 规划每日日程 |
| `task-prioritize` | 任务优先级 - 确定任务优先级 |
| `summary-daily` | 每日总结 - 生成每日工作总结 |
| `reminder-set` | 提醒设置 - 设置重要提醒 |
| `quick-research` | 快速调研 - 快速收集信息 |

---

### 25. 写作助手 (app.writing-assistant)
**目标用户**：作家、编辑、内容工作者

**核心技能**：
| 技能 | 说明 |
|------|------|
| `brainstorm` | 头脑风暴 - 生成创意和想法 |
| `outline-create` | 大纲创建 - 创建文章大纲 |
| `draft-expand` | 草稿扩展 - 扩展简短草稿 |
| `style-edit` | 风格编辑 - 调整写作风格 |
| `grammar-check` | 语法检查 - 检查语法和拼写 |
| `translate-assist` | 翻译辅助 - 翻译和本地化 |

---

## 十、垂直行业类

### 26. 医疗健康 (app.healthcare)
**目标用户**：医疗内容创作者、健康管理师、医学编辑

**核心技能**：
| 技能 | 说明 |
|------|------|
| `health-content` | 健康内容 - 创作健康科普内容 |
| `term-explain` | 术语解释 - 解释医学术语 |
| `reference-check` | 参考核查 - 核查医学引用准确性 |
| `patient-edu` | 患者教育 - 生成患者教育材料 |

---

### 27. 房地产 (app.real-estate)
**目标用户**：房产经纪人、房产营销、物业管理

**核心技能**：
| 技能 | 说明 |
|------|------|
| `listing-write` | 房源描述 - 撰写房源描述 |
| `market-report` | 市场报告 - 生成房产市场报告 |
| `client-match` | 客户匹配 - 分析客户需求匹配房源 |
| `tour-script` | 带看话术 - 生成带看介绍话术 |

---

### 28. 餐饮服务 (app.food-service)
**目标用户**：餐厅老板、外卖运营、美食博主

**核心技能**：
| 技能 | 说明 |
|------|------|
| `menu-design` | 菜单设计 - 设计菜单描述 |
| `recipe-create` | 食谱创作 - 创作和记录食谱 |
| `review-manage` | 评价管理 - 管理和回复评价 |
| `promo-content` | 促销内容 - 创作促销文案 |

---

## 优先级建议

### P0 - 核心场景（立即实现）
1. ✅ 自媒体创作 (app.creator-media) - 已完成
2. ⭐ 智能办公 (app.smart-office) - Word/PPT/Excel/PDF 处理
3. 软件开发 (app.software-dev)
4. 写作助手 (app.writing-assistant)
5. 个人助理 (app.personal-assistant)

### P1 - 高频场景（近期实现）
5. 市场营销 (app.marketing)
6. 数据分析 (app.data-analysis)
7. 项目管理 (app.project-management)
8. 电商运营 (app.ecommerce)

### P2 - 扩展场景（中期实现）
9. 技术写作 (app.tech-writing)
10. 在线教育 (app.online-education)
11. 客户服务 (app.customer-service)
12. 视频脚本 (app.video-script)

### P3 - 专业场景（按需实现）
13. 其他专业领域应用

---

## 下一步行动

1. [ ] 确认优先级和范围
2. [ ] 为 P0 应用详细设计技能规格
3. [ ] 创建应用配置文件 (manifest.json)
4. [ ] 编写应用 AGENTS.md
5. [ ] 实现核心技能 SKILL.md

---

*最后更新: 2026-01-29*
