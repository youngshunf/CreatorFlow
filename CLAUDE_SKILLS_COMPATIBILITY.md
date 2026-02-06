# Claude 技能兼容性

## 目录结构变更

### 修改前
```
~/.creator-flow/global-skills/
```

### 修改后
```
~/.creator-flow/skills/
```

## 兼容性说明

### 1. 与 Claude 技能系统兼容
- 使用相同的目录结构 `~/.creator-flow/skills/`
- 如果用户已经有 Claude 技能，CreatorFlow 会自动识别和加载
- 首次启动时，只复制不存在的技能（不覆盖用户的 Claude 技能）

### 2. 技能加载优先级
```
工作区技能（{workspace}/.creator-flow/skills/）
  ↓ 优先级最高
全局技能（~/.creator-flow/skills/）
  ↓ 包含 CreatorFlow 内置技能 + Claude 技能
内置技能（bundled-skills/）
  ↓ 仅作为复制源
```

### 3. 初始化行为
- 检查 `~/.creator-flow/skills/.meta.json` 是否存在
- 如果不存在，从 `bundled-skills/` 复制技能
- **重要**：如果目标技能目录已存在，跳过复制（保护用户的 Claude 技能）

### 4. 优势
✅ 与 Claude 技能系统无缝集成  
✅ 用户可以在 CreatorFlow 中使用 Claude 技能  
✅ 不会覆盖用户已有的技能  
✅ 统一的技能管理体验

## 测试建议

### 场景 1：全新安装
```bash
# 删除技能目录
rm -rf ~/.creator-flow/skills/

# 启动应用
bun run electron:dev

# 验证
ls ~/.creator-flow/skills/
# 应该看到 33 个内置技能
```

### 场景 2：已有 Claude 技能
```bash
# 模拟已有 Claude 技能
mkdir -p ~/.creator-flow/skills/my-claude-skill
echo "---\nname: My Claude Skill\n---" > ~/.creator-flow/skills/my-claude-skill/SKILL.md

# 启动应用
bun run electron:dev

# 验证
ls ~/.creator-flow/skills/
# 应该看到 my-claude-skill 保持不变，其他内置技能被添加
```

### 场景 3：技能冲突
```bash
# 创建同名技能
mkdir -p ~/.creator-flow/skills/material-organize
echo "---\nname: My Custom Material Organize\n---" > ~/.creator-flow/skills/material-organize/SKILL.md

# 启动应用
bun run electron:dev

# 验证
cat ~/.creator-flow/skills/material-organize/SKILL.md
# 应该保持用户的自定义版本，不被覆盖
```

## 迁移指南

如果用户已经使用了旧版本（`~/.creator-flow/global-skills/`），可以手动迁移：

```bash
# 备份旧目录
mv ~/.creator-flow/global-skills ~/.creator-flow/global-skills.backup

# 重新启动应用，会自动创建新目录
bun run electron:dev

# 如果需要，可以手动复制自定义技能
cp -r ~/.creator-flow/global-skills.backup/my-custom-skill ~/.creator-flow/skills/
```
