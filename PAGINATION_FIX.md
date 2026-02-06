# 分页同步修复

## 🐛 发现的问题

**错误日志**:
```
2026-02-06 08:41:05.410 | ERROR | 请求异常: 请求参数非法: size Input should be less than or equal to 200，输入：1000
```

**问题描述**:
- 服务器 API 限制每次请求最多返回 200 条记录
- 原代码请求 1000 条，导致返回 422 错误
- 同步失败，应用市场无法获取数据

## ✅ 修复方案

### 实现分页获取

**修复前**:
```typescript
// 一次性请求 1000 条（超过限制）
const [skillsResponse, appsResponse, categoriesResponse] = await Promise.all([
  listSkills({ page: 1, size: 1000 }),  // ❌ 超过限制
  listApps({ page: 1, size: 1000 }),    // ❌ 超过限制
  listCategories(),
]);
```

**修复后**:
```typescript
// 分页获取所有数据（每页 200 条）
const PAGE_SIZE = 200;

// 获取技能列表（分页）
const allSkills: MarketplaceSkill[] = [];
let skillPage = 1;
let hasMoreSkills = true;

while (hasMoreSkills) {
  const skillsResponse = await listSkills({ page: skillPage, size: PAGE_SIZE });
  allSkills.push(...skillsResponse.items);
  
  hasMoreSkills = skillPage < skillsResponse.pages;
  skillPage++;
  
  debug(`[syncMarketplaceMetadata] 已获取技能: ${allSkills.length}/${skillsResponse.total}`);
}

// 获取应用列表（分页）
const allApps: MarketplaceApp[] = [];
let appPage = 1;
let hasMoreApps = true;

while (hasMoreApps) {
  const appsResponse = await listApps({ page: appPage, size: PAGE_SIZE });
  allApps.push(...appsResponse.items);
  
  hasMoreApps = appPage < appsResponse.pages;
  appPage++;
  
  debug(`[syncMarketplaceMetadata] 已获取应用: ${allApps.length}/${appsResponse.total}`);
}
```

## 📊 影响分析

### 性能影响
- **请求次数**: 从 2 次增加到 N 次（取决于数据量）
- **单次请求时间**: 减少（每次只获取 200 条）
- **总同步时间**: 略有增加，但在可接受范围内
- **内存占用**: 相同（最终都是加载所有数据）

### 数据量估算
假设有 500 个技能和 100 个应用：
- 技能请求次数: 500 / 200 = 3 次
- 应用请求次数: 100 / 200 = 1 次
- 总请求次数: 4 次（vs 原来的 2 次）

### 优点
- ✅ 符合 API 限制，不会返回 422 错误
- ✅ 可以处理任意数量的数据
- ✅ 提供进度日志，便于调试
- ✅ 失败时可以部分重试（如果需要）

### 缺点
- ⚠️ 同步时间略有增加
- ⚠️ 请求次数增加（但每次更快）

## 🧪 测试建议

### 1. 测试正常同步
```bash
# 清理缓存
rm -rf ~/.creator-flow/marketplace/cache/

# 启动应用，观察日志
bun run electron:dev

# 应该看到类似日志：
# [syncMarketplaceMetadata] 已获取技能: 200/500
# [syncMarketplaceMetadata] 已获取技能: 400/500
# [syncMarketplaceMetadata] 已获取技能: 500/500
# [syncMarketplaceMetadata] 同步完成: 500 个技能, 100 个应用, 10 个分类
```

### 2. 验证缓存内容
```bash
# 检查缓存文件
cat ~/.creator-flow/marketplace/cache/meta.json | jq '.skills | length'
cat ~/.creator-flow/marketplace/cache/meta.json | jq '.apps | length'

# 应该显示完整的数据量
```

### 3. 测试大数据量
如果服务器有超过 200 条数据，验证分页逻辑是否正确工作。

## 📝 修改的文件

- `packages/shared/src/marketplace/sync.ts`
  - `syncMarketplaceMetadata()` - 添加分页逻辑
  - `forceSyncMarketplaceMetadata()` - 添加分页逻辑
  - 添加类型导入 `MarketplaceSkill`, `MarketplaceApp`

## 🎯 修复状态

- ✅ 实现分页获取技能列表
- ✅ 实现分页获取应用列表
- ✅ 添加进度日志
- ✅ 更新强制同步函数
- ✅ 添加类型导入

**准备状态**: ✅ 可以测试
