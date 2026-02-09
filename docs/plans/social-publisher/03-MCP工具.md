# MCP Server 工具设计

Agent 通过 MCP 协议调用社交发布功能。

## 工具列表

### 账号管理

#### social_list_accounts

列出所有已登录的社交平台账号。

```typescript
{
  name: "social_list_accounts",
  description: "列出所有已登录的社交平台账号",
  parameters: {
    type: "object",
    properties: {
      platform: {
        type: "string",
        enum: ["xiaohongshu", "douyin", "wechat_mp", "bilibili", "kuaishou", "weibo", "zhihu"],
        description: "可选，筛选指定平台的账号"
      }
    }
  },
  returns: {
    type: "array",
    items: {
      type: "object",
      properties: {
        id: { type: "string" },
        platform: { type: "string" },
        nickname: { type: "string" },
        avatar: { type: "string" },
        status: { type: "string", enum: ["active", "expired", "banned"] },
        fansCount: { type: "number" }
      }
    }
  }
}
```

#### social_login

登录社交平台账号，会打开浏览器窗口让用户完成登录。

```typescript
{
  name: "social_login",
  description: "登录社交平台账号，会打开浏览器窗口让用户扫码或输入验证码",
  parameters: {
    type: "object",
    properties: {
      platform: {
        type: "string",
        enum: ["xiaohongshu", "douyin", "wechat_mp", "bilibili", "kuaishou", "weibo", "zhihu"],
        description: "要登录的平台"
      }
    },
    required: ["platform"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      account: { type: "object" },
      error: { type: "string" }
    }
  }
}
```

#### social_check_login

检查账号登录状态是否有效。

```typescript
{
  name: "social_check_login",
  description: "检查账号登录状态是否有效",
  parameters: {
    type: "object",
    properties: {
      accountId: {
        type: "string",
        description: "账号 ID"
      }
    },
    required: ["accountId"]
  },
  returns: {
    type: "object",
    properties: {
      valid: { type: "boolean" },
      message: { type: "string" }
    }
  }
}
```

#### social_logout

退出登录并删除账号。

```typescript
{
  name: "social_logout",
  description: "退出登录并删除账号的登录态",
  parameters: {
    type: "object",
    properties: {
      accountId: {
        type: "string",
        description: "账号 ID"
      }
    },
    required: ["accountId"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" }
    }
  }
}
```

### 内容发布

#### social_publish

发布内容到指定平台。

```typescript
{
  name: "social_publish",
  description: "发布内容到指定社交平台",
  parameters: {
    type: "object",
    properties: {
      accountId: {
        type: "string",
        description: "发布账号 ID"
      },
      type: {
        type: "string",
        enum: ["image_text", "video", "article", "text"],
        description: "内容类型"
      },
      title: {
        type: "string",
        description: "标题（部分平台必填）"
      },
      content: {
        type: "string",
        description: "正文/描述"
      },
      mediaFiles: {
        type: "array",
        items: { type: "string" },
        description: "媒体文件的本地路径列表"
      },
      tags: {
        type: "array",
        items: { type: "string" },
        description: "话题标签"
      },
      coverImage: {
        type: "string",
        description: "封面图路径（视频内容可选）"
      },
      scheduledAt: {
        type: "string",
        description: "定时发布时间（ISO 8601 格式），不填则立即发布"
      }
    },
    required: ["accountId", "type", "content"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      postId: { type: "string" },
      postUrl: { type: "string" },
      error: { type: "string" },
      needReview: { type: "boolean" }
    }
  }
}
```

#### social_publish_batch

将同一内容发布到多个平台。

```typescript
{
  name: "social_publish_batch",
  description: "将同一内容发布到多个平台（一键多发）",
  parameters: {
    type: "object",
    properties: {
      accountIds: {
        type: "array",
        items: { type: "string" },
        description: "目标账号 ID 列表"
      },
      type: {
        type: "string",
        enum: ["image_text", "video", "article", "text"],
        description: "内容类型"
      },
      title: { type: "string" },
      content: { type: "string" },
      mediaFiles: {
        type: "array",
        items: { type: "string" }
      },
      tags: {
        type: "array",
        items: { type: "string" }
      },
      coverImage: { type: "string" },
      scheduledAt: { type: "string" }
    },
    required: ["accountIds", "type", "content"]
  },
  returns: {
    type: "array",
    items: {
      type: "object",
      properties: {
        accountId: { type: "string" },
        platform: { type: "string" },
        success: { type: "boolean" },
        postUrl: { type: "string" },
        error: { type: "string" }
      }
    }
  }
}
```

### 内容管理

#### social_get_posts

获取已发布内容列表。

```typescript
{
  name: "social_get_posts",
  description: "获取指定账号的已发布内容列表",
  parameters: {
    type: "object",
    properties: {
      accountId: {
        type: "string",
        description: "账号 ID"
      },
      limit: {
        type: "number",
        description: "返回数量限制，默认 20"
      },
      offset: {
        type: "number",
        description: "偏移量，用于分页"
      }
    },
    required: ["accountId"]
  },
  returns: {
    type: "array",
    items: {
      type: "object",
      properties: {
        postId: { type: "string" },
        title: { type: "string" },
        type: { type: "string" },
        coverUrl: { type: "string" },
        postUrl: { type: "string" },
        publishedAt: { type: "string" },
        stats: {
          type: "object",
          properties: {
            views: { type: "number" },
            likes: { type: "number" },
            comments: { type: "number" }
          }
        }
      }
    }
  }
}
```

#### social_get_post_detail

获取已发布内容的详细信息。

```typescript
{
  name: "social_get_post_detail",
  description: "获取已发布内容的详细信息",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      postId: { type: "string" }
    },
    required: ["accountId", "postId"]
  },
  returns: {
    type: "object",
    properties: {
      postId: { type: "string" },
      title: { type: "string" },
      content: { type: "string" },
      type: { type: "string" },
      mediaUrls: { type: "array", items: { type: "string" } },
      tags: { type: "array", items: { type: "string" } },
      postUrl: { type: "string" },
      publishedAt: { type: "string" },
      stats: { type: "object" }
    }
  }
}
```

#### social_delete_post

删除已发布的内容。

```typescript
{
  name: "social_delete_post",
  description: "删除已发布的内容",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      postId: { type: "string" }
    },
    required: ["accountId", "postId"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      error: { type: "string" }
    }
  }
}
```

### 数据统计

#### social_sync_stats

同步指定内容的数据统计。

```typescript
{
  name: "social_sync_stats",
  description: "同步指定内容的数据统计（阅读量、点赞、评论等）",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      postId: { type: "string" }
    },
    required: ["accountId", "postId"]
  },
  returns: {
    type: "object",
    properties: {
      views: { type: "number" },
      likes: { type: "number" },
      comments: { type: "number" },
      shares: { type: "number" },
      favorites: { type: "number" },
      updatedAt: { type: "string" }
    }
  }
}
```

#### social_sync_all_stats

批量同步账号下所有内容的数据统计。

```typescript
{
  name: "social_sync_all_stats",
  description: "同步账号下所有内容的数据统计",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" }
    },
    required: ["accountId"]
  },
  returns: {
    type: "object",
    properties: {
      updated: { type: "number" },
      failed: { type: "number" },
      errors: { type: "array", items: { type: "string" } }
    }
  }
}
```

#### social_get_account_stats

获取账号整体数据统计。

```typescript
{
  name: "social_get_account_stats",
  description: "获取账号整体数据（粉丝数、总阅读量等）",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" }
    },
    required: ["accountId"]
  },
  returns: {
    type: "object",
    properties: {
      fansCount: { type: "number" },
      followingCount: { type: "number" },
      totalPosts: { type: "number" },
      totalViews: { type: "number" },
      totalLikes: { type: "number" },
      updatedAt: { type: "string" }
    }
  }
}
```

### 评论管理

#### social_get_comments

获取内容的评论列表。

```typescript
{
  name: "social_get_comments",
  description: "获取指定内容的评论列表",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      postId: { type: "string" },
      limit: { type: "number", description: "返回数量限制" }
    },
    required: ["accountId", "postId"]
  },
  returns: {
    type: "array",
    items: {
      type: "object",
      properties: {
        commentId: { type: "string" },
        authorName: { type: "string" },
        authorAvatar: { type: "string" },
        content: { type: "string" },
        likes: { type: "number" },
        createdAt: { type: "string" },
        replies: { type: "array" }
      }
    }
  }
}
```

#### social_reply_comment

回复评论。

```typescript
{
  name: "social_reply_comment",
  description: "回复指定评论",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      commentId: { type: "string" },
      text: { type: "string", description: "回复内容" }
    },
    required: ["accountId", "commentId", "text"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      replyId: { type: "string" },
      error: { type: "string" }
    }
  }
}
```

#### social_delete_comment

删除评论。

```typescript
{
  name: "social_delete_comment",
  description: "删除指定评论",
  parameters: {
    type: "object",
    properties: {
      accountId: { type: "string" },
      commentId: { type: "string" }
    },
    required: ["accountId", "commentId"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      error: { type: "string" }
    }
  }
}
```

### 定时任务

#### social_get_scheduled

获取待发布的定时任务列表。

```typescript
{
  name: "social_get_scheduled",
  description: "获取待发布的定时任务列表",
  parameters: {
    type: "object",
    properties: {
      accountId: {
        type: "string",
        description: "可选，筛选指定账号的任务"
      }
    }
  },
  returns: {
    type: "array",
    items: {
      type: "object",
      properties: {
        taskId: { type: "string" },
        contentId: { type: "string" },
        accountIds: { type: "array", items: { type: "string" } },
        scheduledAt: { type: "string" },
        status: { type: "string" }
      }
    }
  }
}
```

#### social_cancel_scheduled

取消定时发布任务。

```typescript
{
  name: "social_cancel_scheduled",
  description: "取消定时发布任务",
  parameters: {
    type: "object",
    properties: {
      taskId: { type: "string" }
    },
    required: ["taskId"]
  },
  returns: {
    type: "object",
    properties: {
      success: { type: "boolean" },
      error: { type: "string" }
    }
  }
}
```

## Agent 调用示例

### 示例 1：发布图文到小红书

```
用户: 帮我把这篇文章发到小红书，图片在 ~/Desktop/photos 目录下

Agent 思考:
1. 先获取小红书账号
2. 准备内容并发布

Agent 调用:
→ social_list_accounts({ platform: "xiaohongshu" })
← [{ id: "xhs_001", nickname: "创作者小明", status: "active" }]

→ social_publish({
    accountId: "xhs_001",
    type: "image_text",
    title: "今日分享",
    content: "这是文章内容...",
    mediaFiles: [
      "/Users/xxx/Desktop/photos/1.jpg",
      "/Users/xxx/Desktop/photos/2.jpg"
    ],
    tags: ["生活", "分享"]
  })
← { success: true, postUrl: "https://www.xiaohongshu.com/..." }

Agent 回复: 已成功发布到小红书，链接：https://www.xiaohongshu.com/...
```

### 示例 2：一键多平台发布

```
用户: 把这个视频同时发到抖音、B站和快手

Agent 调用:
→ social_list_accounts({})
← [
    { id: "douyin_001", platform: "douyin", nickname: "小明说事" },
    { id: "bili_001", platform: "bilibili", nickname: "小明UP" },
    { id: "ks_001", platform: "kuaishou", nickname: "小明快手" }
  ]

→ social_publish_batch({
    accountIds: ["douyin_001", "bili_001", "ks_001"],
    type: "video",
    title: "今日vlog",
    content: "记录美好生活",
    mediaFiles: ["/path/to/video.mp4"],
    tags: ["vlog", "日常"]
  })
← [
    { accountId: "douyin_001", platform: "douyin", success: true, postUrl: "..." },
    { accountId: "bili_001", platform: "bilibili", success: true, postUrl: "..." },
    { accountId: "ks_001", platform: "kuaishou", success: true, postUrl: "..." }
  ]
```

### 示例 3：查看数据统计

```
用户: 看看我小红书最近发布的内容数据怎么样

Agent 调用:
→ social_list_accounts({ platform: "xiaohongshu" })
← [{ id: "xhs_001", nickname: "创作者小明" }]

→ social_get_posts({ accountId: "xhs_001", limit: 5 })
← [
    { postId: "p1", title: "穿搭分享", stats: { views: 2300, likes: 156 } },
    { postId: "p2", title: "美食探店", stats: { views: 5600, likes: 423 } },
    ...
  ]

Agent 回复:
您最近 5 篇小红书笔记的数据：
1. 《穿搭分享》- 阅读 2.3k，点赞 156
2. 《美食探店》- 阅读 5.6k，点赞 423
...
```
