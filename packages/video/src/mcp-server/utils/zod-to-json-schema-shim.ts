/**
 * zod-to-json-schema shim for Zod v4
 *
 * zod-to-json-schema v3 不支持 Zod v4 的内部结构。
 * 这个 shim 用 Zod v4 自带的 .toJSONSchema() 方法替代。
 * 在 esbuild 打包时通过 alias 注入。
 */

export function zodToJsonSchema(schema: any): Record<string, unknown> {
  if (typeof schema?.toJSONSchema === 'function') {
    // Zod v4: 使用内置的 toJSONSchema()
    const jsonSchema = schema.toJSONSchema();
    // 移除 $schema 字段，保持与 zod-to-json-schema v3 输出一致
    delete jsonSchema.$schema;
    return jsonSchema;
  }

  // Fallback: 返回空 object schema
  return { type: 'object', properties: {}, additionalProperties: false };
}

export default zodToJsonSchema;
