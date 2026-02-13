/**
 * stdio 协议兼容补丁
 *
 * @modelcontextprotocol/sdk 1.x 使用换行分隔 JSON 协议，
 * 而某些客户端（如 Claude Code CLI）使用 Content-Length 头协议。
 *
 * 此补丁让服务端输入端同时支持两种协议：
 * 1. Content-Length 头协议: "Content-Length: N\r\n\r\n{json}"
 * 2. 换行分隔 JSON 协议: "{json}\n"
 *
 * 输出保持原始的换行分隔 JSON 协议（SDK 客户端用 readline 按行读取）。
 */

import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';

/**
 * 应用 stdio 协议兼容补丁
 *
 * 替换 StdioServerTransport 原型上的 ReadBuffer，使输入端支持双协议。
 * 输出端不做修改，保持换行分隔 JSON 协议。
 * 必须在 FastMCP.start() 之前调用。
 */
export function applyStdioCompatPatch(): void {
  // 保存原始的 start 方法
  const originalStart = StdioServerTransport.prototype.start;

  // 替换 start 方法，注入双协议 ReadBuffer
  StdioServerTransport.prototype.start = async function (this: any) {
    // 替换内部的 _readBuffer 为双协议版本
    this._readBuffer = new DualProtocolReadBuffer();

    // 调用原始 start（会绑定 stdin data 事件）
    // 输出端保持原始的换行分隔 JSON 协议（send 方法不修改）
    return originalStart.call(this);
  };
}

/**
 * 双协议 ReadBuffer
 *
 * 自动检测输入是 Content-Length 头协议还是换行分隔 JSON 协议，
 * 并正确解析消息。
 */
class DualProtocolReadBuffer {
  private _buffer: Buffer | undefined;

  append(chunk: Buffer): void {
    this._buffer = this._buffer ? Buffer.concat([this._buffer, chunk]) : chunk;
  }

  readMessage(): any | null {
    if (!this._buffer || this._buffer.length === 0) {
      return null;
    }

    // 检测协议类型：Content-Length 头以 "C" 开头，JSON 以 "{" 开头
    const firstChar = this._buffer[0];

    if (firstChar === 0x43) {
      // 'C' — Content-Length 头协议
      return this._readContentLengthMessage();
    } else if (firstChar === 0x7b) {
      // '{' — 换行分隔 JSON 协议
      return this._readNewlineMessage();
    } else {
      // 跳过空白字符（\r, \n）
      let skip = 0;
      while (skip < this._buffer.length) {
        const b = this._buffer[skip];
        if (b === 0x0a || b === 0x0d || b === 0x20 || b === 0x09) {
          skip++;
        } else {
          break;
        }
      }
      if (skip > 0) {
        this._buffer = this._buffer.subarray(skip);
        if (this._buffer.length === 0) {
          this._buffer = undefined;
          return null;
        }
        return this.readMessage(); // 递归重试
      }
      return null;
    }
  }

  clear(): void {
    this._buffer = undefined;
  }

  /**
   * 解析 Content-Length 头协议消息
   * 格式: "Content-Length: N\r\n\r\n{json}"
   */
  private _readContentLengthMessage(): any | null {
    if (!this._buffer) return null;

    const str = this._buffer.toString('utf8');

    // 查找头部结束标记 \r\n\r\n
    const headerEnd = str.indexOf('\r\n\r\n');
    if (headerEnd === -1) {
      return null; // 头部不完整，等待更多数据
    }

    // 解析 Content-Length
    const header = str.substring(0, headerEnd);
    const match = header.match(/Content-Length:\s*(\d+)/i);
    if (!match) {
      // 无效头部，跳过到下一行
      const nextLine = str.indexOf('\n');
      if (nextLine !== -1) {
        this._buffer = this._buffer.subarray(nextLine + 1);
      } else {
        this._buffer = undefined;
      }
      return null;
    }

    const contentLength = parseInt(match[1]!, 10);
    const headerByteLength = Buffer.byteLength(str.substring(0, headerEnd + 4), 'utf8');

    // 检查 body 是否完整
    if (this._buffer.length < headerByteLength + contentLength) {
      return null; // body 不完整，等待更多数据
    }

    // 提取 body
    const body = this._buffer.subarray(headerByteLength, headerByteLength + contentLength).toString('utf8');
    this._buffer = this._buffer.subarray(headerByteLength + contentLength);
    if (this._buffer.length === 0) {
      this._buffer = undefined;
    }

    try {
      return JSON.parse(body);
    } catch {
      return null;
    }
  }

  /**
   * 解析换行分隔 JSON 协议消息
   * 格式: "{json}\n"
   */
  private _readNewlineMessage(): any | null {
    if (!this._buffer) return null;

    const index = this._buffer.indexOf(0x0a); // '\n'
    if (index === -1) {
      return null; // 消息不完整，等待更多数据
    }

    const line = this._buffer.toString('utf8', 0, index).replace(/\r$/, '');
    this._buffer = this._buffer.subarray(index + 1);
    if (this._buffer.length === 0) {
      this._buffer = undefined;
    }

    try {
      return JSON.parse(line);
    } catch {
      return null;
    }
  }
}
