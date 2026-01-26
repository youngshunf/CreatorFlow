import request from './request';

export interface LoginResult {
  access_token: string;
  access_token_expire_time: string;
  refresh_token: string;
  refresh_token_expire_time: string;
  llm_token: string;
  is_new_user: boolean;
  user: {
    uuid: string;
    username: string;
    nickname: string;
    phone: string;
    email: string | null;
    avatar: string | null;
  };
}

export interface SendCodeResponse {
  success: boolean;
  message: string;
}

export const authApi = {
  /**
   * 发送手机验证码
   * @param phone 手机号
   */
  sendSmsCode: (phone: string) => {
    return request.post<SendCodeResponse>('/auth/send-code', { phone });
  },

  /**
   * 手机号验证码登录
   * @param phone 手机号
   * @param code 验证码
   */
  loginByPhone: (phone: string, code: string) => {
    return request.post<LoginResult>('/auth/phone-login', { phone, code });
  },

  /**
   * 账号密码登录 (可选)
   * @param username 用户名
   * @param password 密码
   */
  loginByPassword: (username: string, password: string) => {
    // 假设密码登录接口为 /auth/login
    // 注意：clound-backend 的密码登录参数通常是 form-data (OAuth2PasswordRequestForm) 
    // 或者 json。这里假设是 JSON，如果报错需要调整。
    // 经查，clound-backend/backend/app/admin/api/v1/auth/auth.py 的 login 接口接收 LoginParam
    return request.post<LoginResult>('/auth/login', { username, password });
  }
};
