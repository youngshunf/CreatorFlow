import request from './request';

export interface UserProfile {
  uuid: string;
  username: string;
  nickname: string;
  avatar: string | null;
  email: string | null;
  phone: string | null;
  gender: string | null;  // male/female/other
  birthday: string | null;  // YYYY-MM-DD
  province: string | null;
  city: string | null;
  district: string | null;
  industry: string | null;
  bio: string | null;
  status: number;
  is_superuser: boolean;
  is_staff: boolean;
  is_multi_login: boolean;
  join_time: string;
  last_login_time: string | null;
}

export interface UpdateUserProfileParams {
  nickname?: string;
  avatar?: string;
  gender?: string;
  birthday?: string;
  province?: string;
  city?: string;
  district?: string;
  industry?: string;
  bio?: string;
}

export interface ChangePasswordParams {
  old_password: string;
  new_password: string;
  confirm_password: string;
}

export interface ChangePasswordBySmsParams {
  code: string;
  new_password: string;
  confirm_password: string;
}

export interface UploadAvatarResponse {
  url: string;
}

export const userApi = {
  /**
   * 获取当前用户信息
   */
  getCurrentUser: () => {
    return request.get<UserProfile>('/sys/users/me');
  },

  /**
   * 更新当前用户资料
   * @param data 用户资料数据
   */
  updateProfile: (data: UpdateUserProfileParams) => {
    return request.put('/sys/users/me/profile', data);
  },

  /**
   * 上传用户头像
   * @param file 头像文件
   */
  uploadAvatar: (file: File) => {
    const formData = new FormData();
    formData.append('file', file);
    return request.upload<UploadAvatarResponse>('/sys/users/me/avatar/upload', formData);
  },

  /**
   * 修改密码（旧密码验证）
   * @param data 密码修改参数
   */
  changePassword: (data: ChangePasswordParams) => {
    return request.put('/sys/users/me/password', data);
  },

  /**
   * 通过短信验证码修改密码
   * @param data 密码修改参数
   */
  changePasswordBySms: (data: ChangePasswordBySmsParams) => {
    return request.put('/sys/users/me/password/sms', data);
  },

  /**
   * 发送短信验证码
   * @param phone 手机号
   */
  sendSmsCode: (phone: string) => {
    return request.post('/auth/send-code', { phone });
  },
};
