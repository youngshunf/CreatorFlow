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

export const userApi = {
  /**
   * 获取当前用户信息
   */
  getCurrentUser: () => {
    return request.get<UserProfile>('/admin/v1/users/me');
  },

  /**
   * 更新当前用户资料
   * @param data 用户资料数据
   */
  updateProfile: (data: UpdateUserProfileParams) => {
    return request.put('/admin/v1/users/me/profile', data);
  },
};
