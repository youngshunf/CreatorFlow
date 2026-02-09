import { debug } from "../utils/debug";
import { getCloudApiUrl } from "../config/environments";

/**
 * 获取当前平台和架构
 */
function getPlatformInfo(): { platform: string; arch: string } {
    const platform = process.platform; // darwin, win32, linux
    const arch = process.arch; // x64, arm64
    return { platform, arch };
}

/**
 * 获取最新版本号
 */
export async function getLatestVersion(): Promise<string | null> {
    try {
        const { platform, arch } = getPlatformInfo();
        const apiUrl = getCloudApiUrl();
        const url = `${apiUrl}/client/version/latest?platform=${platform}&arch=${arch}&app_code=creator-flow`;
        
        debug(`[manifest] Fetching latest version from: ${url}`);
        
        const response = await fetch(url);
        if (!response.ok) {
            debug(`[manifest] Failed to fetch latest version: ${response.status}`);
            return null;
        }
        
        const result = await response.json() as { code: number; data?: string };
        if (result.code === 200 && result.data) {
            debug(`[manifest] Latest version: ${result.data}`);
            return result.data;
        }
        
        debug('[manifest] No published version available');
        return null;
    } catch (error) {
        debug(`[manifest] Error fetching latest version: ${error instanceof Error ? error.message : String(error)}`);
        return null;
    }
}

/**
 * 获取版本清单（包含所有平台的二进制文件信息）
 */
export async function getManifest(version: string): Promise<VersionManifest | null> {
    try {
        const apiUrl = getCloudApiUrl();
        const url = `${apiUrl}/client/version/manifest/${version}?app_code=creator-flow`;
        
        debug(`[manifest] Fetching manifest for version ${version} from: ${url}`);
        
        const response = await fetch(url);
        if (!response.ok) {
            debug(`[manifest] Failed to fetch manifest: ${response.status}`);
            return null;
        }
        
        const result = await response.json() as { code: number; data?: VersionManifest };
        if (result.code === 200 && result.data) {
            debug(`[manifest] Manifest fetched successfully for version ${version}`);
            return result.data as VersionManifest;
        }
        
        debug('[manifest] No manifest data available');
        return null;
    } catch (error) {
        debug(`[manifest] Error fetching manifest: ${error instanceof Error ? error.message : String(error)}`);
        return null;
    }
}


export interface BinaryInfo {
    url: string;
    sha256: string;
    size: number;
    filename?: string;
}

export interface VersionManifest {
    version: string;
    build_time?: string;
    build_timestamp?: number;
    changelog?: string;
    is_force_update?: boolean;
    min_version?: string;
    binaries: Record<string, BinaryInfo>;
}
