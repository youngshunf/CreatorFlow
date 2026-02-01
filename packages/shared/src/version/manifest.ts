import { debug } from "../utils/debug";

// Disabled - external version check service removed
// TODO: Set up self-hosted version check service if needed
// const VERSIONS_URL = 'https://agents.craft.do/electron';

export async function getLatestVersion(): Promise<string | null> {
    // Version check disabled - external service removed
    debug('[manifest] Version check disabled - external service removed');
    return null;
}

export async function getManifest(_version: string): Promise<VersionManifest | null> {
    // Manifest fetch disabled - external service removed
    debug('[manifest] Manifest fetch disabled - external service removed');
    return null;
}


export interface BinaryInfo {
  url: string;
  sha256: string;
  size: number;
  filename?: string;
}

export interface VersionManifest {
  version: string;
  build_time: string;
  build_timestamp: number;
  binaries: Record<string, BinaryInfo>;
}
