/**
 * Remotion configuration file
 * @see https://www.remotion.dev/docs/config
 */
import { Config } from '@remotion/cli/config';

// Set the video image format for rendering
Config.setVideoImageFormat('jpeg');

// Set the output codec
Config.setCodec('h264');

// Configure the output location
Config.setOverwriteOutput(true);

// Set the number of threads for rendering (0 = auto)
Config.setConcurrency(0);

// Enable caching for faster re-renders
Config.setChromiumDisableWebSecurity(true);

// Configure the browser executable (optional, uses bundled Chromium by default)
// Config.setBrowserExecutable('/path/to/chrome');

// Set the log level
Config.setLogLevel('info');

// Configure the public directory for static assets
Config.setPublicDir('./public');

// Set the entry point for the Remotion project
Config.setEntryPoint('./src/index.ts');
