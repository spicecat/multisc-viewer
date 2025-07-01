/**
 * Application configuration
 */
export const appConfig = {
  // Server configuration
  server: {
    port: process.env.PORT || 3000,
    host: process.env.HOST || "localhost",
    cors: {
      origin: process.env.CORS_ORIGIN || "*",
      credentials: true,
    },
  },

  // R Daemon configuration
  daemon: {
    maxParallel: 3,
    timeout: 5 * 60 * 1000, // 5 minutes
    scriptPath: "daemon.r",
  },

  // Cache configuration
  cache: {
    ttl: 30 * 60 * 1000, // 30 minutes
    checkPeriod: 60 * 1000, // 1 minute
  },

  // Dataset configuration
  datasets: {
    basePath: "datasets",
    metaFile: "meta.json",
    requiredFiles: ["data.rds", "genes.json"],
  },

  // Client configuration
  client: {
    defaultGroupBy: "Genotype",
    defaultSplitBy: "CellType",
    maxDatasetSelections: 10,
    chartUpdateDebounce: 1000,
  },

  // Theme configuration
  theme: {
    colors: {
      primary: "#2196f3",
      secondary: "#ff4081",
      surface: "#ffffff",
      background: "#f5f5f5",
      error: "#f44336",
    },
  },
} as const;

// Type definitions for configuration
export type AppConfig = typeof appConfig;
export type ServerConfig = typeof appConfig.server;
export type DaemonConfig = typeof appConfig.daemon;
export type CacheConfig = typeof appConfig.cache;
export type DatasetConfig = typeof appConfig.datasets;
export type ClientConfig = typeof appConfig.client;
export type ThemeConfig = typeof appConfig.theme;
