// Centralized, type-safe app configuration
import { env } from '$env/static/private';

export const CONFIG = {
  publications: {
    basePath: env.PUBLICATIONS_BASE_PATH ?? './publications',
    metaFile: env.PUBLICATIONS_META_FILE ?? 'meta.json'
  },
  datasets: {
    basePath: env.DATASETS_BASE_PATH ?? './datasets',
    metaFile: env.DATASETS_META_FILE ?? 'meta.json',
    requiredFiles: (env.DATASETS_REQUIRED_FILES ?? 'data.rds,genes.json').split(',')
  },
  daemon: {
    server: env.DAEMON_SERVER ?? 'http://localhost',
    ttl: Number(env.DAEMON_TTL ?? 3600),
    ports: (env.DAEMON_PORTS ?? '5001,5002').split(',').map(Number)
  },
  plot: {
    ttl: Number(env.PLOT_TTL ?? 3600)
  }
} as const;
