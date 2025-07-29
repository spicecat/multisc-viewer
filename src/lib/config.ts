import { env } from '$env/dynamic/private';

export const datasetsConfig = {
	dir: env.DATASETS_DIR ?? 'data/datasets',
	meta: env.DATASETS_META ?? 'meta.json',
	requiredFiles: ['data.rds', 'genes.json', 'cluster.colors.rds', 'genotype.colors.rds']
};

export const publicationsConfig = {
	dir: env.PUBLICATIONS_DIR ?? 'data/publications',
	meta: env.PUBLICATIONS_META ?? 'meta.json'
};

export const daemonConfig = {
	server: 'http://localhost',
	ports: env.DAEMON_PORTS ? env.DAEMON_PORTS.split(',').map(Number) : [3000],
	stdTTL: env.DATASET_TTL ? parseInt(env.DATASET_TTL, 10) : 60 * 60 // Default 1 hour
};

export const plotConfig = {
	cacheDir: env.PLOT_DIR ?? 'data/plots',
	stdTTL: env.PLOT_TTL ? parseInt(env.PLOT_TTL, 10) : 60 * 60 // Default 1 hour
};
