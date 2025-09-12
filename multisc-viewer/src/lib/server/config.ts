import { env } from '$env/dynamic/private';

export const datasetsConfig = {
	dir: env.DATASETS_DIR ?? '../data/datasets',
	meta: env.DATASETS_META ?? 'meta.json',
	requiredFiles: ['data.rds', 'genes.json', 'cluster.colors.rds', 'genotype.colors.rds']
};

export const publicationsConfig = {
	dir: env.PUBLICATIONS_DIR ?? '../data/publications',
	meta: env.PUBLICATIONS_META ?? 'meta.json'
};

export const plotConfig = {
	cacheDir: env.PLOT_DIR ?? '../data/plots',
	plotName: env.PLOT_NAME ?? 'plot.png',
	stdTTL: 4 * 60 * 60, // 4 hour
	maxWaitTime: 4 * 60 * 1000 // 4 minutes
};

export const daemonConfig = {
	daemons: env.DAEMONS?.split(',') ?? ['http://127.0.0.1:8000', 'http://127.0.0.1:8001'],
	connectionTimeout: 5 * 1000, // 5 seconds
	timeout: 60 * 1000, // 1 minute
	stdTTL: 4 * 60 * 60 // 4 hour
};
