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
	server: 'http://localhost',
	ports: env.DAEMON_PORTS ? (JSON.parse(env.DAEMON_PORTS) as number[]) : [8000, 8001],
	stdTTL: 4 * 60 * 60 // 4 hour
};
