const daemonPath = 'data/daemon.r';

const { argv } = process;
const env = argv[argv.indexOf('--env') + 1] ?? 'development';

const multiscViewer = {
	name: 'multisc-viewer',
	script: 'vite',
	args: env === 'production' ? 'preview' : 'dev',
	interpreter: 'none',
	watch: false,
	env: {
		NODE_ENV: 'development',
		DATASETS_DIR: 'data/datasets',
		DATASETS_META: 'meta.json',
		PUBLICATIONS_DIR: 'data/publications',
		PUBLICATIONS_META: 'meta.json',
		PLOT_DIR: 'data/plots',
		DAEMON_PORTS: [8001]
	},
	env_production: {
		NODE_ENV: 'production',
		DAEMON_PORTS: [8001]
	}
};

const multiscDaemons = multiscViewer.env.DAEMON_PORTS.map((port) => ({
	name: `multisc-daemon-${port}`,
	script: 'Rscript',
	args: `-e "plumber::pr_run(plumber::pr('${daemonPath}'),port=${port})"`,
	watch: 'daemon.r'
}));

module.exports = {
	apps: [...multiscDaemons, multiscViewer]
};
