const multiscViewer = {
	name: 'multisc-viewer',
	script: 'pnpm',
	args: process.env.NODE_ENV === 'production' ? 'preview' : 'dev',
	interpreter: 'none',
	env: {
		DAEMON_PORTS: [8001, 8002, 8003]
	},
	env_production: {
		DAEMON_PORTS: [
			8001, 8002, 8003, 8004, 8005, 8006, 8007, 8008, 8009, 8010, 8011, 8012, 8013, 8014, 8015,
			8016, 8017, 8018, 8019, 8020
		]
	}
};

const multiscDaemons = multiscViewer.env.DAEMON_PORTS.map((port) => ({
	name: `multisc-daemon-${port}`,
	script: 'Rscript',
	args: `-e "plumber::pr_run(plumber::pr('daemon.r'), port = ${port})"`,
	watch: 'daemon.r'
}));

module.exports = {
	apps: [multiscViewer]
};
