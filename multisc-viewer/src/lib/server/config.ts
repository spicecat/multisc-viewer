import { range } from "lodash-es";
import { env } from "$env/dynamic/private";

export const plotsConfig = {
	/** Directory to save plots */
	plotsDir: env.PLOTS_DIR ?? "../data/plots",
	/** Timeout for plot requests in milliseconds */
	timeout: Number(env.PLOT_TIMEOUT ?? 1000 * 60 * 4), // 4 minutes
	/** Maximum size of the plot directory in bytes */
	maxSize: Number(env.PLOT_MAX_SIZE ?? 1024 * 1024 * 1024 * 4), // 4 GB
};

export const daemonConfig = {
	/** URLs of daemon servers */
	daemonUrls:
		env.DAEMONS?.split(",") ??
		range(8080, 8090).map((p) => `http://localhost:${p}`),
	/** Timeout for connection to daemon servers in milliseconds */
	connectionTimeout: Number(env.CONNECTION_TIMEOUT ?? 1000 * 5), // 5 seconds
	/** Timeout for daemon requests in milliseconds */
	timeout: plotsConfig.timeout, // 4 minutes
	/** Maximum size for daemon datasets cache */
	maxSize: Number(env.DAEMON_MAX_SIZE ?? 1024 * 1024 * 1024 * 4), // 4 GB,
	/** Time-to-live for daemon datasets cache in milliseconds */
	ttl: Number(env.DAEMON_TTL ?? 1000 * 60 * 10), // 10 minutes
	/** Interval to refresh daemon registration in milliseconds */
	refresh: Number(env.DAEMON_REFRESH ?? 1000 * 60 * 10), // 10 minutes
};
