import type { paths } from '$lib/types/api';
import type { Datasets, DEGs, PlotsParams } from '$lib/types/daemon';
import { LRUCache } from 'lru-cache';
import createClient, { type Client, type Middleware } from 'openapi-fetch';
import { daemonConfig } from './config';

const { daemonUrls, connectionTimeout, timeout, maxSize, ttl, refresh } = daemonConfig;
const defaultSize = 1024 * 1024; // default size for datasets with unknown size 1 MB

export const daemons: Daemon[] = [];

/**
 * Class representing a connection to a daemon server.
 */
class Daemon {
	private client: Client<paths>;
	private baseUrl: string;
	private cache: LRUCache<string, string>;
	requestLoad = 0;

	/**
	 * Create a Daemon instance.
	 * @param baseUrl - URL of the daemon server
	 */
	constructor(baseUrl: string) {
		this.baseUrl = baseUrl;
		const logger: Middleware = {
			async onRequest({ request }) {
				console.log(`Request ${request.method} ${request.url}`);
			},
			async onResponse({ response }) {
				console.log(`Response ${response.status} ${response.url}`);
			},
			async onError({ request, error }) {
				console.error(`Error ${request.url} - ${(<Error>error).message}`);
			}
		};

		const requestCache: LRUCache<string, Response> = new LRUCache({ ttl, ttlAutopurge: true });
		const cacheMiddleware: Middleware = {
			async onRequest({ request, schemaPath }) {
				if (request.method === 'GET' && schemaPath !== '/health' && requestCache.has(request.url)) {
					console.log(`Cache hit for ${request.url}`);
					return Response.json(requestCache.get(request.url));
				}
			},
			async onResponse({ request, response, schemaPath }) {
				if (request.method === 'GET' && schemaPath !== '/health')
					requestCache.set(request.url, await response.clone().json());
			}
		};

		this.client = createClient<paths>({ baseUrl });
		this.client.use(logger);
		this.client.use(cacheMiddleware);

		this.cache = new LRUCache<string, string>({
			maxSize,
			ttl,
			dispose: this.unload
		});
	}

	/** Priority for loading datasets.
	 * Lower priority value means higher priority.
	 * Priority is based on current cache size and request load.
	 * @param ds - dataset id
	 * @returns Priority value
	 */
	priority = (ds: string) =>
		this.cache.calculatedSize + this.requestLoad - (this.cache.has(ds) ? this.cache.maxSize : 0);

	/**
	 * Check the health of the daemon server.
	 * @returns Health status of the daemon
	 */
	health = async () => {
		try {
			const { data, error } = await this.client.GET('/health', {
				signal: AbortSignal.timeout(connectionTimeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return { status: 'unhealthy' };
		}
	};

	/**
	 * Get the list of datasets currently loaded in the daemon.
	 * @returns List of all loaded dataset ids
	 */
	loaded = async () => {
		try {
			const { data, error } = await this.client.GET('/loaded', {
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			const datasets = await this.datasets(data);
			for (const ds of data) this.cache.set(ds, ds, { size: datasets[ds].size ?? defaultSize });
			return data;
		} catch {
			return [];
		}
	};

	/**
	 * Unload a dataset from the daemon.
	 * @param ds - dataset ids to unload
	 * @returns List of all loaded dataset ids
	 */
	private unload = async (ds: string) => {
		try {
			const { data, error } = await this.client.POST('/unload', {
				body: { ds: [ds] },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			this.cache.delete(ds);
			return data;
		} catch {
			return [];
		}
	};

	/**
	 * Load datasets into the daemon.
	 * @param ds - list of dataset ids to load
	 * @returns List of all loaded dataset ids
	 */
	load = async (ds: string[]) => {
		try {
			const { data, error } = await this.client.POST('/load', {
				body: { ds },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			const datasets = await this.datasets(data);
			for (const ds of data) this.cache.set(ds, ds, { size: datasets[ds].size ?? defaultSize });
			return data;
		} catch {
			return [];
		}
	};

	/**
	 * Fetch datasets metadata.
	 * @param ds - list of dataset ids to fetch; if `undefined`, fetch all datasets
	 * @returns Dictionary mapping dataset ids to metadata
	 */
	datasets = async (ds?: string[]) => {
		try {
			const { data, error } = await this.client.GET('/datasets', {
				params: { query: { ds } },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return {};
		}
	};

	/**
	 * Fetch publications metadata.
	 * @param pub - list of publication ids to fetch; if `undefined`, fetch all publications
	 * @returns Dictionary mapping publication ids to metadata
	 */
	publications = async (pub?: string[]) => {
		try {
			const { data, error } = await this.client.GET('/publications', {
				params: { query: { pub } },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return {};
		}
	};

	/**
	 * Fetch genes for datasets.
	 * @param ds - list of dataset ids to fetch
	 * @returns Dictionary mapping dataset ids to genes
	 */
	genes = async (ds: string[]) => {
		try {
			const { data, error } = await this.client.GET('/genes', {
				params: { query: { ds } },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return {};
		}
	};

	/**
	 * Fetch differentially expressed genes for datasets.
	 * @param ds - list of dataset ids to fetch
	 * @returns Dictionary mapping dataset ids to dictionary mapping deg id to differentially expressed genes
	 */
	degs = async (ds: string[]): Promise<DEGs> => {
		try {
			const { data, error } = await this.client.GET('/degs', {
				params: { query: { ds } },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return {};
		}
	};

	/**
	 * Fetch rows for gene datasets and differentially expressed genes.
	 * @param ds - list of dataset ids to fetch
	 * @returns List of gene rows
	 * */
	genesRows = async (ds: string[]) => {
		try {
			const { data, error } = await this.client.GET('/genes-rows', {
				params: { query: { ds } },
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return [];
		}
	};

	/** Render plots.
	 * @param plotsParams - parameters for plot rendering
	 * @returns List of rendered plots paths
	 */
	plots = async (plotsParams: PlotsParams) => {
		try {
			this.requestLoad += maxSize;
			const datasets = await this.datasets(plotsParams.ds);
			const { data, error } = await this.client.POST('/plots', {
				body: {
					...plotsParams,
					ds: plotsParams.ds.sort(
						(a, b) => (datasets[a].size ?? defaultSize) - (datasets[b].size ?? defaultSize)
					)
				},
				signal: AbortSignal.timeout(timeout)
			});
			if (error) throw error;
			return data;
		} catch {
			return [];
		} finally {
			this.requestLoad -= maxSize;
		}
	};
}

const registerDaemons = async () => {
	daemons.length = 0;
	await Promise.all(
		daemonUrls.map(async (daemonUrl) => {
			const daemon = new Daemon(daemonUrl);
			const { status } = await daemon.health();
			if (status !== 'unhealthy') daemons.push(daemon);
		})
	);
};

setInterval(registerDaemons, refresh);
await registerDaemons();

/**
 * Map daemons to target dataset ids.
 * @param ds - list of dataset ids; if empty, target all datasets
 * @returns Map of Daemon to list of target dataset ids
 */
export const getDaemonTargets = async (ds?: string[], load: boolean = false) => {
	const daemonTargets = new Map<Daemon, string[]>();
	if (ds?.length === 0) return daemonTargets;

	const datasetDaemons = new Map<string, Daemon[]>();
	const daemonDatasets = new Map<Daemon, Datasets>(
		await Promise.all(
			daemons.map(async (daemon): Promise<[Daemon, Datasets]> => {
				const datasets = await daemon.datasets(ds);
				for (const d of Object.keys(datasets)) {
					if (!datasetDaemons.has(d)) datasetDaemons.set(d, []);
					datasetDaemons.get(d)!.push(daemon);
				}
				return [daemon, datasets];
			})
		)
	);

	for (const [d] of datasetDaemons) {
		const daemon = datasetDaemons
			.get(d)!
			.reduce((minDaemon, currDaemon) =>
				currDaemon.priority(d) < minDaemon.priority(d) ? currDaemon : minDaemon
			);
		if (!daemonTargets.has(daemon)) daemonTargets.set(daemon, []);
		daemonTargets.get(daemon)!.push(d);
		daemon.requestLoad += (load ? 1 : -1) * (daemonDatasets.get(daemon)![d].size ?? defaultSize);
	}
	for (const [daemon, ds] of daemonTargets)
		for (const d of ds)
			daemon.requestLoad -= (load ? 1 : -1) * (daemonDatasets.get(daemon)![d].size ?? defaultSize);

	return daemonTargets;
};
