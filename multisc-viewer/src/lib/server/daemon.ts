import type {
	Datasets,
	DEGs,
	Genes,
	Health,
	Loaded,
	Load,
	Plots,
	PlotsParams,
	Publications,
	Unload
} from '$lib/types/daemon';
import axios from 'axios';
import { setupCache, type AxiosCacheInstance } from 'axios-cache-interceptor';
import { LRUCache } from 'lru-cache';
import { daemonConfig } from './config';

const { daemonUrls, connectionTimeout, timeout, maxSize, ttl, refresh } = daemonConfig;
const defaultSize = 1024 * 1024; // default size for datasets with unknown size 1 MB

export const daemons: Daemon[] = [];

/**
 * Class representing a connection to a daemon server.
 */
class Daemon {
	private axiosInstance: AxiosCacheInstance;
	private cache: LRUCache<string, string>;
	requestLoad = 0;

	/**
	 * Create a Daemon instance.
	 * @param daemonUrl - URL of the daemon server
	 */
	constructor(daemonUrl: string) {
		this.axiosInstance = setupCache(
			axios.create({
				baseURL: daemonUrl,
				timeout,
				paramsSerializer: { indexes: null }
			}),
			{ ttl }
		);
		this.cache = new LRUCache<string, string>({
			maxSize,
			ttl,
			dispose: this.unload
		});
	}

	toString = () => `Daemon(${this.axiosInstance.getUri()})`;

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
			const { data } = await this.axiosInstance.get<Health>('/health', {
				cache: false,
				signal: AbortSignal.timeout(connectionTimeout)
			});
			console.log(`Health ${this}: status ${data.status}`);
			return data;
		} catch (error) {
			console.error(`Health ${this}: error ${error}`);
			return { status: 'unhealthy' };
		}
	};

	/**
	 * Get the list of datasets currently loaded in the daemon.
	 * @returns List of all loaded dataset ids
	 */
	loaded = async () => {
		try {
			const { data } = await this.axiosInstance.get<Loaded>('/loaded', { cache: false });
			const datasets = await this.datasets(data);
			for (const ds of data) this.cache.set(ds, ds, { size: datasets[ds].size ?? defaultSize });
			console.log(`Loaded ${this}: loaded ${data}`);
			return data;
		} catch (error) {
			console.error(`Loaded ${this}: error ${error}`);
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
			const { data } = await this.axiosInstance.post<Unload>('/unload', { ds }, { cache: false });
			this.cache.delete(ds);
			console.log(`Unload ${this} ds=${ds}: loaded ${data}`);
			return data;
		} catch (error) {
			console.error(`Unload ${this} ds=${ds}: error ${error}`);
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
			const { data } = await this.axiosInstance.post<Load>('/load', { ds }, { cache: false });
			const datasets = await this.datasets(data);
			for (const ds of data) this.cache.set(ds, ds, { size: datasets[ds].size ?? defaultSize });
			console.log(`Load ${this} ds=${ds}: loaded ${data}`);
			return data;
		} catch (error) {
			console.error(`Load ${this} ds=${ds}: error ${error}`);
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
			const { cached, data } = await this.axiosInstance.get<Datasets>('/datasets', {
				params: { ds }
			});
			if (!cached) console.log(`Datasets ${this} ds=${ds}: keys ${Object.keys(data)}`);
			return data;
		} catch (error) {
			console.error(`Datasets ${this} ds=${ds}: error ${error}`);
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
			const { cached, data } = await this.axiosInstance.get<Publications>('/publications', {
				params: { pub }
			});
			if (!cached) console.log(`Publications ${this} pub=${pub}: key ${Object.keys(data)}`);
			return data;
		} catch (error) {
			console.error(`Publications ${this} pub=${pub}: error ${error}`);
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
			const { data } = await this.axiosInstance.get<Genes>('/genes', { params: { ds } });
			console.log(`Genes ${this} ds=${ds}: keys ${Object.keys(data)}`);
			return data;
		} catch (error) {
			console.error(`Genes ${this} ds=${ds}: error ${error}`);
			return {};
		}
	};

	/**
	 * Fetch differentially expressed genes for datasets.
	 * @param ds - list of dataset ids to fetch
	 * @returns Dictionary mapping dataset ids to differentially expressed genes
	 */
	degs = async (ds: string[]): Promise<DEGs> => {
		try {
			const { data } = await this.axiosInstance.get<DEGs>('/degs', { params: { ds } });
			console.log(`DEGs ${this} ds=${ds}: keys ${Object.keys(data)}`);
			return data;
		} catch (error) {
			console.error(`DEGs ${this} ds=${ds}: error ${error}`);
			return {};
		}
	};

	/** Render plots.
	 * @param plotsParams - parameters for plot rendering
	 * @returns List of rendered plots paths
	 */
	plots = async (plotsParams: PlotsParams) => {
		this.requestLoad += maxSize;
		try {
			const { data } = await this.axiosInstance.post<Plots>('/plots', plotsParams);
			console.log(`Plots ${this} params ${Object.values(plotsParams)}: plots ${data}`);
			return data;
		} catch (error) {
			console.error(`Plots ${this} params ${Object.values(plotsParams)}: error ${error}`);
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
registerDaemons();

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
