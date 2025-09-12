import type {
	DatasetsResponse,
	LoadResponse,
	PublicationsResponse,
	RenderResponse,
	StatusResponse,
	UnloadResponse
} from '$lib/types/daemon';
import type { PlotParams } from '$lib/types/plot';
import axios, { type AxiosInstance } from 'axios';
import { setupCache } from 'axios-cache-interceptor';
import NodeCache from 'node-cache';
import { daemonConfig } from './config';
import { datasets } from './data';

type Daemon = AxiosInstance;

const { daemons, connectionTimeout, timeout, stdTTL } = daemonConfig;

const datasetCache = new NodeCache({ stdTTL }); // Dataset to daemon
datasetCache.on('expired', (ds: string, daemon: Daemon) => _unloadDataset(ds, daemon));

const daemonLoad = new Map<Daemon, number>(); // Daemon to load

const removeDaemon = (daemon: Daemon) => {
	daemonLoad.delete(daemon);
	datasetCache.keys().forEach((ds) => {
		const d = datasetCache.get<Daemon>(ds);
		if (d === daemon) datasetCache.del(ds);
	});
	console.log(`Removed daemon ${daemon.getUri()}`);
};

const registerDaemons = async () =>
	await Promise.all(
		daemons.map(async (daemonUrl) => {
			const daemon: Daemon = setupCache(axios.create({ baseURL: daemonUrl, timeout }));
			try {
				const { data } = await daemon.get<StatusResponse>('/status', {
					timeout: connectionTimeout
				});
				const load = data.datasets.reduce((l: number, ds: string) => {
					datasetCache.set(ds, daemon);
					return datasets[ds].size + l;
				}, 0);
				daemonLoad.set(daemon, load);
				console.log(`Registered ${daemon.getUri()} with datasets ${data.datasets}, load ${load}.`);
			} catch (error) {
				console.warn(`Daemon ${daemon.getUri()} is not responding. ${error}`);
				removeDaemon(daemon);
			}
		})
	);

setInterval(registerDaemons, 10 * 60 * 1000); // refresh daemons every 10 minutes
registerDaemons();

export const getDatasets = async (): Promise<DatasetsResponse[]> => {
	const results = await Promise.allSettled(
		Array.from(daemonLoad.keys()).map(async (daemon) => {
			const { data } = await daemon.get<DatasetsResponse>('/datasets');
			return data;
		})
	);
	return results
		.filter((r): r is PromiseFulfilledResult<DatasetsResponse> => r.status === 'fulfilled')
		.map((r) => r.value);
};

export const getPublications = async (): Promise<PublicationsResponse[]> => {
	const results = await Promise.allSettled(
		Array.from(daemonLoad.keys()).map(async (daemon) => {
			const { data } = await daemon.get<PublicationsResponse>('/publications');
			return data;
		})
	);
	return results
		.filter((r): r is PromiseFulfilledResult<PublicationsResponse> => r.status === 'fulfilled')
		.map((r) => r.value);
};

const _unloadDataset = async (ds: string, daemon: Daemon) => {
	try {
		const { data } = await daemon.post<UnloadResponse>('/unload', { datasets: [ds] });
		const load = (daemonLoad.get(daemon) ?? 0) - datasets[ds].size;
		if (daemonLoad.has(daemon)) daemonLoad.set(daemon, load);
		console.log(`Unloaded ${ds} from ${daemon.getUri()}; datasets ${data.datasets}, load ${load}.`);
	} catch (error) {
		console.warn(`Failed to unload ${ds} from ${daemon.getUri()}: ${error}`);
	} finally {
		datasetCache.del(ds);
	}
};

export const loadDatasets = async (dsList: string[]): Promise<Map<Daemon, Set<string>>> => {
	dsList = dsList.filter((ds) => ds in datasets);

	const cachedDatasets = datasetCache.mget<Daemon>(dsList);

	const daemonDatasets = new Map<Daemon, Set<string>>();
	const daemonUnloadedDatasets = new Map<Daemon, Set<string>>();

	if (daemonLoad.size)
		dsList.forEach((ds) => {
			const cached = cachedDatasets[ds];
			const useCached = !!cached && daemonLoad.has(cached);
			const daemon = useCached
				? cached
				: [...daemonLoad.entries()].reduce((minDaemon, currDaemon) =>
						currDaemon[1] < minDaemon[1] ? currDaemon : minDaemon
					)[0]; // daemon with minimum load

			if (!daemonDatasets.has(daemon)) {
				daemonDatasets.set(daemon, new Set<string>());
				daemonUnloadedDatasets.set(daemon, new Set<string>());
			}
			daemonDatasets.get(daemon)!.add(ds);
			if (!useCached) daemonUnloadedDatasets.get(daemon)!.add(ds);

			const load = datasets[ds].size;
			daemonLoad.set(daemon, (daemonLoad.get(daemon) ?? 0) + load);
		});

	await Promise.all(
		Array.from(daemonUnloadedDatasets.entries()).map(async ([daemon, dsSet]) => {
			const dsList = Array.from(dsSet);
			try {
				const { data } = await daemon.post<LoadResponse>('/load', { datasets: dsList });
				for (const [ds, loaded] of Object.entries(data)) {
					if (loaded) datasetCache.set(ds, daemon);
					else {
						// rollback load increment
						daemonLoad.set(daemon, (daemonLoad.get(daemon) ?? 0) - datasets[ds].size);
						console.error(`Failed to load ${ds} on ${daemon.getUri()}`);
					}
				}
			} catch (error) {
				// rollback all increments for daemon
				for (const ds of dsList)
					daemonLoad.set(daemon, (daemonLoad.get(daemon) ?? 0) - datasets[ds].size);
				console.error(`Error loading datasets ${dsList} on ${daemon.getUri()}: ${error}`);
			}
		})
	);
	return daemonDatasets;
};

export const render = async (plotParams: PlotParams): Promise<RenderResponse> => {
	try {
		const daemonDatasets = await loadDatasets(plotParams.datasets);
		const plots = await Promise.all(
			Array.from(daemonDatasets.entries()).map(async ([daemon, dsSet]) => {
				const dsList = Array.from(dsSet);
				try {
					const { data } = await daemon.post<RenderResponse>('/render', {
						...plotParams,
						datasets: dsList
					});
					console.log(`Rendered ${dsList} on ${daemon.getUri()}`);
					return data;
				} catch (error) {
					console.error(`Error rendering ${dsList} on ${daemon.getUri()}: ${error}`);
					return [];
				}
			})
		);
		return plots.flat();
	} catch (error) {
		console.error(`Render pipeline error: ${error}`);
		return [];
	}
};
