import type {
	LoadResponse,
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

const { daemons, timeout, stdTTL } = daemonConfig;

const datasetCache = new NodeCache({ stdTTL }); // Dataset to daemon
datasetCache.on('expired', (ds: string, daemon: Daemon) => _unloadDataset(ds, daemon));

const daemonLoad = new Map<Daemon, number>(); // Daemon to load

// register daemons
await Promise.all(
	daemons.map(async (daemonUrl) => {
		const daemon: Daemon = setupCache(axios.create({ baseURL: daemonUrl }));
		try {
			const { data } = await daemon.get<StatusResponse>('/status', { timeout });
			const load = data.datasets.reduce((l: number, ds: string) => {
				datasetCache.set(ds, daemon);
				return datasets[ds].size + l;
			}, 0);
			daemonLoad.set(daemon, load);
			console.log(`Registered ${daemon.getUri()} with datasets ${data.datasets}, load ${load}.`);
		} catch (error) {
			console.log(`Daemon ${daemon.getUri()} is not responding. ${error}`);
			daemonLoad.set(daemon, 0);
		}
	})
);
const _unloadDataset = async (ds: string, daemon: Daemon) => {
	const { data } = await daemon.post<UnloadResponse>('/unload', { datasets: [ds] });
	datasetCache.del(ds);
	const load = daemonLoad.get(daemon)! - datasets[ds].size;
	daemonLoad.set(daemon, load);
	console.log(`Unloaded ${ds} from ${daemon.getUri()}; datasets ${data.datasets}, load ${load}.`);
};

export const loadDatasets = async (dsList: string[]): Promise<Map<Daemon, Set<string>>> => {
	if (daemonLoad.size === 0) throw new Error('No daemons registered');
	dsList = dsList.filter((ds) => ds in datasets);

	const cachedDatasets = datasetCache.mget<Daemon>(dsList);

	const daemonDatasets = new Map<Daemon, Set<string>>();
	const daemonUnloadedDatasets = new Map<Daemon, Set<string>>();
	dsList.forEach((ds) => {
		const daemon =
			cachedDatasets[ds] ??
			[...daemonLoad.entries()].reduce((minDaemon, currDaemon) =>
				currDaemon[1] < minDaemon[1] ? currDaemon : minDaemon
			)[0]; // daemon with minimum load

		if (!daemonDatasets.has(daemon)) {
			daemonDatasets.set(daemon, new Set<string>());
			daemonUnloadedDatasets.set(daemon, new Set<string>());
		}
		daemonDatasets.get(daemon)!.add(ds);
		if (!(ds in cachedDatasets)) daemonUnloadedDatasets.get(daemon)!.add(ds);

		const load = datasets[ds].size;
		daemonLoad.set(daemon, daemonLoad.get(daemon)! + load);
	});

	await Promise.all(
		daemonUnloadedDatasets.entries().map(async ([daemon, dsSet]) => {
			const dsList = Array.from(dsSet);
			const { data } = await daemon.post<LoadResponse>('/load', { datasets: dsList });
			for (const [ds, loaded] of Object.entries(data)) {
				if (loaded) datasetCache.set(ds, daemon);
				else {
					daemonLoad.set(daemon, daemonLoad.get(daemon)! - datasets[ds].size);
					console.error(`Failed to load ${ds} on ${daemon.getUri()}`);
					throw new Error(`Failed to load ${ds} on ${daemon.getUri()}`);
				}
			}
		})
	);
	return daemonDatasets;
};

export const render = async (plotParams: PlotParams): Promise<RenderResponse> => {
	const daemonDatasets = await loadDatasets(plotParams.datasets);
	const plots = await Promise.all(
		daemonDatasets.entries().map(async ([daemon, dsSet]) => {
			const dsList = Array.from(dsSet);
			const { data } = await daemon.post<RenderResponse>('/render', {
				...plotParams,
				datasets: dsList
			});
			console.log(`Rendered ${dsList} on ${daemon.getUri()}`);
			return data;
		})
	);
	return plots.flat();
};
