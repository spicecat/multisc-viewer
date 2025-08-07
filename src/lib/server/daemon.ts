import { daemonConfig } from '$lib/config';
import type { PlotParams } from '$lib/types/plot';
import axios from 'axios';
import NodeCache from 'node-cache';
import { datasets } from './data';

type Daemon = string;

const { server, ports, stdTTL } = daemonConfig;

const daemonLoad = new Map<Daemon, number>(); // Daemon to load
const datasetCache = new NodeCache({ stdTTL }); // Dataset to daemon
datasetCache.on('expired', (ds: string, daemon: Daemon) => _unloadDataset(ds, daemon));

await Promise.all(
	ports.map(async (port) => {
		const daemon: Daemon = `${server}:${port}`;
		try {
			const response = await axios.get(`${daemon}/status`);
			const load = response.data.datasets.reduce((l: number, ds: string) => {
				datasetCache.set(ds, daemon);
				return datasets.get(ds)!.size + l;
			}, 0);
			daemonLoad.set(daemon, load);
		} catch (error) {
			console.log(`Daemon ${daemon} is not responding. ${error}`);
		}
	})
);

const _unloadDataset = async (ds: string, daemon: Daemon) => {
	await axios.post(`${daemon}/unload`, { ds });
	datasetCache.del(ds);
	daemonLoad.set(daemon, daemonLoad.get(daemon)! - datasets.get(ds)!.size);
};

export const loadDataset = async (ds: string): Promise<Daemon> => {
	const cacheDaemon = datasetCache.get<Daemon>(ds);
	if (cacheDaemon) return cacheDaemon;
	if (daemonLoad.size === 0) throw new Error('No daemons available');
	const [daemon] = [...daemonLoad.entries()].reduce((minDaemon, currDaemon) =>
		currDaemon[1] < minDaemon[1] ? currDaemon : minDaemon
	);
	await axios.post(`${daemon}/load`, { ds });
	datasetCache.set(ds, daemon);
	const load = datasets.get(ds)!.size;
	daemonLoad.set(daemon, daemonLoad.get(daemon)! + load);
	return daemon;
};

export const render = async (params: PlotParams): Promise<string> => {
	const daemon = await loadDataset(params.ds);
	const response = await axios.post(`${daemon}/render`, params);
	console.log(response);
	return response.data.key;
};
