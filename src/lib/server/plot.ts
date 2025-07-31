import { plotConfig } from '$lib/config';
import type { Plot, PlotResults, PlotsParams } from '$lib/types/plot';
import NodeCache from 'node-cache';
import { existsSync, mkdirSync, rmSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { render } from './daemon';

const { cacheDir, stdTTL } = plotConfig;

if (existsSync(cacheDir)) rmSync(cacheDir, { recursive: true, force: true });
mkdirSync(cacheDir, { recursive: true });

const plotCache = new NodeCache({ stdTTL });
plotCache.on('expired', (key: string) => {
	const [dataset, , ,] = key.split(':');
	const dir = `${cacheDir}/${dataset}/${key}`;
	if (existsSync(dir)) rmSync(dir, { recursive: true, force: true });
});

const _cacheKeys = ({ datasets, gene, groupBy, splitBy }: PlotsParams): string[] =>
	datasets.map((ds) => `${ds}:${gene}:${groupBy}:${splitBy}`);

export const plot = async (params: PlotsParams): Promise<PlotResults> => {
	const keys = _cacheKeys(params);
	const plots = plotCache.mget<Plot>(keys);
	await Promise.all(
		params.datasets.map(async (ds) => {
			if (plots[ds]) return;
			const key = await render({ ds, ...params });
			const clusteringPath = `${cacheDir}/${ds}/${key}/umap.png`;
			const violinPath = `${cacheDir}/${ds}/${key}/vln.png`;
			const featurePath = `${cacheDir}/${ds}/${key}/feat.png`;
			const clustering = 'data:image/png;base64,' + (await readFile(clusteringPath, 'base64'));
			const violin = 'data:image/png;base64,' + (await readFile(violinPath, 'base64'));
			const feature = 'data:image/png;base64,' + (await readFile(featurePath, 'base64'));
			plots[ds] = { clustering, violin, feature };
		})
	);
	return plots;
};
