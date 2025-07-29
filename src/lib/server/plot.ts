import { plotConfig } from '$lib/config';
import type { PlotParams, PlotResult } from '$lib/types/plot';
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

const _cacheKey = ({ ds, gene, groupBy, splitBy }: PlotParams): string =>
	`${ds}:${gene}:${groupBy}:${splitBy}`;

export const plot = async (params: PlotParams): Promise<PlotResult> => {
	const key = _cacheKey(params);
	const cachedResult = plotCache.get<PlotResult>(key);
	if (cachedResult) return cachedResult;
	await render(params);
	const clusteringPath = `${cacheDir}/${params.ds}/${key}/umap.png`;
	const violinPath = `${cacheDir}/${params.ds}/${key}/vln.png`;
	const featurePath = `${cacheDir}/${params.ds}/${key}/feat.png`;
	const clustering = 'data:image/png;base64,' + (await readFile(clusteringPath, 'base64'));
	const violin = 'data:image/png;base64,' + (await readFile(violinPath, 'base64'));
	const feature = 'data:image/png;base64,' + (await readFile(featurePath, 'base64'));
	const result = { clustering, violin, feature };
	plotCache.set(key, result);
	return result;
};
