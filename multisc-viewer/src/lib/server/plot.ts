import type { PlotParams, PlotResults, PlotType } from '$lib/types/plot';
import { uniq } from 'lodash-es';
import NodeCache from 'node-cache';
import { mkdirSync, rmSync, watch } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { basename, dirname, join, resolve } from 'node:path';
import { plotConfig } from './config';
import { render } from './daemon';

type PlotRequest = { resolve: (plotPath: string) => void; reject: (error: Error) => void };

const { cacheDir, plotName, stdTTL, maxWaitTime } = plotConfig;

// clear cache
rmSync(cacheDir, { recursive: true, force: true });
mkdirSync(cacheDir, { recursive: true });

const plotCache = new NodeCache({ stdTTL });
plotCache.on('expired', async (_plotId: string, plotPath: string) =>
	rmSync(plotPath, { recursive: true, force: true })
);

const requestCache = new NodeCache({ stdTTL: maxWaitTime });
requestCache.on('expired', async (plotId: string, request: PlotRequest) =>
	request.reject(new Error(`Timeout waiting for plot ${plotId}`))
);

watch(cacheDir, { recursive: true }, async (event, plotPath) => {
	if (event !== 'rename' || !plotPath || basename(plotPath) !== plotName) return;
	const plotId = dirname(plotPath);
	plotPath = resolve(cacheDir, plotPath);
	plotCache.set<string>(plotId, plotPath);
	requestCache.get<PlotRequest>(plotId)?.resolve(plotPath);
});

export const getPlotId = (
	ds: string,
	gene: string,
	groupBy: string,
	splitBy: string,
	pt: PlotType
): string => join(ds, gene, `${groupBy}:${splitBy}`, pt);

const getPlotIds = ({ datasets, genes, groupBy, splitBy, plotTypes }: PlotParams): Set<string> =>
	new Set(
		datasets.flatMap((ds) =>
			genes.flatMap((gene) => plotTypes.map((pt) => getPlotId(ds, gene, groupBy, splitBy, pt)))
		)
	);

const plotIdParts = (plotId: string) => {
	const [ds, gene, grouping, pt] = plotId.split('/');
	const [groupBy, splitBy] = grouping.split(':');
	return { ds, gene, groupBy, splitBy, pt };
};

const queryPlots = (plotIds: Set<string>): PlotResults => {
	const queryPlot = async (plotId: string): Promise<string> => {
		const plotPath =
			plotCache.get<string>(plotId) ??
			(await new Promise<string>((resolve, reject) => {
				requestCache.set<PlotRequest>(plotId, { resolve, reject });
			}));
		plotIds.delete(plotId);
		const plotData = await readFile(plotPath, 'base64');
		return 'data:image/png;base64,' + plotData;
	};
	return Object.fromEntries([...plotIds].map((plotId) => [plotId, queryPlot(plotId)]));
};

export const plot = (plotParams: PlotParams): PlotResults => {
	const keys = getPlotIds(plotParams);
	const plots = queryPlots(keys);
	const missingDatasets = uniq([...keys].map((plotId) => plotIdParts(plotId).ds));
	render({ ...plotParams, datasets: missingDatasets })
		.then((renderedIds) => {
			const rendered = new Set(renderedIds);
			// reject pending requests if not included
			for (const plotId of keys) {
				if (!rendered.has(plotId)) {
					const req = requestCache.get<PlotRequest>(plotId);
					if (req) {
						requestCache.del(plotId);
						req.reject(new Error(`Daemon failed to render plot ${plotId}`));
					}
				}
			}
		})
		.catch((error) => {
			// reject all pending requests
			for (const plotId of keys) {
				const req = requestCache.get<PlotRequest>(plotId);
				if (req) {
					requestCache.del(plotId);
					req.reject(new Error(`Daemon render error for ${plotId}: ${String(error)}`));
				}
			}
		});
	return plots;
};
