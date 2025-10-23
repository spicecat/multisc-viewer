import type { PlotRequests, PlotsParams } from '$lib/types/daemon';
import { uniq } from 'lodash-es';
import { LRUCache } from 'lru-cache';
import { mkdirSync, rmSync, statSync, watch } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { basename, dirname, join, resolve } from 'node:path';
import { plotsConfig } from './config';
import { getDaemonTargets } from './daemon';

type PlotRequest = { resolve: (plotPath: string) => void; reject: (reason: string) => void };

const { plotsDir, timeout, maxSize } = plotsConfig;
const plotName = 'plot.png';

// clear plots directory
rmSync(plotsDir, { force: true, recursive: true });
mkdirSync(plotsDir, { recursive: true });

/** Cache mapping requested plot ids to plot requests */
const requestCache = new LRUCache<string, PlotRequest>({
	ttl: timeout,
	ttlAutopurge: true,
	dispose: (request, plotId, reason) => {
		if (reason === 'evict') request.reject(`Timeout waiting for plot ${plotId}`);
	}
});

/** Cache mapping rendered plot ids to path to directory with plot */
const plotsCache = new LRUCache<string, string>({
	maxSize,
	sizeCalculation: (path) => {
		try {
			return Math.max(statSync(path).size, 1);
		} catch {
			return 1;
		}
	},
	dispose: (path) => rmSync(path, { force: true, recursive: true })
});

watch(plotsDir, { recursive: true }, async (event, plotPath) => {
	if (event !== 'rename' || !plotPath || basename(plotPath) !== plotName) return;
	const plotId = dirname(plotPath);
	plotPath = resolve(plotsDir, plotPath);
	requestCache.get(plotId)?.resolve(plotPath);
	requestCache.delete(plotId);
	plotsCache.set(plotId, plotPath);
});

/**
 * Get plot id from plot parameters.
 * @param ds Dataset id
 * @param gene Gene
 * @param pt Plot type
 * @param groupBy Grouping variable
 * @param splitBy Splitting variable
 * @returns Plot id
 */
export const getPlotId = (
	d: string,
	g: string,
	p: PlotsParams['pt'][number],
	groupBy: PlotsParams['groupBy'],
	splitBy: PlotsParams['splitBy']
): string => join(d, g, p, `${groupBy}:${splitBy}`);

/**
 * Get plot ids from plot parameters.
 * @param plotsParams Plot parameters
 * @returns Set of plot ids
 */
const getPlotIds = ({ ds, gene, pt, groupBy, splitBy }: PlotsParams) =>
	new Set(
		ds.flatMap((d) => gene.flatMap((g) => pt.map((p) => getPlotId(d, g, p, groupBy, splitBy))))
	);

/**
 * Get plot parameters from plot id.
 * @param plotId Plot id
 * @returns Plot parameters
 */
const getPlotParams = (plotId: string) => {
	const [ds, gene, pt, grouping] = plotId.split('/');
	const [groupBy, splitBy] = grouping.split(':');
	return { ds, gene, pt, groupBy, splitBy };
};

/**
 * Wait for plots by plot ids. Plot ids are removed from the set as they are found.
 * @param plotIds Set of plot ids
 * @returns Dictionary mapping plot ids to promises resolving to paths to plots
 */
const queryPlots = (plotIds: Set<string>): PlotRequests => {
	const queryPlot = async (plotId: string): Promise<string> => {
		const plotPath =
			plotsCache.get(plotId) ??
			(await new Promise<string>((resolve, reject) => {
				requestCache.set(plotId, {
					resolve,
					reject: (reason: string) => {
						console.warn(reason);
						reject(reason);
					}
				});
			}));
		plotIds.delete(plotId);
		const plotData = await readFile(plotPath, 'base64');
		return 'data:image/png;base64,' + plotData;
	};
	return Object.fromEntries(Array.from(plotIds).map((plotId) => [plotId, queryPlot(plotId)]));
};

/** Request rendering of plots by parameters.
 * @param plotsParams Plot parameters
 * @returns Dictionary mapping plot ids to promises resolving to paths to plots
 */
export const plots = (plotsParams: PlotsParams): PlotRequests => {
	if (!plotsParams.ds.length || !plotsParams.gene.length || !plotsParams.pt.length) return {};

	plotsParams.ds = uniq(plotsParams.ds);
	plotsParams.gene = uniq(plotsParams.gene);
	plotsParams.pt = uniq(plotsParams.pt);
	const keys = getPlotIds(plotsParams);
	const plotsRequests = queryPlots(keys);
	const missingDatasets = uniq(Array.from(keys).map((plotId) => getPlotParams(plotId).ds));

	getDaemonTargets(missingDatasets, true).then(async (daemonTargets) => {
		for (const [daemon, ds] of daemonTargets)
			daemon.plots({ ...plotsParams, ds }).then((plots) =>
				getPlotIds({ ...plotsParams, ds }).forEach((plotId) => {
					if (!plots.includes(plotId)) {
						requestCache.get(plotId)?.reject(`Daemon failed to render plot ${plotId}`);
						requestCache.delete(plotId);
					}
				})
			);
	});

	return plotsRequests;
};
