import { getDatasets } from '$lib/server/data';
import { plots } from '$lib/server/plots';
import type { PlotsParams } from '$lib/types/daemon';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const ds = url.searchParams.getAll('ds');
	const genesParam = url.searchParams.getAll('gene');

	const gene = (
		genesParam.length === 0
			? Object.values(await getDatasets(ds)).map((d) => d.defaultGene)
			: genesParam
	).filter((g) => g !== undefined);
	const pt = <PlotsParams['pt']>url.searchParams.getAll('pt');
	if (pt.length === 0) pt.push(...(<PlotsParams['pt']>['umap', 'vln', 'feat']));
	const groupBy = <PlotsParams['groupBy']>url.searchParams.get('groupBy') ?? 'CellType';
	const splitBy = <PlotsParams['splitBy']>url.searchParams.get('splitBy') ?? 'Genotype';
	const plotsResults = plots({ ds, gene, pt, groupBy, splitBy });
	return json(plotsResults);
};

export const POST: RequestHandler = async ({ request }) => {
	const {
		ds,
		gene,
		pt,
		groupBy = 'CellType',
		splitBy = 'Genotype'
	}: PlotsParams = await request.json();
	const plotsResults = plots({ ds, gene, pt, groupBy, splitBy });
	return json(plotsResults);
};
