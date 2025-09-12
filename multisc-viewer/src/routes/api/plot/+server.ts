import { datasets } from '$lib/server/data';
import { plot } from '$lib/server/plot';
import { Grouping, PlotType } from '$lib/types/plot';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const dsList = url.searchParams.getAll('ds');
	const geneList = url.searchParams.getAll('gene');
	const genes = geneList.length === 0 ? datasets[dsList[0]]?.defaultGenes : geneList;
	const groupBy = (url.searchParams.get('groupBy') ?? Grouping.CellType) as Grouping;
	const splitBy = (url.searchParams.get('splitBy') ?? Grouping.Genotype) as Grouping;
	const plots = (url.searchParams.getAll('pt') as PlotType[]) ?? Object.values(PlotType);
	const plotResults = plot({ datasets: dsList, genes, groupBy, splitBy, plotTypes: plots });
	return json(plotResults);
};

export const POST: RequestHandler = async ({ request }) => {
	const { datasets, genes, groupBy, splitBy, plots } = await request.json();
	const plotResults = plot({ datasets, genes, groupBy, splitBy, plotTypes: plots });
	return json(plotResults);
};
