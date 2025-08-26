import { datasets, getGenes } from '$lib/server/data';
import { plot } from '$lib/server/plot';
import { Grouping, PlotType } from '$lib/types/plot';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const dsList = url.searchParams.getAll('ds');
	const gene =
		url.searchParams.get('gene') ?? datasets[dsList[0]]?.defaultGene ?? (await getGenes(dsList))[0];
	const groupBy = (url.searchParams.get('groupBy') ?? Grouping.CellType) as Grouping;
	const splitBy = (url.searchParams.get('splitBy') ?? Grouping.Genotype) as Grouping;
	const plots = (url.searchParams.getAll('pt') as PlotType[]) ?? Object.values(PlotType);
	const plotResults = await plot({ datasets: dsList, gene, groupBy, splitBy, plotTypes: plots });
	return json(plotResults);
};

export const POST: RequestHandler = async ({ request }) => {
	const { datasets, gene, groupBy, splitBy, plots } = await request.json();
	const plotResults = await plot({ datasets, gene, groupBy, splitBy, plotTypes: plots });
	return json(plotResults);
};
