import { datasets, getGenes, publications } from '$lib/server/data';
import { getPlotId, plot } from '$lib/server/plot';
import { Grouping, PlotType } from '$lib/types/plot';
import { error, redirect } from '@sveltejs/kit';
import type { PageServerLoad } from '../$types';

export const load: PageServerLoad = async ({ url, params }) => {
	const publication = publications[params.publicationId];
	if (params.publicationId && !publication)
		error(404, `Publication ${params.publicationId} not found`);

	const dsList = url.searchParams.getAll('ds').filter((ds) => ds in datasets);
	const datasetsData = dsList.map((ds) => datasets[ds]);
	const genes = await getGenes(dsList);
	const paramGene = url.searchParams.get('gene') ?? '';
	const gene = genes.includes(paramGene) ? paramGene : (datasetsData[0]?.defaultGene ?? genes[0]);
	const groupBy = (url.searchParams.get('groupBy') ?? Grouping.CellType) as Grouping;
	const splitBy = (url.searchParams.get('splitBy') ?? Grouping.Genotype) as Grouping;
	const plotTypes = url.searchParams.getAll('pt') as PlotType[];
	if (plotTypes.length === 0) plotTypes.push(...Object.values(PlotType));

	const plotParams = { datasets: dsList, gene, groupBy, splitBy, plotTypes };

	const searchParams = new URLSearchParams();
	searchParams.delete('ds');
	dsList.forEach((ds) => searchParams.append('ds', ds));
	searchParams.set('gene', gene);
	searchParams.set('groupBy', groupBy);
	searchParams.set('splitBy', splitBy);
	searchParams.delete('pt');
	plotTypes.forEach((pt) => searchParams.append('pt', pt));

	if (searchParams.toString() !== url.searchParams.toString())
		throw redirect(303, `${url.pathname}?${searchParams.toString()}`);

	const plotIds = Object.fromEntries(
		dsList.map((ds) => [
			ds,
			Object.fromEntries(plotTypes.map((pt) => [pt, getPlotId(ds, gene, groupBy, splitBy, pt)]))
		])
	);
	const plots = plot(plotParams);

	return {
		publication,
		datasets: datasetsData,
		genes,
		plotParams,
		plotIds,
		plots,
		searchParams: searchParams.toString()
	};
};
