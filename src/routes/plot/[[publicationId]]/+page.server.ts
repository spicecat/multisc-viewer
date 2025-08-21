import { datasets, getGenes, publications } from '$lib/server/data';
import { getPlotId, plot } from '$lib/server/plot';
import { Grouping, PlotType } from '$lib/types/plot';
import { error, fail, redirect } from '@sveltejs/kit';
import type { Actions, PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ url, params }) => {
	const publication = publications[params.publicationId ?? ''];
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
		throw redirect(301, `${url.pathname}?${searchParams.toString()}`);

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

export const actions = {
	plot: async ({ request }) => {
		const data = await request.formData();
		const publicationId = data.get('publicationId');
		const datasets = data.getAll('datasets');
		const gene = data.get('gene');
		const groupBy = data.get('groupBy');
		const splitBy = data.get('splitBy');
		const plotTypes = data.getAll('pt');
		if (datasets.length === 0) return fail(400, { error: 'No datasets selected' });
		const searchParams = new URLSearchParams();
		datasets.forEach((ds) => searchParams.append('ds', ds.toString()));
		if (gene) searchParams.set('gene', gene.toString());
		if (groupBy) searchParams.set('groupBy', groupBy.toString());
		if (splitBy) searchParams.set('splitBy', splitBy.toString());
		plotTypes.forEach((pt) => searchParams.append('pt', pt.toString()));
		if (publicationId)
			redirect(303, `/plot/${publicationId.toString()}?${searchParams.toString()}`);
		else redirect(303, `/plot?${searchParams.toString()}`);
	}
} satisfies Actions;
