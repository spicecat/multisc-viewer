import { datasets, getGenes, getDEGs, publications } from '$lib/server/data';
import { getPlotId, plot } from '$lib/server/plot';
import { Grouping, type PlotParams, PlotType } from '$lib/types/plot';
import { error, fail, redirect } from '@sveltejs/kit';
import type { Actions, PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ url, params }) => {
	const publication = publications[params.publicationId ?? ''];
	if (params.publicationId && !publication)
		error(404, `Publication ${params.publicationId} not found`);

	const dsList = url.searchParams.getAll('ds').filter((ds) => ds in datasets);
	const allGenes = await getGenes(dsList);
	const degs = await getDEGs(dsList);

	const datasetsData = dsList.map((ds) => datasets[ds]);
	const geneList = url.searchParams.getAll('gene');
	const genes = geneList.length === 0 ? datasetsData[0]?.defaultGenes : geneList;
	const groupBy = (url.searchParams.get('groupBy') ?? Grouping.CellType) as Grouping;
	const splitBy = groupBy === Grouping.CellType ? Grouping.Genotype : Grouping.CellType;
	const plotTypes = url.searchParams.getAll('pt') as PlotType[];
	if (plotTypes.length === 0) plotTypes.push(...Object.values(PlotType));

	const plotParams: PlotParams = { datasets: dsList, genes, groupBy, splitBy, plotTypes };

	const searchParams = new URLSearchParams();
	searchParams.delete('ds');
	dsList.forEach((ds) => searchParams.append('ds', ds));
	genes.forEach((gene) => searchParams.append('gene', gene));
	searchParams.set('groupBy', groupBy);
	searchParams.set('splitBy', splitBy);
	searchParams.delete('pt');
	plotTypes.forEach((pt) => searchParams.append('pt', pt));

	if (searchParams.toString() !== url.searchParams.toString())
		throw redirect(301, `${url.pathname}?${searchParams.toString()}`);

	const plotIds = Object.fromEntries(
		dsList.map((ds) => [
			ds,
			Object.fromEntries(
				genes.map((gene) => [
					gene,
					Object.fromEntries(plotTypes.map((pt) => [pt, getPlotId(ds, gene, groupBy, splitBy, pt)]))
				])
			)
		])
	);
	const plots = plot(plotParams);

	return {
		publication,
		datasets: datasetsData,
		genes: allGenes,
		degs,
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
		const genes = data.getAll('gene');
		const groupBy = data.get('groupBy');
		const plotTypes = data.getAll('pt');
		if (datasets.length === 0) return fail(400, { error: 'No datasets selected' });
		const searchParams = new URLSearchParams();
		datasets.forEach((ds) => searchParams.append('ds', ds.toString()));
		genes.forEach((gene) => searchParams.append('gene', gene.toString()));
		if (groupBy) searchParams.set('groupBy', groupBy.toString());
		plotTypes.forEach((pt) => searchParams.append('pt', pt.toString()));
		redirect(
			303,
			`/plot${publicationId ? '/' + publicationId.toString() : ''}?${searchParams.toString()}`
		);
	}
} satisfies Actions;
