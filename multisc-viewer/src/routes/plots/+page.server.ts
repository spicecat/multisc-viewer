import { getDatasets, getPublications, getGenesRows } from '$lib/server/data';
import { getPlotId, plots } from '$lib/server/plots';
import type { Datasets, PlotsParams, Publication } from '$lib/types/daemon';
import { error, redirect } from '@sveltejs/kit';
import { uniq } from 'lodash-es';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ url }) => {
	let publication: Publication | undefined, datasets: Datasets | undefined;

	const pub = url.searchParams.get('pub');
	const ds = uniq(url.searchParams.getAll('ds'));
	const gene = uniq(url.searchParams.getAll('gene'));
	const pt = uniq(<PlotsParams['pt']>url.searchParams.getAll('pt'));
	const groupBy = <PlotsParams['groupBy']>url.searchParams.get('groupBy') ?? 'CellType';
	const splitBy = <PlotsParams['splitBy']>url.searchParams.get('splitBy') ?? 'Genotype';

	if (gene.length === 0) {
		if (pub) publication = (await getPublications([pub]))[pub];
		if (pub && !publication) error(404, `Publication ${pub} not found`);
		datasets = await getDatasets(publication ? publication.datasets : undefined);

		const genesSet = new Set<string>();
		for (const d of ds) {
			const defaultGene = datasets[d]?.defaultGene;
			if (defaultGene) genesSet.add(defaultGene);
		}
		gene.push(...genesSet);
	}
	if (pt.length === 0) pt.push(...(<PlotsParams['pt']>['umap', 'vln', 'feat']));

	const searchParams = new URLSearchParams();
	if (pub) searchParams.set('pub', pub);
	searchParams.delete('ds');
	for (const d of ds) searchParams.append('ds', d);
	searchParams.delete('gene');
	for (const g of gene) searchParams.append('gene', g);
	searchParams.delete('pt');
	for (const p of pt) searchParams.append('pt', p);
	searchParams.set('groupBy', groupBy);
	searchParams.set('splitBy', splitBy);

	if (searchParams.toString() !== url.searchParams.toString())
		throw redirect(301, `${url.pathname}?${searchParams.toString()}`);

	if (pub && !publication) publication = (await getPublications([pub]))[pub];
	if (pub && !publication) error(404, `Publication ${pub} not found`);
	if (!datasets) datasets = await getDatasets(publication ? publication.datasets : ds);

	const plotsIdMap = Object.fromEntries(
		ds.map((d) => [
			d,
			Object.fromEntries(
				gene.map((g) => [
					g,
					Object.fromEntries(pt.map((p) => [p, getPlotId(d, g, p, groupBy, splitBy)]))
				])
			)
		])
	);

	const plotsParams: PlotsParams = { ds, gene, pt, groupBy, splitBy };
	const plotsResults = plots(plotsParams);

	const genesRows = await getGenesRows(Object.keys(datasets));

	const name = publication
		? (publication.name ?? publication._id)
		: `datasets ${Object.values(datasets)
				.map((d) => d.name ?? d._id)
				.join(', ')}`;
	return {
		datasets,
		publication,
		genesRows,
		plotsIdMap,
		plotsParams,
		plotsResults,
		meta: {
			title: `MultiSC-Viewer - Plots for ${name}`,
			description: `Explore visualizations for ${name} with MultiSC-Viewer, a web application for visualizing and comparing gene expression across multiple datasets, brain regions, disease conditions, and species.`
		}
	};
};
