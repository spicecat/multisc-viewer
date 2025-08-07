import { datasets, getGenes } from '$lib/server/data';
import { plot } from '$lib/server/plot';
import type { Grouping } from '$lib/types/plot';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ url }) => {
	const ds = url.searchParams.getAll('ds');
	const datasetsData = ds.map((id) => datasets.get(id)!).filter(Boolean);
	const genes = await getGenes(ds);
	const plotParams = {
		datasets: ds,
		gene: url.searchParams.get('gene') ?? datasetsData[0]?.defaultGene ?? genes[0],
		groupBy: (url.searchParams.get('groupBy') ?? 'CellType') as Grouping,
		splitBy: (url.searchParams.get('splitBy') ?? 'Genotype') as Grouping
	};
	url.searchParams.set('gene', plotParams.gene);
	url.searchParams.set('groupBy', plotParams.groupBy);
	url.searchParams.set('splitBy', plotParams.splitBy);
	return {
		datasets: datasetsData,
		plots: plot(plotParams),
		genes
	};
};
