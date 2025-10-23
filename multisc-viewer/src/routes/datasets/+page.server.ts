import { getDatasets } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	datasets: await getDatasets(),
	meta: {
		title: 'MultiSC-Viewer - Datasets',
		description:
			'Explore datasets with MultiSC-Viewer, a web application for visualizing and comparing gene expression across multiple datasets, brain regions, disease conditions, and species.'
	}
});
