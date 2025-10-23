import { getDatasets, getPublications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	datasets: await getDatasets(),
	publications: await getPublications(),
	meta: {
		title: 'MultiSC-Viewer - Home',
		description:
			'MultiSC-Viewer is a web application for visualizing and comparing gene expression in multiple single cell/nucleus datasets across different brain regions, disease conditions, and species.'
	}
});
