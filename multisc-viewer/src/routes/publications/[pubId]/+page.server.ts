import { getDatasets, getPublications } from '$lib/server/data';
import { error } from '@sveltejs/kit';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ params }) => {
	const publications = await getPublications([params.pubId]);
	const publication = publications[params.pubId];
	if (!publication) error(404, `Publication ${params.pubId} not found`);
	const datasets = publication.datasets.length === 0 ? {} : await getDatasets(publication.datasets);
	const name = publication.name ?? publication._id;
	return {
		datasets,
		publication,
		meta: {
			title: `MultiSC-Viewer - Publication ${name}`,
			description: `Explore publication "${name}" with MultiSC-Viewer, a web application for visualizing and comparing gene expression across multiple datasets, brain regions, disease conditions, and species.`
		}
	};
};
