import { getPublications, getDatasets } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => {
	const publications = await getPublications();
	const ds = Object.values(publications).flatMap((pub) => pub.datasets);
	const datasets = await getDatasets(ds);
	return {
		datasets,
		publications
	};
};
