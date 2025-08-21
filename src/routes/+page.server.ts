import { datasets, publications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	datasets: Object.values(datasets),
	publications: Object.values(publications)
});
