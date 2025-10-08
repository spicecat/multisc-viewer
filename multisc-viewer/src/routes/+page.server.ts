import { getDatasets, getPublications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	datasets: await getDatasets(),
	publications: await getPublications()
});
