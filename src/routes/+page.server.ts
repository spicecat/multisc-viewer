import { datasets, publications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	datasets: Array.from(datasets.values()),
	publications: Array.from(publications.values())
});

export const actions = {
	plot: async ({ request }) => {
		const data = await request.formData();
		console.log(data);
	}
};
