import type { Actions, PageServerLoad } from './$types';
import { fail, redirect } from '@sveltejs/kit';
import { datasets, publications } from '$lib/server/data';

export const load: PageServerLoad = async () => ({
	datasets: Array.from(datasets.values()),
	publications: Array.from(publications.values())
});

export const actions = {
	plot: async ({ request }) => {
		const data = await request.formData();
		const selected = data.getAll('selected');
		if (selected.length === 0) return fail(400, { error: 'No datasets selected' });
		const urlParams = new URLSearchParams();
		selected.forEach((ds) => urlParams.append('ds', ds.toString()));
		redirect(303, `/plot?${urlParams.toString()}`);
	}
} satisfies Actions;
