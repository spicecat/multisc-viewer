import type { Actions, PageServerLoad } from './$types';
import { fail, redirect } from '@sveltejs/kit';
import { publications } from '$lib/server/data';
import { error } from '@sveltejs/kit';

export const load: PageServerLoad = async ({ params }) => {
	const publication = publications.get(params.slug);
	if (publication) return { publication };
	error(404, `Publication ${params.slug} not found`);
};

export const actions = {
	plot: async ({ request }) => {
		const data = await request.formData();
		const selected = data.getAll('selected');
		console.log(data, selected);
		if (selected.length === 0) return fail(400, { error: 'No datasets selected' });
		const urlParams = new URLSearchParams();
		selected.forEach((ds) => urlParams.append('ds', ds.toString()));
		redirect(303, `/plot?${urlParams.toString()}`);
	}
} satisfies Actions;
