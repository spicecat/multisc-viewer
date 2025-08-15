import { publications } from '$lib/server/data';
import { error, fail, redirect } from '@sveltejs/kit';
import type { Actions, PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ params }) => {
	const publication = publications[params.publicationId];
	if (publication) return { publication };
	error(404, `Publication ${params.publicationId} not found`);
};

export const actions = {
	plot: async ({ request, params }) => {
		const data = await request.formData();
		const selected = data.getAll('selected');
		console.log(data, selected);
		if (selected.length === 0) return fail(400, { error: 'No datasets selected' });
		const urlParams = new URLSearchParams();
		selected.forEach((ds) => urlParams.append('ds', ds.toString()));
		redirect(303, `/plot/${params.publicationId}?${urlParams.toString()}`);
	}
} satisfies Actions;
