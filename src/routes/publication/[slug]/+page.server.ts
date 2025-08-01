import { publications } from '$lib/server/data';
import { error } from '@sveltejs/kit';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ params }) => {
	const publication = publications.get(params.slug);
	if (publication) return { publication };
	error(404, `Publication ${params.slug} not found`);
};
