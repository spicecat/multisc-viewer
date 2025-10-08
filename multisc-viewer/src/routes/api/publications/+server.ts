import { getPublications } from '$lib/server/data';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const pub = url.searchParams.getAll('pub');
	return json(await getPublications(pub));
};
