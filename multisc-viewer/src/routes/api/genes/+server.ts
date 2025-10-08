import { getGenes } from '$lib/server/data';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const ds = url.searchParams.getAll('ds');
	return json(await getGenes(ds));
};
