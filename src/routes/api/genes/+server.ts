import { getGenes } from '$lib/server/data';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const datasets = url.searchParams.get('datasets')?.split(',') || [];
	return json(getGenes(datasets));
};
