import { getDEGs } from '$lib/server/data';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const datasets = url.searchParams.getAll('ds');
	return json(await getDEGs(datasets));
};
