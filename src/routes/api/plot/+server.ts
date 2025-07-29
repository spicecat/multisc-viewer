import { type RequestHandler, json } from '@sveltejs/kit';
import { render } from '$lib/server/plot';

export const GET: RequestHandler = async ({ url }) => {
	const dataset = url.searchParams.get('dataset');
	const gene = url.searchParams.get('gene');
	const groupBy = url.searchParams.get('groupBy');
	const splitBy = url.searchParams.get('splitBy');
	const render = await render(dataset, gene, groupBy, splitBy);
	return json(render);
};
