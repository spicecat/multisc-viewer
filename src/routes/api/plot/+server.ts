import { plot } from '$lib/server/plot';
import type { Grouping } from '$lib/types/plot';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async ({ url }) => {
	const datasets = url.searchParams.getAll('ds');
	const gene = url.searchParams.get('gene') ?? undefined;
	const groupBy = (url.searchParams.get('groupBy') ?? undefined) as Grouping;
	const splitBy = (url.searchParams.get('splitBy') ?? undefined) as Grouping;
	const plotResults = await plot({ datasets, gene, groupBy, splitBy });
	return json(plotResults);
};

export const POST: RequestHandler = async ({ request }) => {
	const { datasets, gene, groupBy, splitBy } = await request.json();
	const plotResults = await plot({ datasets, gene, groupBy, splitBy });
	return json(plotResults);
};
