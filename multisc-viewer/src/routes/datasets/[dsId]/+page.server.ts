import { getDatasets } from '$lib/server/data';
import { error } from '@sveltejs/kit';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ params }) => {
	const datasets = await getDatasets([params.dsId]);
	const dataset = datasets[params.dsId];
	if (dataset) return { dataset };
	error(404, `Dataset ${params.dsId} not found`);
};
