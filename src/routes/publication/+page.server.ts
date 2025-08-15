import { publications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const prerender = true;

export const load: PageServerLoad = async () => ({
	publications: Object.values(publications)
});
