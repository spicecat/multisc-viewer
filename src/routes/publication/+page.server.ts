import { publications } from '$lib/server/data';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => ({
	publications: Array.from(publications.values())
});
