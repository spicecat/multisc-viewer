import { error } from '@sveltejs/kit';
import type { PageServerLoad } from './$types';

const SPEC_URL =
	'https://git.jasonxu.dev/JasonXu/plot-viewer/raw/branch/main/multisc-daemon/openapi.json';

export const load: PageServerLoad = async ({ fetch }) => {
	const res = await fetch(SPEC_URL, { headers: { accept: 'application/json' } });
	if (!res.ok) {
		throw error(res.status, `Failed to fetch OpenAPI spec.`);
	}
	const spec = await res.json();
	return { spec };
};
