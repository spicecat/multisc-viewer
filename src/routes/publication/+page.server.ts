import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ fetch }) => {
	const res = await fetch('/api/publications');
	if (!res.ok) throw new Error('Failed to load publications');
	const publications = await res.json();
	return { publications };
};
