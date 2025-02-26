export function makeTitle(dataset: string): string {
	const [species, , author, disease, abbrev] = dataset.split('_').slice(0, -1);

	return `${species} ${author} ${disease} ${abbrev}`;
}

