export function makeTitle(dataset: string): string {
	const [species, , author, disease, abbrev] = dataset.split('_');

	return `${species} ${author} ${disease} ${abbrev}`;
}
