import type { Dataset, Gene, Publication } from '$lib/types/data';

export type Data = Gene | Dataset | Publication;

export type Columns = {
	key: string | ((d: Data) => string);
	label: string;
	href?: (d: Data) => string;
}[];

export type Select = 'checkbox' | 'radio';

export const geneColumns: Columns = [{ key: '', label: 'Gene' }];

export const datasetColumns: Columns = [
	{ key: 'year', label: 'Year' },
	{ key: 'author', label: 'Author' },
	{ key: 'region', label: 'Region' },
	{ key: 'disease', label: 'Disease' },
	{ key: 'cellType', label: 'Cell Type' },
	{ key: 'PMID', label: 'PMID' }
];

export const publicationColumns: Columns = [
	{
		key: 'title',
		label: 'Title',
		href: (d: Data) => `/publication/${(d as Publication).id}`
	},
	{ key: 'year', label: 'Year' },
	{ key: 'author', label: 'Author' },
	{ key: 'journal', label: 'Journal' },
	{ key: 'PMID', label: 'PMID' },
	{ key: (d: Data) => `${(d as Publication).datasets.length}`, label: 'Datasets' }
];
