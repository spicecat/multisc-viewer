import type { Author, Dataset, Gene, Publication } from '$lib/types/data';

export type GeneData = {
	id: Gene;
	[key: string]: Gene;
};

export type DatasetData = Omit<Dataset, 'size' | 'defaultGene'> & {
	[key: string]: string | number | string[] | Gene | Author;
};

export type PublicationData = Omit<Publication, 'datasets'> & {
	datasets: string[];
	href: string | URL;
	[key: string]: string | number | Author[] | string[] | URL;
};

export type Data = GeneData | DatasetData | PublicationData;

export type Columns = {
	key: string;
	label: string;
	href?: string;
}[];

export type Select = 'checkbox' | 'radio';

export const geneColumns: Columns = [{ key: 'id', label: 'Gene' }];

export const datasetColumns = [
	{ key: 'title', label: 'Dataset' },
	{ key: 'year', label: 'Year' },
	{ key: 'authors', label: 'Authors' },
	{ key: 'region', label: 'Region' },
	{ key: 'disease', label: 'Disease' },
	{ key: 'cellType', label: 'Cell Type' },
	{ key: 'PMID', label: 'PMID' }
];

export const publicationColumns: Columns = [
	{ key: 'title', label: 'Title', href: 'href' },
	{ key: 'year', label: 'Year' },
	{ key: 'authors', label: 'Authors' },
	{ key: 'journal', label: 'Journal' },
	{ key: 'PMID', label: 'PMID' },
	{ key: 'datasets', label: 'Datasets' }
];
