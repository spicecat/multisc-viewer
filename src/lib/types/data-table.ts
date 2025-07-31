export type Columns = {
	key: string;
	label: string;
	sortable?: boolean;
	href?: string;
}[];

export const datasetColumns: Columns = [
	{ key: 'id', label: 'Dataset' },
	{ key: 'year', label: 'Year' },
	{ key: 'author', label: 'Author' },
	{ key: 'PMID', label: 'PMID' },
	{ key: 'region', label: 'Region' },
	{ key: 'disease', label: 'Disease' },
	{ key: 'cellType', label: 'Cell Type' }
];
