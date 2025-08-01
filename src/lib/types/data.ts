export type Author = string;
export type Gene = string;

export interface Dataset {
	id: string;
	year: number;
	author: Author;
	PMID: string[];
	region: string[];
	disease: string[];
	cellType: string;
	size: number;
	defaultGene?: Gene;
	[key: string]: string | number | string[] | Gene | Author | undefined;
}

export interface Publication {
	id: string;
	title: string;
	year: number;
	author: Author[];
	journal: string;
	PMID: string;
	abstract: string;
	datasets: Dataset[];
	[key: string]: string | number | Author[] | Dataset[] | undefined;
}
