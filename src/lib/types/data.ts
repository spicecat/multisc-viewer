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
}

export interface Publication {
	id: string;
	year: number;
	author: Author[];
	PMID: string;
	journal: string;
	title: string;
	abstract: string;
	datasets: Dataset[];
}

export type Data = Gene | Dataset | Publication;
