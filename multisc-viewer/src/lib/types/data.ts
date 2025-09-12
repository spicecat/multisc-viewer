export type Author = string;
export type Gene = string;
export type DEGs = Record<string, Gene[]>;

export interface Dataset {
	id: string;
	title: string;
	year: number;
	authors: Author[];
	PMID: string[];
	region: string[];
	disease: string[];
	cellType: string;
	defaultGenes: Gene[];
	size: number;
}

export interface Publication {
	id: string;
	title: string;
	year: number;
	authors: Author[];
	journal: string;
	PMID: string;
	abstract: string;
	datasets: Dataset[];
}
