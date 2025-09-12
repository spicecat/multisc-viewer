export enum Grouping {
	Genotype = 'Genotype',
	CellType = 'CellType'
}
export enum PlotType {
	UMAP = 'umap',
	Violin = 'vln',
	Feature = 'feat'
}

export interface PlotParams {
	datasets: string[];
	genes: string[];
	groupBy: Grouping;
	splitBy: Grouping;
	plotTypes: PlotType[];
}

export type PlotResults = Record<string, Promise<string>>; // Base64 encoded image
