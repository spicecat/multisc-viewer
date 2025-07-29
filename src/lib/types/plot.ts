export type Grouping = 'Genotype' | 'CellType';

export type PlotParams = {
	ds: string;
	gene: string;
	groupBy: Grouping;
	splitBy: Grouping;
};

export interface PlotResult {
	clustering: string;
	violin: string;
	feature: string;
}
