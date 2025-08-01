export type Grouping = 'Genotype' | 'CellType';

export interface PlotParams {
	ds: string;
	gene?: string;
	groupBy?: Grouping;
	splitBy?: Grouping;
}

export type PlotsParams = Omit<PlotParams, 'ds'> & {
	datasets: string[];
};

export interface Plot {
	clustering: string;
	violin: string;
	feature: string;
}

export type PlotResults = {
	[ds: string]: Plot;
};
