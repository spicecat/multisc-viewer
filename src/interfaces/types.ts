export interface ChartResult {
  clustering: string;
  violin: string;
}

export interface Dataset {
  name: string;
  year: number;
  region: string[];
  PMID: string;
  species: string;
  author: string;
  disease: string[];
  size: number;
  cellType: string;
}

export interface Publication {
  publicationId: string;
  name: string;
  description: string;
  datasets: Dataset[];
}