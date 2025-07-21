// Shared types for the app
export interface Publication {
  PMID: number;
  title: string;
  // ...add other fields as needed
}

export interface Dataset {
  name: string;
  size: number;
  // ...add other fields as needed
}

export interface ChartResult {
  clustering: string;
  violin: string;
  feature: string;
}
