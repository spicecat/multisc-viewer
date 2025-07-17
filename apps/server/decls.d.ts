interface Author {
  name: string;
  affiliation: string;
}

interface Publication {
  PMID: number;
  title: string;
  abstract: string;
  authors: Author[];
  journal: string;
  doi: string;
  year: number;
  datasets: string[];
}

interface Dataset {
  name: string;
  year: number;
  region: string[];
  PMID: string;
  species: string;
  author: string;
  disease: string[];
  size: number;
  cellType: string;
  defaultGene: string;
}

interface ChartResult {
  clustering: string;
  violin: string;
  feature: string;
}

type IndexProps = {
  datasets: Dataset[];
  token: string;
};

type CompareProps = {
  order: string[];
  genes: string[];
  gene: string;
};
