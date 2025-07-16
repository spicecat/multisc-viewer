type Author = {
  name: string;
  affiliation: string;
};

type Publication = {
  PMID: number;
  title: string;
  abstract: string;
  authors: Author[];
  journal: string;
  doi: string;
  year: number;
  datasets: string[];
};

type Dataset = {
  name: string;
  year: number;
  region: string[];
  PMID: string;
  species: string;
  author: string;
  disease: string[];
  size: number;
  cellType: string;
};

type ChartResult = {
  key: string;
  path: string;
};

type IndexProps = {
  datasets: Dataset[];
  token: string;
};

type CompareProps = {
  order: string[];
  genes: string[];
  gene: string;
};
