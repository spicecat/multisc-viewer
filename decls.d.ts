type Publication = {
  publicationId: string;
  title: string;
  journal: string;
  year: number;
  authors: string;
  datasets: string[];
};

type Dataset = {
  id: string;
  name: string;
  size: number;
  publicationId: string;
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
