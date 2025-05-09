/// <reference types="svelte" />
/// <reference types="vite/client" />

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
}

type Gene = string;

interface Publication {
  publicationId: string;
  title: string;
  journal: string;
  abstract: string;
  datasets: Dataset[];
  year?: number[];
  PMID?: string[];
  species?: string[];
  author?: string[];
  disease?: string[];
  cellType?: string[];
}

type Data = Dataset | Gene | Publication;

type Grouping = "Genotype" | "CellType";

interface User {
  id: number;
  name: string;
}

type NicePrimitive =
  | number
  | string
  | boolean
  | null
  | undefined
  | Date
  | NiceObject;
interface NiceObject {
  [k: string]: NicePrimitive | NicePrimitive[];
}

interface BasePageProps<T extends NiceObject = any> {
  user?: User | null;
  __meta?: T;
}

type PageProps<T extends NiceObject = any, M extends NiceObject = any> = T &
  BasePageProps<M>;

declare module "*.svelte" {
  const component: ATypedSvelteComponent;
  export default component;
}

declare module "$meta" {
  interface State {
    route: string;
    path: string;
    params: Record<string, any>;
    query: Record<string, any>;
    user: User;
    extra: any;
  }

  export const __state: State;
  const state: State;
  export default state;
}

interface IndexProps extends PageProps {
  token: string;
}

interface PublicationProps extends PageProps {
  publication: Publication;
}

interface CompareProps extends PageProps {
  genes: string[];
  gene: string;
  order: string[];
}

interface PlotConfig {
  selectedGene: string;
  groupBy: string;
  splitBy: string;
}

interface RenderResult {
  violin: string;
  clustering: string;
}

interface ChartResult {
  clustering: string;
  violin: string;
}

interface Column {
  key: string;
  label: string;
  url?: boolean;
  sortable?: boolean;
}
