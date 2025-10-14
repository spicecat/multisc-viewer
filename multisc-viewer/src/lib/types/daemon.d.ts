import type { components, paths } from './api';

export type Datasets = paths['/datasets']['get']['responses']['200']['content']['application/json'];
export type Publications =
	paths['/publications']['get']['responses']['200']['content']['application/json'];
export type Genes = paths['/genes']['get']['responses']['200']['content']['application/json'];
export type DEGs = paths['/degs']['get']['responses']['200']['content']['application/json'];

export type PlotsParams = paths['/plots']['post']['requestBody']['content']['application/json'];
export type Plots = paths['/plots']['post']['responses']['200']['content']['application/json'];
export type Plot = Plots[number];
export type PlotRequests = Record<string, Promise<Plot>>;

export type Dataset = Datasets[string];
export type Publication = Publications[string];
export type Gene = Genes[string];
export type DEG = DEGs[string];
export type Author = components['schemas']['Author'];
export type OntologyTerm = components['schemas']['OntologyTerm'];
