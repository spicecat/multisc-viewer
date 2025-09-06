import type { Dataset, Publication } from './data';

export interface DatasetsResponse {
	datasets: Dataset[];
}

export interface PublicationsResponse {
	publications: Publication[];
}

export interface StatusResponse {
	datasets: string[];
}

export type UnloadResponse = StatusResponse;

export interface LoadResponse {
	[ds: string]: boolean;
}

export type RenderResponse = string[]; // rendered cache keys
