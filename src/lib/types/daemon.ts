export interface StatusResponse {
	datasets: string[];
}

export type UnloadResponse = StatusResponse;

export interface LoadResponse {
	[ds: string]: boolean;
}

export type RenderResponse = string[]; // rendered cache keys
