import type { paths } from "./api";

export type DatasetsIndex =
	paths["/datasets-index"]["get"]["responses"]["200"]["content"]["application/json"];

type Datasets =
	paths["/datasets"]["get"]["responses"]["200"]["content"]["application/json"];
export type Dataset = Datasets[number];
export type DatasetsMap = Record<string, Dataset>;

type Publications =
	paths["/publications"]["get"]["responses"]["200"]["content"]["application/json"];
export type Publication = Publications[number];
export type PublicationsMap = Record<string, Publication>;

type Genes =
	paths["/genes"]["get"]["responses"]["200"]["content"]["application/json"];
export type Gene = Genes[number];
export type GenesMap = Record<string, Gene>;

type DEGs =
	paths["/degs"]["get"]["responses"]["200"]["content"]["application/json"];
export type DEG = DEGs[number];
export type DEGsMap = Record<string, DEG>;

type GenesRows =
	paths["/genes-rows"]["get"]["responses"]["200"]["content"]["application/json"];
export type GenesRow = GenesRows[number];
export type GenesRowsMap = Record<string, GenesRow>;

export type PlotsParams = paths["/plots"]["post"]["parameters"]["query"];
export type Plots =
	paths["/plots"]["post"]["responses"]["200"]["content"]["application/json"];
type Plot = Plots[number];
export type PlotRequests = Record<string, Promise<Plot>>;

export type Author = Dataset["author"];
export type OntologyTerm = Dataset["healthCondition"];
