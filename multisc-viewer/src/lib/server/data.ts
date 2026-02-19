import { merge } from "lodash-es";
import type {
	DatasetsMap,
	DEGsMap,
	GenesMap,
	GenesRowsMap,
	PublicationsMap,
} from "$lib/types/daemon";
import { getDaemonTargets } from "./daemon";

/**
 * Fetch datasets metadata.
 * @param ds - list of dataset ids to fetch; if `undefined`, fetch all datasets
 * @returns Dictionary mapping dataset ids to metadata
 */
export const getDatasets = async (ds?: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const datasets = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.datasets(dds)),
	);
	return datasets.reduce(merge<DatasetsMap, DatasetsMap>, {});
};

/**
 * Fetch PublicationsMap metadata.
 * @param pub - list of publication ids to fetch; if `undefined`, fetch all PublicationsMap
 * @returns Dictionary mapping publication ids to metadata
 */
export const getPublications = async (pub?: string[]) => {
	const daemonTargets = await getDaemonTargets();
	const PublicationsMap = await Promise.all(
		daemonTargets.entries().map(async ([daemon]) => daemon.publications(pub)),
	);
	return PublicationsMap.reduce(merge<PublicationsMap, PublicationsMap>, {});
};

/**
 * Fetch GenesMap for datasets.
 * @param ds - list of dataset ids to fetch
 * @returns Dictionary mapping dataset ids to GenesMap
 */
export const getGenes = async (ds: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const GenesMap = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.genes(dds)),
	);
	return GenesMap.reduce(merge<GenesMap, GenesMap>, {});
};

/**
 * Fetch differentially expressed GenesMap for datasets.
 * @param ds - list of dataset ids to fetch
 * @returns Dictionary mapping dataset ids to differentially expressed GenesMap
 */
export const getDEGs = async (ds: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const DEGsMap = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.degs(dds)),
	);
	return DEGsMap.reduce(merge<DEGsMap, DEGsMap>, {});
};

/**
 * Fetch rows for gene datasets and differentially expressed GenesMap.
 * @param ds - list of dataset ids to fetch
 * @returns List of gene rows
 * */
export const getGenesRows = async (
	ds: string[],
	limit?: number,
	offset?: number,
) => {
	const daemonTargets = await getDaemonTargets(ds);
	const GenesRowsMap = await Promise.all(
		daemonTargets
			.entries()
			.map(async ([daemon, dds]) => daemon.genesRows(dds, limit, offset)),
	);
	return GenesRowsMap.reduce(merge<GenesRowsMap, GenesRowsMap>, {});
};
