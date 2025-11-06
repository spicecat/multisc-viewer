import type { Datasets, DEGs, Genes, GenesRows, Publications } from '$lib/types/daemon';
import { merge } from 'lodash-es';
import { getDaemonTargets } from './daemon';

/**
 * Fetch datasets metadata.
 * @param ds - list of dataset ids to fetch; if `undefined`, fetch all datasets
 * @returns Dictionary mapping dataset ids to metadata
 */
export const getDatasets = async (ds?: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const datasets = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.datasets(dds))
	);
	return datasets.reduce(merge<Datasets, Datasets>, {});
};

/**
 * Fetch publications metadata.
 * @param pub - list of publication ids to fetch; if `undefined`, fetch all publications
 * @returns Dictionary mapping publication ids to metadata
 */
export const getPublications = async (pub?: string[]) => {
	const daemonTargets = await getDaemonTargets();
	const publications = await Promise.all(
		daemonTargets.entries().map(async ([daemon]) => daemon.publications(pub))
	);
	return publications.reduce(merge<Publications, Publications>, {});
};

/**
 * Fetch genes for datasets.
 * @param ds - list of dataset ids to fetch
 * @returns Dictionary mapping dataset ids to genes
 */
export const getGenes = async (ds: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const genes = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.genes(dds))
	);
	return genes.reduce(merge<Genes, Genes>, {});
};

/**
 * Fetch differentially expressed genes for datasets.
 * @param ds - list of dataset ids to fetch
 * @returns Dictionary mapping dataset ids to differentially expressed genes
 */
export const getDEGs = async (ds: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const degs = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.degs(dds))
	);
	return degs.reduce(merge<DEGs, DEGs>, {});
};

/**
 * Fetch rows for gene datasets and differentially expressed genes.
 * @param ds - list of dataset ids to fetch
 * @returns List of gene rows
 * */
export const getGenesRows = async (ds: string[]) => {
	const daemonTargets = await getDaemonTargets(ds);
	const genesRows = await Promise.all(
		daemonTargets.entries().map(async ([daemon, dds]) => daemon.genesRows(dds))
	);
	return genesRows.reduce(merge<GenesRows, GenesRows>, []);
};
