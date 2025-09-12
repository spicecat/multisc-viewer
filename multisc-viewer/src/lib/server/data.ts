import type { Dataset, DEGs, Gene, Publication } from '$lib/types/data';
import { uniq } from 'lodash-es';
import { existsSync, readFileSync, statSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { datasetsConfig, publicationsConfig } from './config';

type DatasetMeta = Omit<Dataset, 'size'>;

type PublicationMeta = Omit<Publication, 'datasets'> & { datasets: string[] };

const {
	dir: datasetsDir,
	meta: datasetsMeta,
	requiredFiles: datasetsRequiredFiles
} = datasetsConfig;

const { dir: publicationsDir, meta: publicationsMeta } = publicationsConfig;

export const datasets: Record<string, Dataset> = Object.fromEntries(
	existsSync(`${datasetsDir}/${datasetsMeta}`)
		? JSON.parse(readFileSync(`${datasetsDir}/${datasetsMeta}`, 'utf-8'))
				.filter(({ id }: DatasetMeta) =>
					datasetsRequiredFiles.every((file) => existsSync(`${datasetsDir}/${id}/${file}`))
				)
				.map((ds: DatasetMeta) => [
					ds.id,
					{ ...ds, size: statSync(`${datasetsDir}/${ds.id}/data.rds`).size }
				])
		: []
);

export const publications: Record<string, Publication> = Object.fromEntries(
	existsSync(`${publicationsDir}/${publicationsMeta}`)
		? JSON.parse(readFileSync(`${publicationsDir}/${publicationsMeta}`, 'utf-8')).map(
				(pub: PublicationMeta) => [
					pub.id,
					{
						...pub,
						datasets: pub.datasets.filter((ds) => ds in datasets).map((ds) => datasets[ds])
					}
				]
			)
		: []
);

export const getGenes = async (datasets: string[]): Promise<Gene[]> => {
	const genes: Gene[][] = await Promise.all(
		datasets.map(async (ds) =>
			JSON.parse(await readFile(`${datasetsDir}/${ds}/genes.json`, 'utf-8'))
		)
	);
	return uniq(genes.flat()).toSorted();
};

export const getDEGs = async (datasets: string[]): Promise<DEGs> => {
	const degs: DEGs[] = await Promise.all(
		datasets.map(async (ds) =>
			existsSync(`${datasetsDir}/${ds}/DEGs.json`)
				? JSON.parse(await readFile(`${datasetsDir}/${ds}/DEGs.json`, 'utf-8'))
				: {}
		)
	);
	return Object.assign({}, ...degs);
};
