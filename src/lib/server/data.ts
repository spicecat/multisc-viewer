import { datasetsConfig, publicationsConfig } from '$lib/config';
import type { Dataset, Gene, Publication } from '$lib/types/data';
import { existsSync, readFileSync, statSync } from 'node:fs';
import { readFile } from 'node:fs/promises';

type DatasetMeta = Omit<Dataset, 'size'>;

type PublicationMeta = Omit<Publication, 'datasets'> & { datasets: string[] };

const {
	dir: datasetsDir,
	meta: datasetsMeta,
	requiredFiles: datasetsRequiredFiles
} = datasetsConfig;

const { dir: publicationsDir, meta: publicationsMeta } = publicationsConfig;

export const datasets = new Map<string, Dataset>(
	JSON.parse(readFileSync(`${datasetsDir}/${datasetsMeta}`, 'utf-8'))
		.filter(({ id }: DatasetMeta) =>
			datasetsRequiredFiles.every((file: string) => existsSync(`${datasetsDir}/${id}/${file}`))
		)
		.map((ds: DatasetMeta) => [
			ds.id,
			{ ...ds, size: statSync(`${datasetsDir}/${ds.id}/data.rds`).size }
		])
);

export const publications = new Map<string, Publication>(
	JSON.parse(readFileSync(`${publicationsDir}/${publicationsMeta}`, 'utf-8')).map(
		(pub: PublicationMeta) => [
			pub.id,
			{ ...pub, datasets: pub.datasets.map((id: string) => datasets.get(id)).filter(Boolean) }
		]
	)
);

export const getGenes = async (datasets: string[]): Promise<Gene[]> => {
	const genes = await Promise.all(
		datasets.map(async (id) =>
			JSON.parse(await readFile(`${datasetsDir}/${id}/genes.json`, 'utf-8'))
		)
	);
	return Array.from(new Set(genes.flat())).toSorted();
};
