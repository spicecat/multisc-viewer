import { Injectable } from '@nestjs/common';
import { spawn } from 'child_process';
import { randomBytes } from 'crypto';
import { existsSync, readdirSync, readFileSync, rmSync } from 'fs';

export interface ChartResult {
	clustering: string;
	violin: string;
}

export interface Dataset {
	name: string;
	year: number;
	region: string[];
	PMID: string;
	species: string;
	author: string;
	disease: string[];
	size: number;
}

@Injectable()
export class AppService {
	public static readonly META: Dataset[] = JSON.parse(readFileSync('datasets/meta.json').toString());

	private readonly cache: Map<string, Promise<ChartResult>>;
	private readonly expirations: Map<string, NodeJS.Timeout>;

	public constructor() {
		this.cache = new Map();
		this.expirations = new Map();
	}

	// TODO: probably don't need all the checks given the go script
	public getDatasets(): Dataset[] {
		return readdirSync('datasets')
			.filter((dataset) => existsSync(`datasets/${dataset}/data.rds`) && existsSync(`datasets/${dataset}/genes.json`))
			.map((name) => AppService.META.find((dataset) => dataset.name === name))
			.filter((dataset) => !!dataset);
	}

	public getGenes(datasets: string[]): string[] {
		const geneSets = datasets.map<string[]>((dataset) => JSON.parse(readFileSync(`datasets/${dataset}/genes.json`).toString()));

		return geneSets.slice(1).reduce((curr, next) => curr.filter((gene) => next.includes(gene)), geneSets[0]);
	}

	public async render(dataset: string, gene: string, groupBy: string, splitBy: string): Promise<ChartResult> {
		const key = this._cacheKey(dataset, gene, groupBy, splitBy);

		if (this.cache.has(key)) {
			if (this.expirations.has(key)) {
				clearTimeout(this.expirations.get(key)!);
				this.expirations.set(
					key,
					setTimeout(() => {
						this.cache.delete(key);
						this.expirations.delete(key);
					}, 30 * 60 * 1000)
				);
			}

			return this.cache.get(key)!;
		} else {
			const result = new Promise<ChartResult>((resolve, reject) => {
				const clusteringFile = randomBytes(16).toString('hex') + '.png';
				const violinFile = randomBytes(16).toString('hex') + '.png';

				const proc = spawn('Rscript', ['../../plot.r', gene, groupBy, splitBy, clusteringFile, violinFile], { cwd: `datasets/${dataset}` });

				let stdout = '',
					stderr = '';

				proc.stdout.on('data', (chunk) => (stdout += chunk));
				proc.stderr.on('data', (chunk) => (stderr += chunk));

				proc.on('error', () => reject(new Error(`Error reading genes of '${dataset}'`)));
				proc.on('exit', () => {
					try {
						const clustering = 'data:image/png;base64,' + readFileSync(`datasets/${dataset}/${clusteringFile}`).toString('base64');
						const violin = 'data:image/png;base64,' + readFileSync(`datasets/${dataset}/${violinFile}`).toString('base64');

						rmSync(`datasets/${dataset}/${clusteringFile}`);
						rmSync(`datasets/${dataset}/${violinFile}`);

						resolve({ clustering, violin });
					} catch (e) {
						console.error(e);
						console.error('stdout', stdout);
						console.error('stderr', stderr);
						reject(new Error(`Unable to find plot of '${dataset}'`));
					}
				});
			});

			result.then(() =>
				this.expirations.set(
					key,
					setTimeout(() => {
						this.cache.delete(key);
						this.expirations.delete(key);
					}, 30 * 60 * 1000)
				)
			);

			this.cache.set(key, result);
			return result;
		}
	}

	public preload(datasets: string[], gene: string, groupBy: string, splitBy: string): void {
		datasets.forEach((dataset) => this.render(dataset, gene, groupBy, splitBy));
	}

	public async generate(datasets: string[], gene: string, groupBy: string, splitBy: string): Promise<Record<string, ChartResult>> {
		return Promise.all(datasets.map((dataset) => this.render(dataset, gene, groupBy, splitBy))).then((results) =>
			Object.fromEntries(results.map((result, i) => [datasets[i], result]))
		);
	}

	private _cacheKey(dataset: string, gene: string, groupBy: string, splitBy: string): string {
		return `${dataset}:${gene}:${groupBy}:${splitBy}`;
	}
}

