import { existsSync, readFileSync } from "node:fs";
import type { Publication, Dataset } from "$lib/types";
import { CONFIG } from "$lib/config";

export class DataService {
  publications: Map<number, Publication>;
  datasets: Map<string, Dataset>;
  datasetsBasePath: string;

  constructor(config: typeof CONFIG) {
    const publicationsMeta = `${config.publications.basePath}/${config.publications.metaFile}`;
    this.publications = new Map(
      JSON.parse(readFileSync(publicationsMeta, "utf-8")).map(
        (pub: Publication) => [pub.PMID, pub]
      )
    );
    this.datasetsBasePath = config.datasets.basePath;
    const datasetsMeta = `${this.datasetsBasePath}/${config.datasets.metaFile}`;
    this.datasets = new Map(
      JSON.parse(readFileSync(datasetsMeta, "utf-8"))
        .filter(({ name }: Dataset) =>
          config.datasets.requiredFiles.every((file: string) =>
            existsSync(`${this.datasetsBasePath}/${name}/${file}`)
          )
        )
        .map((ds: Dataset) => [ds.name, ds])
    );
  }

  getGenes(datasets: string[]): string[] {
    return Array.from(
      new Set(
        datasets
          .map((ds) =>
            JSON.parse(
              readFileSync(
                `${this.datasetsBasePath}/${ds}/genes.json`,
                "utf-8"
              )
            )
          )
          .flat()
      )
    ).toSorted();
  }
}

export const dataService = new DataService(CONFIG);
