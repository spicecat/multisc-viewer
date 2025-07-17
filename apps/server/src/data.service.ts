import { ConfigService } from "@nestjs/config";
import { Injectable } from "@nestjs/common";
import { existsSync, readFileSync } from "node:fs";

@Injectable()
export class DataService {
  readonly publications: Map<number, Publication>;
  readonly datasets: Map<string, Dataset>;
  readonly datasetsBasePath: string;

  constructor(private readonly configService: ConfigService) {
    const publicationsMeta = `${this.configService.get<string>("publications.basePath")!}/${this.configService.get<string>("publications.metaFile")!}`;
    this.publications = new Map(
      JSON.parse(readFileSync(publicationsMeta).toString()).map(
        (pub: Publication) => [pub.PMID, pub]
      )
    );

    this.datasetsBasePath =
      this.configService.get<string>("datasets.basePath")!;
    const datasetsMeta = `${this.datasetsBasePath}/${this.configService.get<string>("datasets.metaFile")!}`;
    this.datasets = new Map(
      JSON.parse(readFileSync(datasetsMeta).toString())
        .filter(({ name }: Dataset) =>
          this.configService
            .get<string[]>("datasets.requiredFiles")!
            .every((file) =>
              existsSync(`${this.datasetsBasePath}/${name}/${file}`)
            )
        )
        .map((ds: Dataset) => [ds.name, ds])
    );
  }

  public getGenes(datasets: string[]): string[] {
    return Array.from(
      new Set(
        datasets
          .map((ds) =>
            JSON.parse(
              readFileSync(
                `${this.datasetsBasePath}/${ds}/genes.json`
              ).toString()
            )
          )
          .flat()
      )
    ).toSorted();
  }
}
