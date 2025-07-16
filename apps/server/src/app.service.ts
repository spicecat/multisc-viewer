import { Injectable } from "@nestjs/common";
import { ConfigService } from "@nestjs/config";
import NodeCache from "node-cache";
import { existsSync, readFileSync, rmSync } from "node:fs";
import { DaemonService } from "./daemon.service";

export interface ChartResult {
  clustering: string;
  violin: string;
  feature: string;
}

@Injectable()
export class AppService {
  private readonly plotCache: NodeCache;
  private readonly datasetsBasePath: string;
  readonly publications: Map<number, Publication>;
  readonly datasets: Map<string, Dataset>;

  public constructor(
    private readonly configService: ConfigService,
    private readonly daemon: DaemonService
  ) {
    this.plotCache = new NodeCache({
      stdTTL: this.configService.get<number>("plot.ttl"),
    });

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

  public async render(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): Promise<ChartResult> {
    const key = this._cacheKey(dataset, gene, groupBy, splitBy);
    const cachedResult = this.plotCache.get<Promise<ChartResult>>(key);
    if (cachedResult) return cachedResult;

    const result = this.daemon
      .render(dataset, gene, groupBy, splitBy)
      .then((opid) => {
        try {
          const clusteringPath = `${this.datasetsBasePath}/${dataset}/${opid}_umap.png`;
          const violinPath = `${this.datasetsBasePath}/${dataset}/${opid}_vln.png`;
          const featurePath = `${this.datasetsBasePath}/${dataset}/${opid}_feat.png`;

          const clustering =
            "data:image/png;base64," +
            readFileSync(clusteringPath).toString("base64");
          const violin =
            "data:image/png;base64," +
            readFileSync(violinPath).toString("base64");
          const feature =
            "data:image/png;base64," +
            readFileSync(featurePath).toString("base64");

          rmSync(clusteringPath);
          rmSync(violinPath);
          rmSync(featurePath);

          return { clustering, violin, feature } as ChartResult;
        } catch (e) {
          console.error(e);
          this.plotCache.del(key);
          throw new Error(`Unable to find plot for '${dataset}'`);
        }
      });

    this.plotCache.set(key, result);
    return result;
  }

  public async load(datasets: string[]): Promise<void> {
    return this.daemon.load(datasets);
  }

  private _cacheKey(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): string {
    return `${dataset}:${gene}:${groupBy}:${splitBy}`;
  }
}
