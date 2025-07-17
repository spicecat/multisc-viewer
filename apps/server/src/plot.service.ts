import { Injectable } from "@nestjs/common";
import { ConfigService } from "@nestjs/config";
import NodeCache from "node-cache";
import { readFileSync, rmSync } from "node:fs";
import { DaemonService } from "./daemon.service";

@Injectable()
export class PlotService {
  private readonly datasetsBasePath: string;
  private readonly plotCache: NodeCache;

  public constructor(
    private readonly configService: ConfigService,
    private readonly daemon: DaemonService
  ) {
    this.datasetsBasePath =
      this.configService.get<string>("datasets.basePath")!;

    this.plotCache = new NodeCache({
      stdTTL: this.configService.get<number>("plot.ttl"),
    });
  }

  private _cacheKey(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): string {
    return `${dataset}:${gene}:${groupBy}:${splitBy}`;
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

    const opid = await this.daemon.render(dataset, gene, groupBy, splitBy);
    console.log(`Rendering plots for ${dataset} with opid ${opid}`);
    const clusteringPath = `${this.datasetsBasePath}/${dataset}/${opid}_umap.png`;
    const violinPath = `${this.datasetsBasePath}/${dataset}/${opid}_vln.png`;
    const featurePath = `${this.datasetsBasePath}/${dataset}/${opid}_feat.png`;

    const clustering =
      "data:image/png;base64," +
      readFileSync(clusteringPath).toString("base64");
    const violin =
      "data:image/png;base64," + readFileSync(violinPath).toString("base64");
    const feature =
      "data:image/png;base64," + readFileSync(featurePath).toString("base64");

    rmSync(clusteringPath);
    rmSync(violinPath);
    rmSync(featurePath);

    const result: ChartResult = {
      clustering,
      violin,
      feature,
    };
    this.plotCache.set(key, result);
    return result;
  }
}
