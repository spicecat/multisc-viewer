import NodeCache from "node-cache";
import { readFileSync, rmSync } from "node:fs";
import { daemonService } from "./daemon.service";
import { CONFIG } from "$lib/config";
import type { ChartResult } from "$lib/types";

class PlotService {
  datasetsBasePath: string;
  plotCache: NodeCache;

  constructor(config: typeof CONFIG) {
    this.datasetsBasePath = config.datasets.basePath;
    this.plotCache = new NodeCache({ stdTTL: config.plot.ttl });
  }

  _cacheKey(dataset: string, gene: string, groupBy: string, splitBy: string): string {
    return `${dataset}:${gene}:${groupBy}:${splitBy}`;
  }

  async render(dataset: string, gene: string, groupBy: string, splitBy: string): Promise<ChartResult> {
    const key = this._cacheKey(dataset, gene, groupBy, splitBy);
    const cachedResult = this.plotCache.get<ChartResult>(key);
    if (cachedResult) return cachedResult;
    const opid = await daemonService.render(dataset, gene, groupBy, splitBy);
    const clusteringPath = `${this.datasetsBasePath}/${dataset}/${opid}_umap.png`;
    const violinPath = `${this.datasetsBasePath}/${dataset}/${opid}_vln.png`;
    const featurePath = `${this.datasetsBasePath}/${dataset}/${opid}_feat.png`;
    const clustering = "data:image/png;base64," + readFileSync(clusteringPath, "base64");
    const violin = "data:image/png;base64," + readFileSync(violinPath, "base64");
    const feature = "data:image/png;base64," + readFileSync(featurePath, "base64");
    rmSync(clusteringPath);
    rmSync(violinPath);
    rmSync(featurePath);
    const result: ChartResult = { clustering, violin, feature };
    this.plotCache.set(key, result);
    return result;
  }
}

export const plotService = new PlotService(CONFIG);
