import { Injectable } from "@nestjs/common";
import { existsSync, readdirSync, readFileSync, rmSync } from "fs";
import { publications } from "../src/config/publications";
import { DaemonService } from "./daemon.service";

@Injectable()
export class AppService {
  public static readonly META: Dataset[] = JSON.parse(
    readFileSync("datasets/meta.json").toString(),
  );

  private readonly cache: Map<string, Promise<ChartResult>>;
  private readonly expirations: Map<string, NodeJS.Timeout>;
  private readonly activity: Map<string, NodeJS.Timeout>;

  public constructor(private readonly daemon: DaemonService) {
    this.cache = new Map();
    this.expirations = new Map();
    this.activity = new Map();
  }

  public getPublications(): Publication[] {
    return publications.map((publication) => ({
      ...publication,
      datasets: this.getDatasets(), // TODO: add datasets to publication
    }));
  }

  public getDatasets(): Dataset[] {
    return readdirSync("datasets")
      .filter(
        (ds) =>
          existsSync(`datasets/${ds}/data.rds`) &&
          existsSync(`datasets/${ds}/genes.json`),
      )
      .map((name) => AppService.META.find((ds) => ds.name === name))
      .filter((ds) => !!ds);
  }

  public getGenes(datasets: string[]): string[] {
    const geneSets = datasets.map<string[]>((ds) =>
      JSON.parse(readFileSync(`datasets/${ds}/genes.json`).toString()),
    );
    return geneSets
      .slice(1)
      .reduce(
        (curr, next) => curr.filter((gene) => next.includes(gene)),
        geneSets[0],
      );
  }

  public async render(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string,
    token: string | null,
  ): Promise<ChartResult> {
    const key = this._cacheKey(dataset, gene, groupBy, splitBy);

    if (token) this._refresh(token);

    if (this.cache.has(key)) {
      if (this.expirations.has(key)) this.expirations.get(key)!.refresh();
      return this.cache.get(key)!;
    } else {
      const result = this.daemon
        .render(dataset, gene, groupBy, splitBy)
        .then((opid) => {
          try {
            const clustering =
              "data:image/png;base64," +
              readFileSync(`datasets/${dataset}/${opid}_umap.png`).toString(
                "base64",
              );
            const violin =
              "data:image/png;base64," +
              readFileSync(`datasets/${dataset}/${opid}_vln.png`).toString(
                "base64",
              );

            rmSync(`datasets/${dataset}/${opid}_umap.png`);
            rmSync(`datasets/${dataset}/${opid}_vln.png`);

            return { clustering, violin };
          } catch (e) {
            console.error(e);
            throw new Error(`Unable to find plot of '${dataset}'`);
          }
        });

      result
        .then(() =>
          this.expirations.set(
            key,
            setTimeout(
              () => {
                this.cache.delete(key);
                this.expirations.delete(key);
              },
              30 * 60 * 1000,
            ),
          ),
        )
        .then(() => {
          if (token) this._refresh(token);
        });

      this.cache.set(key, result);
      return result;
    }
  }

  public async preload(
    token: string | null,
    datasets: string[],
  ): Promise<void> {
    if (token !== null) {
      this._refresh(token);
    }

    return this.daemon.batch(token, datasets).then(() => {
      if (token) this._refresh(token);
    });
  }

  private _cacheKey(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string,
  ): string {
    return `${dataset}:${gene}:${groupBy}:${splitBy}`;
  }

  private _refresh(token: string) {
    if (this.activity.has(token)) this.activity.get(token)!.refresh();
    else {
      this.activity.set(
        token,
        setTimeout(() => {
          this.activity.delete(token);
          this.daemon.unload(token);
        }, 60_000),
      );
    }
  }
}
