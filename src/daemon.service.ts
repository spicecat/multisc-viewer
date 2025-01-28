import { Injectable } from "@nestjs/common";
import { AppService } from "./app.service";
import { Daemon } from "./utils/Daemon";
import { approximateRenderTime } from "./utils/utils";

@Injectable()
export class DaemonService {
  private readonly procs: Daemon[];

  public constructor() {
    this.procs = new Array(10).fill(null).map(() => new Daemon());
  }

  public async batch(loader: string | null, datasets: string[]): Promise<void> {
    return Promise.all(datasets.map((ds) => this.load(loader, ds))).then();
  }

  public async load(loader: string | null, ds: string): Promise<void> {
    const proc = this.procs.find((proc) =>
      proc.datasets.some((dataset) => dataset.ds === ds),
    );

    if (!proc) {
      if (loader !== null) {
        const availableProc = this.procs.find(
          (proc) => !proc.datasets.some((set) => set.loaders.has(loader)),
        );

        if (!availableProc) {
          const newProc = this._spawn();

          return newProc.load(loader, ds);
        } else {
          return availableProc.load(loader, ds);
        }
      } else {
        const bestProc = this.procs.reduce(
          (best, proc) =>
            proc.datasets.length < best.datasets.length ? proc : best,
          this.procs[0],
        );

        return bestProc.load(loader, ds).then(() => {
          setTimeout(
            () => {
              if (
                bestProc.datasets.find((dataset) => dataset.ds === ds)?.loaders
                  .size === 0
              ) {
                if (bestProc.datasets.length === 1 && this.procs.length > 10) {
                  bestProc.kill();
                  this.procs.splice(this.procs.indexOf(bestProc), 1);
                } else {
                  bestProc.unload(ds);
                }
              }
            },
            approximateRenderTime(
              AppService.META.find((dataset) => dataset.name === ds)!.size,
            ) * 30,
          );
        });
      }
    } else {
      if (loader !== null) {
        proc.datasets.find((dataset) => dataset.ds === ds)!.loaders.add(loader);
      }

      // TODO: not precise (could be an operation later than the actual end of the load of this dataset)
      // consider changing to Promise.resolve if all operations dependent on this case are queued operations
      return proc.idle;
    }
  }

  public async unload(loader: string): Promise<void> {
    const promises: Promise<void>[] = [];

    this.procs.forEach((proc) => {
      proc.datasets.forEach((ds) => {
        if (ds.loaders.has(loader)) {
          if (ds.loaders.size === 1) {
            if (proc.datasets.length === 1 && this.procs.length > 10) {
              proc.kill();
              this.procs.splice(this.procs.indexOf(proc), 1);
            } else {
              promises.push(proc.unload(ds.ds));

              proc.datasets.splice(proc.datasets.indexOf(ds), 1);
            }
          } else {
            ds.loaders.delete(loader);
          }
        }
      });
    });

    return Promise.all(promises).then();
  }

  public async render(
    ds: string,
    gene: string,
    groupBy: string,
    splitBy: string,
  ): Promise<string> {
    const proc = this.procs.find((proc) =>
      proc.datasets.some((dataset) => dataset.ds === ds),
    );

    if (proc) {
      return proc.render(ds, gene, groupBy, splitBy);
    } else {
      return this.load(null, ds).then(() =>
        this.render(ds, gene, groupBy, splitBy),
      );
    }
  }

  private _spawn(): Daemon {
    const proc = new Daemon();

    this.procs.push(proc);

    return proc;
  }
}
