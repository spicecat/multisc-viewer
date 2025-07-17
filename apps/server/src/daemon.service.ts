import { Injectable, OnModuleInit } from "@nestjs/common";
import { ConfigService } from "@nestjs/config";
import axios from "axios";
import NodeCache from "node-cache";
import { DataService } from "./data.service";

interface Daemon {
  port: number;
  load: number;
  datasets: Set<string>;
}

@Injectable()
export class DaemonService implements OnModuleInit {
  private readonly DAEMON_SERVER: string;
  private readonly datasetCache: NodeCache;
  private readonly procs: Map<number, Daemon> = new Map();

  public constructor(
    private readonly configService: ConfigService,
    private readonly dataService: DataService
  ) {
    this.DAEMON_SERVER = this.configService.get<string>("daemon.server")!;
    this.datasetCache = new NodeCache({
      stdTTL: this.configService.get<number>("daemon.ttl"),
    });

    this.datasetCache.on("expired", (ds: string, daemon: Daemon) =>
      this._unload(daemon, ds)
    );
  }

  async onModuleInit(): Promise<void> {
    console.log("Initializing DaemonService...");
    await new Promise((resolve) => setTimeout(resolve, 2000)); // Wait for daemons to start
    await Promise.all(
      this.configService.get<number[]>("daemon.ports")!.map(async (port) => {
        try {
          const response = await axios.get(
            `${this.DAEMON_SERVER}:${port}/status`
          );
          if (response.data)
            this.procs.set(port, {
              port,
              load: 0,
              datasets: new Set<string>(),
            });
        } catch (error) {
          console.error(`Daemon on port ${port} not responding`);
        }
      })
    );
  }

  private async _unload(daemon: Daemon, ds: string): Promise<void> {
    try {
      await axios.post(`${this.DAEMON_SERVER}:${daemon.port}/unload`, { ds });
      daemon.datasets.delete(ds);
      const datasetInfo = this.dataService.datasets.get(ds);
      if (datasetInfo) {
        daemon.load -= datasetInfo.size;
      }
    } catch (error: any) {
      console.error(
        `Failed to unload dataset ${ds} from port ${daemon.port}:`,
        error.message
      );
    }
  }

  public async load(datasets: string | string[]): Promise<void> {
    if (this.procs.size === 0)
      throw new Error("No daemons available to load dataset.");

    if (!Array.isArray(datasets)) datasets = [datasets];

    await Promise.all(
      datasets.map(async (ds) => {
        if (this.datasetCache.has(ds)) this.datasetCache.ttl(ds);
        else {
          const datasetInfo = this.dataService.datasets.get(ds);
          if (!datasetInfo) {
            console.error(`Dataset ${ds} not found.`);
            return;
          }

          const daemon = [...this.procs.values()].reduce((a, b) =>
            a.load < b.load ? a : b
          ); // daemon with the least load
          daemon.load += datasetInfo.size;

          try {
            await axios.post(`${this.DAEMON_SERVER}:${daemon.port}/load`, {
              ds,
            });
            daemon.datasets.add(ds);
            this.datasetCache.set(ds, daemon);
          } catch (error: any) {
            daemon.load -= datasetInfo.size;
            console.error(
              `Failed to load dataset ${ds} on port ${daemon.port}:`,
              error.message
            );
          }
        }
      })
    );

    console.log("Datasets in cache:", this.datasetCache.keys());
  }

  public async render(
    ds: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): Promise<string> {
    await this.load(ds);
    const daemon = this.datasetCache.get<Daemon>(ds);
    if (!daemon)
      throw new Error(`Dataset ${ds} is not loaded`);

    try {
      const response = await axios.post(
        `${this.DAEMON_SERVER}:${daemon.port}/render`,
        {
          ds,
          gene,
          groupBy,
          splitBy,
        }
      );
      return response.data.opid;
    } catch (error: any) {
      console.error(
        `Failed to render dataset ${ds} on port ${daemon.port}:`,
        error.message
      );
      throw new Error(`Failed to render plot for dataset ${ds}`);
    }
  }
}
