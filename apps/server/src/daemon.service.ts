import { Injectable, OnModuleInit } from "@nestjs/common";
import axios from "axios";
import NodeCache from "node-cache";

const DAEMON_SERVER = "http://localhost";
const DAEMON_PORT_START = 8001;
const DAEMON_COUNT = 10;

interface Daemon {
  port: number;
  load: number;
  datasets: Set<string>;
}

@Injectable()
export class DaemonService implements OnModuleInit {
  readonly #datasetCache = new NodeCache({ stdTTL: 30 * 60 }); // 30 minutes TTL
  readonly #procs: Map<number, Daemon> = new Map();

  public constructor() {
    this.#datasetCache.on("expired", (ds: string, daemon: Daemon) => this._unload(daemon, ds));
  }

  async onModuleInit(): Promise<void> {
    console.log("Initializing DaemonService...");
    await new Promise((resolve) => setTimeout(resolve, 2000)); // Wait for daemons to start

    await Promise.all(
      Array.from({ length: DAEMON_COUNT }, async (_, i) => {
        const port = DAEMON_PORT_START + i;
        try {
          const response = await axios.get(`${DAEMON_SERVER}:${port}/status`);
          if (response.data)
            this.#procs.set(port, {
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

  public async load(datasets: string | string[]): Promise<void> {
    if (!Array.isArray(datasets)) datasets = [datasets];

    await Promise.all(
      datasets.map(async (ds, index) => {
        // if (this.#datasetCache.has(ds)) this.#datasetCache.ttl(ds);
        // else {
        //   const daemon = this.#procs[index % this.#procs.length];
        //   try {
        //     await axios.post(`${DAEMON_SERVER}:${daemon.port}/load`, { ds });
        //     daemon.datasets.add(ds);
        //     this.#datasetCache.set(ds, daemon);
        //   } catch (error: any) {
        //     console.error(
        //       `Failed to load dataset ${ds} on port ${daemon.port}:`,
        //       error.message
        //     );
        //     // Handle error, maybe try another daemon or throw
        //   }
        // }
      })
    );

    console.log("Datasets in cache:", this.#datasetCache.keys());
  }

  public async render(
    ds: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): Promise<string> {
    await this.load(ds);
    const daemon = this.#datasetCache.get<Daemon>(ds);
    if (!daemon) {
      throw new Error(`Dataset ${ds} is not loaded`);
    }

    try {
      const response = await axios.post(
        `${DAEMON_SERVER}:${daemon.port}/render`,
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

  private async _unload(daemon: Daemon, ds: string): Promise<void> {
    try {
      await axios.post(`${DAEMON_SERVER}:${daemon.port}/unload`, { ds });
      daemon.datasets.delete(ds);
    } catch (error: any) {
      console.error(
        `Failed to unload dataset ${ds} from port ${daemon.port}:`,
        error.message
      );
    }
  }
}
