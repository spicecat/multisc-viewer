import axios from "axios";
import NodeCache from "node-cache";
import { dataService } from "./data.service";
import { CONFIG } from "$lib/config";

interface Daemon {
  port: number;
  load: number;
  datasets: Set<string>;
}

class DaemonService {
  DAEMON_SERVER: string;
  datasetCache: NodeCache;
  procs: Map<number, Daemon> = new Map();

  constructor(config: typeof CONFIG) {
    this.DAEMON_SERVER = config.daemon.server;
    this.datasetCache = new NodeCache({ stdTTL: config.daemon.ttl });
    this.datasetCache.on("expired", (ds: string, daemon: Daemon) => this._unload(daemon, ds));
    this.init(config.daemon.ports);
  }

  async init(ports: number[]) {
    await Promise.all(
      ports.map(async (port) => {
        try {
          const response = await axios.get(`${this.DAEMON_SERVER}:${port}/status`);
          if (response.data)
            this.procs.set(port, { port, load: 0, datasets: new Set<string>() });
        } catch (error) {
          // Daemon not responding
        }
      })
    );
  }

  async _unload(daemon: Daemon, ds: string) {
    try {
      await axios.post(`${this.DAEMON_SERVER}:${daemon.port}/unload`, { ds });
      daemon.datasets.delete(ds);
      const datasetInfo = dataService.datasets.get(ds);
      if (datasetInfo) daemon.load -= datasetInfo.size;
    } catch (error) {}
  }

  async load(datasets: string | string[]) {
    if (this.procs.size === 0) throw new Error("No daemons available to load dataset.");
    if (!Array.isArray(datasets)) datasets = [datasets];
    await Promise.all(
      datasets.map(async (ds) => {
        if (this.datasetCache.has(ds)) this.datasetCache.ttl(ds);
        else {
          const datasetInfo = dataService.datasets.get(ds);
          if (!datasetInfo) return;
          const daemon = [...this.procs.values()].reduce((a, b) => (a.load < b.load ? a : b));
          daemon.load += datasetInfo.size;
          try {
            await axios.post(`${this.DAEMON_SERVER}:${daemon.port}/load`, { ds });
            daemon.datasets.add(ds);
            this.datasetCache.set(ds, daemon);
          } catch (error) {
            daemon.load -= datasetInfo.size;
          }
        }
      })
    );
  }

  async render(ds: string, gene: string, groupBy: string, splitBy: string): Promise<string> {
    await this.load(ds);
    const daemon = this.datasetCache.get<Daemon>(ds);
    if (!daemon) throw new Error(`Dataset ${ds} is not loaded`);
    const response = await axios.post(`${this.DAEMON_SERVER}:${daemon.port}/render`, { ds, gene, groupBy, splitBy });
    return response.data.opid;
  }
}

export const daemonService = new DaemonService(CONFIG);
