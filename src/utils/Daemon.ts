import type { ChildProcessWithoutNullStreams } from "node:child_process";
import { spawn } from "node:child_process";
import { randomBytes } from "node:crypto";

interface DSRecord {
  ds: string;
  loaders: Set<string>;
}

interface ProcRecord {
  proc: ChildProcessWithoutNullStreams;
  acks: ((opid: string) => boolean)[];
  datasets: DSRecord[];
  idle: Promise<void>;
}

export class Daemon {
  private readonly proc: ChildProcessWithoutNullStreams;

  private _acks: ((opid: string) => boolean)[];
  private _datasets: DSRecord[];
  private _idle: Promise<void>;

  public constructor() {
    this.proc = spawn("Rscript", ["daemon.r"]);

    this._acks = [];
    this._datasets = [];
    this._idle = Promise.resolve();

    let data = "";
    this.proc.stdout.on("data", (chunk) => {
      data += chunk;

      const match = /^ack ([^\n]+)\n/.exec(data);
      if (match !== null) {
        const opid = match[1];
        this._acks = this._acks.filter((ack) => ack(opid));
        data = data.slice(match[0].length);
      }
    });
  }

  public get idle(): Promise<void> {
    return this._idle;
  }

  public get datasets(): DSRecord[] {
    return this._datasets;
  }

  public exec(
    executor: (proc: ChildProcessWithoutNullStreams, opid: string) => void
  ): Promise<void> {
    const opid = randomBytes(16).toString("hex");

    const promise = new Promise<void>((resolve, reject) => {
      this._idle.then(() => {
        executor(this.proc, opid);

        this._acks.push((id) => {
          if (id === opid) {
            resolve();
            return true;
          } else return false;
        });
      });
    });

    this._idle = promise;

    return promise;
  }

  public async load(loader: string | null, ds: string): Promise<void> {
    const dataset = {
      ds,
      loaders: loader !== null ? new Set([loader]) : new Set<string>(),
    };

    this._datasets.push(dataset);

    return this.exec((proc, opid) => proc.stdin.write(`${opid} load ${ds}\n`));
  }

  public async render(
    ds: string,
    gene: string,
    groupBy: string,
    splitBy: string
  ): Promise<string> {
    if (!this._datasets.some((dataset) => dataset.ds === ds))
      throw new Error(`Attempt to render unloaded dataset ${ds}`);

    // TODO: kinda sus lol
    let id: string;
    return this.exec((proc, opid) => {
      proc.stdin.write(`${opid} render ${ds} ${gene} ${groupBy} ${splitBy}\n`);
      id = opid;
    }).then(() => id!);
  }

  public async unload(ds: string): Promise<void> {
    if (!this._datasets.some((dataset) => dataset.ds === ds))
      throw new Error(`Attempt to unload unloaded dataset ${ds}`);

    this._datasets.splice(
      this._datasets.findIndex((dataset) => dataset.ds === ds),
      1
    );

    return this.exec((proc, opid) =>
      proc.stdin.write(`${opid} unload ${ds}\n`)
    );
  }

  public kill(): void {
    this.proc.kill();
  }
}
