"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Daemon = void 0;
const node_child_process_1 = require("node:child_process");
const node_crypto_1 = require("node:crypto");
class Daemon {
    proc;
    _acks;
    _datasets;
    _idle;
    constructor() {
        this.proc = (0, node_child_process_1.spawn)("Rscript", ["daemon.r"]);
        this.proc.stdin.on("error", (error) => {
            if (error.code === "EPIPE")
                console.error("EPIPE error: child process pipe is closed.", error);
            else
                console.error("Error on stdin:", error);
        });
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
        this.proc.stderr.on("data", (chunk) => console.error(chunk.toString()));
    }
    get idle() {
        return this._idle;
    }
    get datasets() {
        return this._datasets;
    }
    exec(executor) {
        const opid = (0, node_crypto_1.randomBytes)(16).toString("hex");
        const promise = new Promise((resolve, reject) => {
            this._idle.then(() => {
                executor(this.proc, opid);
                this._acks.push((id) => {
                    if (id === opid) {
                        resolve();
                        return true;
                    }
                    else
                        return false;
                });
            });
        });
        this._idle = promise;
        return promise;
    }
    async load(loader, ds) {
        const dataset = {
            ds,
            loaders: loader !== null ? new Set([loader]) : new Set(),
        };
        this._datasets.push(dataset);
        return this.exec((proc, opid) => proc.stdin.write(`${opid} load ${ds}\n`));
    }
    async render(ds, gene, groupBy, splitBy) {
        if (!this._datasets.some((dataset) => dataset.ds === ds))
            throw new Error(`Attempt to render unloaded dataset ${ds}`);
        // TODO: kinda sus lol
        let id;
        return this.exec((proc, opid) => {
            proc.stdin.write(`${opid} render ${ds} ${gene} ${groupBy} ${splitBy}\n`);
            id = opid;
        }).then(() => id);
    }
    async unload(ds) {
        if (!this._datasets.some((dataset) => dataset.ds === ds))
            throw new Error(`Attempt to unload unloaded dataset ${ds}`);
        this._datasets.splice(this._datasets.findIndex((dataset) => dataset.ds === ds), 1);
        return this.exec((proc, opid) => proc.stdin.write(`${opid} unload ${ds}\n`));
    }
    kill() {
        this.proc.kill();
    }
}
exports.Daemon = Daemon;
//# sourceMappingURL=Daemon.js.map