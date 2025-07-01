"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
var __metadata = (this && this.__metadata) || function (k, v) {
    if (typeof Reflect === "object" && typeof Reflect.metadata === "function") return Reflect.metadata(k, v);
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.DaemonService = void 0;
const common_1 = require("@nestjs/common");
const app_service_1 = require("./app.service");
const Daemon_1 = require("./utils/Daemon");
const utils_1 = require("./utils/utils");
let DaemonService = class DaemonService {
    procs;
    constructor() {
        this.procs = new Array(10).fill(null).map(() => new Daemon_1.Daemon());
    }
    async batch(loader, datasets) {
        return Promise.all(datasets.map((ds) => this.load(loader, ds))).then();
    }
    async load(loader, ds) {
        const proc = this.procs.find((proc) => proc.datasets.some((dataset) => dataset.ds === ds));
        if (proc) {
            if (loader)
                proc.datasets.find((dataset) => dataset.ds === ds).loaders.add(loader);
            // TODO: not precise (could be an operation later than the actual end of the load of this dataset)
            // consider changing to Promise.resolve if all operations dependent on this case are queued operations
            return proc.idle;
        }
        else {
            if (loader) {
                const availableProc = this.procs.find((proc) => !proc.datasets.some((set) => set.loaders.has(loader)));
                if (availableProc)
                    return availableProc.load(loader, ds);
                else {
                    const newProc = this._spawn();
                    return newProc.load(loader, ds);
                }
            }
            else {
                const bestProc = this.procs.reduce((best, proc) => proc.datasets.length < best.datasets.length ? proc : best, this.procs[0]);
                return bestProc.load(loader, ds).then(() => {
                    setTimeout(() => {
                        if (bestProc.datasets.find((dataset) => dataset.ds === ds)?.loaders
                            .size === 0) {
                            if (bestProc.datasets.length === 1 && this.procs.length > 10) {
                                bestProc.kill();
                                this.procs.splice(this.procs.indexOf(bestProc), 1);
                            }
                            else
                                bestProc.unload(ds);
                        }
                    }, (0, utils_1.approximateRenderTime)(app_service_1.AppService.META.find((dataset) => dataset.name === ds).size) * 30);
                });
            }
        }
    }
    async unload(loader) {
        const promises = [];
        this.procs.forEach((proc) => {
            proc.datasets.forEach((ds) => {
                if (ds.loaders.has(loader)) {
                    if (ds.loaders.size === 1) {
                        if (proc.datasets.length === 1 && this.procs.length > 10) {
                            proc.kill();
                            this.procs.splice(this.procs.indexOf(proc), 1);
                        }
                        else {
                            promises.push(proc.unload(ds.ds));
                            proc.datasets.splice(proc.datasets.indexOf(ds), 1);
                        }
                    }
                    else
                        ds.loaders.delete(loader);
                }
            });
        });
        return Promise.all(promises).then();
    }
    async render(ds, gene, groupBy, splitBy) {
        const proc = this.procs.find((proc) => proc.datasets.some((dataset) => dataset.ds === ds));
        if (proc)
            return proc.render(ds, gene, groupBy, splitBy);
        else
            return this.load(null, ds).then(() => this.render(ds, gene, groupBy, splitBy));
    }
    _spawn() {
        const proc = new Daemon_1.Daemon();
        this.procs.push(proc);
        return proc;
    }
};
exports.DaemonService = DaemonService;
exports.DaemonService = DaemonService = __decorate([
    (0, common_1.Injectable)(),
    __metadata("design:paramtypes", [])
], DaemonService);
//# sourceMappingURL=daemon.service.js.map