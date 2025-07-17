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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.DaemonService = void 0;
const common_1 = require("@nestjs/common");
const config_1 = require("@nestjs/config");
const axios_1 = __importDefault(require("axios"));
const node_cache_1 = __importDefault(require("node-cache"));
const data_service_1 = require("./data.service");
let DaemonService = class DaemonService {
    configService;
    dataService;
    DAEMON_SERVER;
    datasetCache;
    procs = new Map();
    constructor(configService, dataService) {
        this.configService = configService;
        this.dataService = dataService;
        this.DAEMON_SERVER = this.configService.get("daemon.server");
        this.datasetCache = new node_cache_1.default({
            stdTTL: this.configService.get("daemon.ttl"),
        });
        this.datasetCache.on("expired", (ds, daemon) => this._unload(daemon, ds));
    }
    async onModuleInit() {
        console.log("Initializing DaemonService...");
        await new Promise((resolve) => setTimeout(resolve, 2000));
        await Promise.all(this.configService.get("daemon.ports").map(async (port) => {
            try {
                const response = await axios_1.default.get(`${this.DAEMON_SERVER}:${port}/status`);
                if (response.data)
                    this.procs.set(port, {
                        port,
                        load: 0,
                        datasets: new Set(),
                    });
            }
            catch (error) {
                console.error(`Daemon on port ${port} not responding`);
            }
        }));
    }
    async _unload(daemon, ds) {
        try {
            await axios_1.default.post(`${this.DAEMON_SERVER}:${daemon.port}/unload`, { ds });
            daemon.datasets.delete(ds);
            const datasetInfo = this.dataService.datasets.get(ds);
            if (datasetInfo) {
                daemon.load -= datasetInfo.size;
            }
        }
        catch (error) {
            console.error(`Failed to unload dataset ${ds} from port ${daemon.port}:`, error.message);
        }
    }
    async load(datasets) {
        if (this.procs.size === 0)
            throw new Error("No daemons available to load dataset.");
        if (!Array.isArray(datasets))
            datasets = [datasets];
        await Promise.all(datasets.map(async (ds) => {
            if (this.datasetCache.has(ds))
                this.datasetCache.ttl(ds);
            else {
                const datasetInfo = this.dataService.datasets.get(ds);
                if (!datasetInfo) {
                    console.error(`Dataset ${ds} not found.`);
                    return;
                }
                const daemon = [...this.procs.values()].reduce((a, b) => a.load < b.load ? a : b);
                daemon.load += datasetInfo.size;
                try {
                    await axios_1.default.post(`${this.DAEMON_SERVER}:${daemon.port}/load`, {
                        ds,
                    });
                    daemon.datasets.add(ds);
                    this.datasetCache.set(ds, daemon);
                }
                catch (error) {
                    daemon.load -= datasetInfo.size;
                    console.error(`Failed to load dataset ${ds} on port ${daemon.port}:`, error.message);
                }
            }
        }));
        console.log("Datasets in cache:", this.datasetCache.keys());
    }
    async render(ds, gene, groupBy, splitBy) {
        await this.load(ds);
        const daemon = this.datasetCache.get(ds);
        if (!daemon)
            throw new Error(`Dataset ${ds} is not loaded`);
        try {
            const response = await axios_1.default.post(`${this.DAEMON_SERVER}:${daemon.port}/render`, {
                ds,
                gene,
                groupBy,
                splitBy,
            });
            return response.data.opid;
        }
        catch (error) {
            console.error(`Failed to render dataset ${ds} on port ${daemon.port}:`, error.message);
            throw new Error(`Failed to render plot for dataset ${ds}`);
        }
    }
};
exports.DaemonService = DaemonService;
exports.DaemonService = DaemonService = __decorate([
    (0, common_1.Injectable)(),
    __metadata("design:paramtypes", [config_1.ConfigService,
        data_service_1.DataService])
], DaemonService);
//# sourceMappingURL=daemon.service.js.map