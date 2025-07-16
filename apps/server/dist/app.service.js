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
exports.AppService = void 0;
const common_1 = require("@nestjs/common");
const config_1 = require("@nestjs/config");
const node_cache_1 = __importDefault(require("node-cache"));
const node_fs_1 = require("node:fs");
const daemon_service_1 = require("./daemon.service");
let AppService = class AppService {
    configService;
    daemon;
    plotCache;
    datasetsBasePath;
    publications;
    datasets;
    constructor(configService, daemon) {
        this.configService = configService;
        this.daemon = daemon;
        this.plotCache = new node_cache_1.default({
            stdTTL: this.configService.get("plot.ttl"),
        });
        const publicationsMeta = `${this.configService.get("publications.basePath")}/${this.configService.get("publications.metaFile")}`;
        this.publications = new Map(JSON.parse((0, node_fs_1.readFileSync)(publicationsMeta).toString()).map((pub) => [pub.PMID, pub]));
        this.datasetsBasePath =
            this.configService.get("datasets.basePath");
        const datasetsMeta = `${this.datasetsBasePath}/${this.configService.get("datasets.metaFile")}`;
        this.datasets = new Map(JSON.parse((0, node_fs_1.readFileSync)(datasetsMeta).toString())
            .filter(({ name }) => this.configService
            .get("datasets.requiredFiles")
            .every((file) => (0, node_fs_1.existsSync)(`${this.datasetsBasePath}/${name}/${file}`)))
            .map((ds) => [ds.name, ds]));
    }
    getGenes(datasets) {
        return Array.from(new Set(datasets
            .map((ds) => JSON.parse((0, node_fs_1.readFileSync)(`${this.datasetsBasePath}/${ds}/genes.json`).toString()))
            .flat())).toSorted();
    }
    async render(dataset, gene, groupBy, splitBy) {
        const key = this._cacheKey(dataset, gene, groupBy, splitBy);
        const cachedResult = this.plotCache.get(key);
        if (cachedResult)
            return cachedResult;
        const result = this.daemon
            .render(dataset, gene, groupBy, splitBy)
            .then((opid) => {
            try {
                const clusteringPath = `${this.datasetsBasePath}/${dataset}/${opid}_umap.png`;
                const violinPath = `${this.datasetsBasePath}/${dataset}/${opid}_vln.png`;
                const featurePath = `${this.datasetsBasePath}/${dataset}/${opid}_feat.png`;
                const clustering = "data:image/png;base64," +
                    (0, node_fs_1.readFileSync)(clusteringPath).toString("base64");
                const violin = "data:image/png;base64," +
                    (0, node_fs_1.readFileSync)(violinPath).toString("base64");
                const feature = "data:image/png;base64," +
                    (0, node_fs_1.readFileSync)(featurePath).toString("base64");
                (0, node_fs_1.rmSync)(clusteringPath);
                (0, node_fs_1.rmSync)(violinPath);
                (0, node_fs_1.rmSync)(featurePath);
                return { clustering, violin, feature };
            }
            catch (e) {
                console.error(e);
                this.plotCache.del(key);
                throw new Error(`Unable to find plot for '${dataset}'`);
            }
        });
        this.plotCache.set(key, result);
        return result;
    }
    async load(datasets) {
        return this.daemon.load(datasets);
    }
    _cacheKey(dataset, gene, groupBy, splitBy) {
        return `${dataset}:${gene}:${groupBy}:${splitBy}`;
    }
};
exports.AppService = AppService;
exports.AppService = AppService = __decorate([
    (0, common_1.Injectable)(),
    __metadata("design:paramtypes", [config_1.ConfigService,
        daemon_service_1.DaemonService])
], AppService);
//# sourceMappingURL=app.service.js.map