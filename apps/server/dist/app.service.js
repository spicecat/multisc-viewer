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
var AppService_1;
Object.defineProperty(exports, "__esModule", { value: true });
exports.AppService = void 0;
const common_1 = require("@nestjs/common");
const fs_1 = require("fs");
const publications_1 = require("./config/publications");
const daemon_service_1 = require("./daemon.service");
let AppService = class AppService {
    static { AppService_1 = this; }
    daemon;
    static META = JSON.parse((0, fs_1.readFileSync)('datasets/meta.json').toString());
    cache;
    expirations;
    activity;
    constructor(daemon) {
        this.daemon = daemon;
        this.cache = new Map();
        this.expirations = new Map();
        this.activity = new Map();
    }
    getPublications() {
        return publications_1.publications.map((publication) => ({
            ...publication,
            datasets: this.getDatasets() // TODO: add datasets to publication
        }));
    }
    getDatasets() {
        return (0, fs_1.readdirSync)('datasets')
            .filter((ds) => (0, fs_1.existsSync)(`datasets/${ds}/data.rds`) && (0, fs_1.existsSync)(`datasets/${ds}/genes.json`))
            .map((name) => AppService_1.META.find((ds) => ds.name === name))
            .filter((ds) => !!ds);
    }
    getGenes(datasets) {
        const geneSets = datasets.map((ds) => JSON.parse((0, fs_1.readFileSync)(`datasets/${ds}/genes.json`).toString()));
        return geneSets.slice(1).reduce((curr, next) => curr.filter((gene) => next.includes(gene)), geneSets[0]);
    }
    async render(dataset, gene, groupBy, splitBy, token) {
        const key = this._cacheKey(dataset, gene, groupBy, splitBy);
        if (token)
            this._refresh(token);
        if (this.cache.has(key)) {
            if (this.expirations.has(key))
                this.expirations.get(key).refresh();
            return this.cache.get(key);
        }
        else {
            const result = this.daemon.render(dataset, gene, groupBy, splitBy).then((opid) => {
                try {
                    const clustering = 'data:image/png;base64,' + (0, fs_1.readFileSync)(`datasets/${dataset}/${opid}_umap.png`).toString('base64');
                    const violin = 'data:image/png;base64,' + (0, fs_1.readFileSync)(`datasets/${dataset}/${opid}_vln.png`).toString('base64');
                    const feature = 'data:image/png;base64,' + (0, fs_1.readFileSync)(`datasets/${dataset}/${opid}_feat.png`).toString('base64');
                    (0, fs_1.rmSync)(`datasets/${dataset}/${opid}_umap.png`);
                    (0, fs_1.rmSync)(`datasets/${dataset}/${opid}_vln.png`);
                    (0, fs_1.rmSync)(`datasets/${dataset}/${opid}_feat.png`);
                    return { clustering, violin, feature };
                }
                catch (e) {
                    console.error(e);
                    throw new Error(`Unable to find plot of '${dataset}'`);
                }
            });
            result
                .then(() => this.expirations.set(key, setTimeout(() => {
                this.cache.delete(key);
                this.expirations.delete(key);
            }, 30 * 60 * 1000)))
                .then(() => {
                if (token)
                    this._refresh(token);
            });
            this.cache.set(key, result);
            return result;
        }
    }
    async preload(token, datasets) {
        if (token !== null) {
            this._refresh(token);
        }
        return this.daemon.batch(token, datasets).then(() => {
            if (token)
                this._refresh(token);
        });
    }
    _cacheKey(dataset, gene, groupBy, splitBy) {
        return `${dataset}:${gene}:${groupBy}:${splitBy}`;
    }
    _refresh(token) {
        if (this.activity.has(token))
            this.activity.get(token).refresh();
        else {
            this.activity.set(token, setTimeout(() => {
                this.activity.delete(token);
                this.daemon.unload(token);
            }, 60_000));
        }
    }
};
exports.AppService = AppService;
exports.AppService = AppService = AppService_1 = __decorate([
    (0, common_1.Injectable)(),
    __metadata("design:paramtypes", [daemon_service_1.DaemonService])
], AppService);
//# sourceMappingURL=app.service.js.map