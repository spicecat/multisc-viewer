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
var __param = (this && this.__param) || function (paramIndex, decorator) {
    return function (target, key) { decorator(target, key, paramIndex); }
};
var _a, _b;
Object.defineProperty(exports, "__esModule", { value: true });
exports.AppController = void 0;
const common_1 = require("@nestjs/common");
const crypto_1 = require("crypto");
const app_service_1 = require("./app.service");
const page_decorator_1 = require("./utils/decorators/page.decorator");
let AppController = class AppController {
    service;
    constructor(service) {
        this.service = service;
    }
    getPublications() {
        return this.service.getPublications();
    }
    getGenes(datasets) {
        return this.service.getGenes(datasets.split(","));
    }
    getDatasets() {
        return this.service.getDatasets();
    }
    async getPlot(dataset, gene, groupBy, splitBy, token = null) {
        return this.service.render(dataset, gene, groupBy, splitBy, token);
    }
    async preload(token = null, datasets) {
        return this.service.preload(token, datasets.split(","));
    }
    index() {
        const datasets = this.service.getDatasets();
        return { datasets, token: (0, crypto_1.randomBytes)(32).toString("hex") };
    }
    publication() {
        return {
            publications: this.service.getPublications(),
        };
    }
    publicationId(publicationId) {
        const publication = this.service
            .getPublications()
            .find((publication) => publication.publicationId === publicationId);
        if (!publication)
            throw new common_1.BadRequestException("Unknown publication requested");
        return { publication };
    }
    comparison(datasets = "", gene = "", groupBy = "", splitBy = "", token = null) {
        const knownDatasets = Object.fromEntries(this.getDatasets().map((ds) => [ds.name, ds]));
        if (!datasets?.split(",").every((ds) => ds in knownDatasets))
            throw new common_1.BadRequestException("Unknown dataset requested");
        const genes = this.getGenes(datasets);
        if (!gene) {
            const defaultGene = genes[0];
            if (!defaultGene)
                throw new common_1.BadRequestException("No common genes between datasets");
            gene = defaultGene;
        }
        if (!genes.includes(gene))
            throw new common_1.BadRequestException("Selected gene is not in all datasets");
        if (!groupBy && !splitBy) {
            groupBy = "Genotype";
            splitBy = "CellType";
        }
        else if (!groupBy)
            groupBy = splitBy === "Genotype" ? "CellType" : "Genotype";
        else if (!splitBy)
            splitBy = groupBy === "Genotype" ? "CellType" : "Genotype";
        this.service.preload(token, datasets.split(","));
        return {
            order: datasets
                .split(",")
                .sort((a, b) => knownDatasets[a].size - knownDatasets[b].size),
            genes,
            gene,
        };
    }
};
exports.AppController = AppController;
__decorate([
    (0, common_1.Get)("/publications"),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", []),
    __metadata("design:returntype", Array)
], AppController.prototype, "getPublications", null);
__decorate([
    (0, common_1.Get)("/genes"),
    __param(0, (0, common_1.Query)("datasets")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String]),
    __metadata("design:returntype", Array)
], AppController.prototype, "getGenes", null);
__decorate([
    (0, common_1.Get)("/datasets"),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", []),
    __metadata("design:returntype", Array)
], AppController.prototype, "getDatasets", null);
__decorate([
    (0, common_1.Get)("/plot"),
    __param(0, (0, common_1.Query)("dataset")),
    __param(1, (0, common_1.Query)("gene")),
    __param(2, (0, common_1.Query)("groupBy")),
    __param(3, (0, common_1.Query)("splitBy")),
    __param(4, (0, common_1.Query)("token")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String, String, String, String, Object]),
    __metadata("design:returntype", Promise)
], AppController.prototype, "getPlot", null);
__decorate([
    (0, common_1.Post)("/preload"),
    __param(0, (0, common_1.Query)("token")),
    __param(1, (0, common_1.Query)("datasets")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [Object, String]),
    __metadata("design:returntype", Promise)
], AppController.prototype, "preload", null);
__decorate([
    (0, page_decorator_1.Page)(),
    (0, common_1.Get)("/"),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", []),
    __metadata("design:returntype", typeof (_a = typeof IndexProps !== "undefined" && IndexProps) === "function" ? _a : Object)
], AppController.prototype, "index", null);
__decorate([
    (0, page_decorator_1.Page)(),
    (0, common_1.Get)("/publication"),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", []),
    __metadata("design:returntype", void 0)
], AppController.prototype, "publication", null);
__decorate([
    (0, page_decorator_1.Page)(),
    (0, common_1.Get)("/publication/:publicationId"),
    __param(0, (0, common_1.Param)("publicationId")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String]),
    __metadata("design:returntype", void 0)
], AppController.prototype, "publicationId", null);
__decorate([
    (0, page_decorator_1.Page)(),
    (0, common_1.Get)("/compare"),
    __param(0, (0, common_1.Query)("datasets")),
    __param(1, (0, common_1.Query)("gene")),
    __param(2, (0, common_1.Query)("groupBy")),
    __param(3, (0, common_1.Query)("splitBy")),
    __param(4, (0, common_1.Query)("token")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String, String, String, String, Object]),
    __metadata("design:returntype", typeof (_b = typeof CompareProps !== "undefined" && CompareProps) === "function" ? _b : Object)
], AppController.prototype, "comparison", null);
exports.AppController = AppController = __decorate([
    (0, common_1.Controller)(),
    __metadata("design:paramtypes", [app_service_1.AppService])
], AppController);
//# sourceMappingURL=app.controller.js.map