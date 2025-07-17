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
Object.defineProperty(exports, "__esModule", { value: true });
exports.AppController = void 0;
const common_1 = require("@nestjs/common");
const data_service_1 = require("./data.service");
const plot_service_1 = require("./plot.service");
const daemon_service_1 = require("./daemon.service");
let AppController = class AppController {
    dataService;
    plotService;
    daemonService;
    constructor(dataService, plotService, daemonService) {
        this.dataService = dataService;
        this.plotService = plotService;
        this.daemonService = daemonService;
    }
    getPublications() {
        return [...this.dataService.publications.values()];
    }
    getDatasets() {
        return [...this.dataService.datasets.values()];
    }
    getGenes(datasets = "") {
        return this.dataService.getGenes(datasets.split(","));
    }
    async load(datasets = "") {
        return this.daemonService.load(datasets.split(","));
    }
    async getPlot(dataset, gene, groupBy, splitBy) {
        return this.plotService.render(dataset, gene, groupBy, splitBy);
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
    (0, common_1.Get)("/datasets"),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", []),
    __metadata("design:returntype", Array)
], AppController.prototype, "getDatasets", null);
__decorate([
    (0, common_1.Get)("/genes"),
    __param(0, (0, common_1.Query)("datasets")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String]),
    __metadata("design:returntype", Array)
], AppController.prototype, "getGenes", null);
__decorate([
    (0, common_1.Post)("/load"),
    __param(0, (0, common_1.Body)("datasets")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String]),
    __metadata("design:returntype", Promise)
], AppController.prototype, "load", null);
__decorate([
    (0, common_1.Get)("/plot"),
    __param(0, (0, common_1.Query)("dataset")),
    __param(1, (0, common_1.Query)("gene")),
    __param(2, (0, common_1.Query)("groupBy")),
    __param(3, (0, common_1.Query)("splitBy")),
    __metadata("design:type", Function),
    __metadata("design:paramtypes", [String, String, String, String]),
    __metadata("design:returntype", Promise)
], AppController.prototype, "getPlot", null);
exports.AppController = AppController = __decorate([
    (0, common_1.Controller)(),
    __metadata("design:paramtypes", [data_service_1.DataService,
        plot_service_1.PlotService,
        daemon_service_1.DaemonService])
], AppController);
//# sourceMappingURL=app.controller.js.map