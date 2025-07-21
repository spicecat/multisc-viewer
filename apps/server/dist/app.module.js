"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.AppModule = void 0;
const common_1 = require("@nestjs/common");
const config_1 = require("@nestjs/config");
const serve_static_1 = require("@nestjs/serve-static");
const app_controller_1 = require("./app.controller");
const configuration_1 = __importDefault(require("./config/configuration"));
const daemon_service_1 = require("./daemon.service");
const data_service_1 = require("./data.service");
const plot_service_1 = require("./plot.service");
let AppModule = class AppModule {
};
exports.AppModule = AppModule;
exports.AppModule = AppModule = __decorate([
    (0, common_1.Module)({
        imports: [
            config_1.ConfigModule.forRoot({
                load: [configuration_1.default],
            }),
            serve_static_1.ServeStaticModule.forRoot({
                rootPath: "../../viewer/build/client",
                serveRoot: "/",
            }),
        ],
        controllers: [app_controller_1.AppController],
        providers: [data_service_1.DataService, plot_service_1.PlotService, daemon_service_1.DaemonService],
    })
], AppModule);
//# sourceMappingURL=app.module.js.map