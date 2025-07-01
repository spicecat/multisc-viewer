"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.viteNodeApp = void 0;
const common_1 = require("@nestjs/common");
const core_1 = require("@nestjs/core");
const cookie_parser_1 = __importDefault(require("cookie-parser"));
const dotenv_1 = require("dotenv");
const app_module_1 = require("./app.module");
const template_engine_js_1 = require("./client/template-engine.js");
const error_page_filter_1 = require("./utils/filters/error-page.filter");
const redirect_filter_1 = require("./utils/filters/redirect.filter");
const routing_interceptor_1 = require("./utils/interceptors/routing.interceptor");
(0, dotenv_1.config)();
async function bootstrap() {
    const app = await core_1.NestFactory.create(app_module_1.AppModule);
    app.engine("svelte", template_engine_js_1.svelte);
    app.setViewEngine("svelte");
    app.setBaseViewsDir("src/client/routes");
    app
        .use((0, cookie_parser_1.default)())
        .useGlobalPipes(new common_1.ValidationPipe({
        transform: true,
        transformOptions: { enableImplicitConversion: true },
    }))
        .useGlobalFilters(new error_page_filter_1.ErrorPageFilter(app.get(core_1.HttpAdapterHost).httpAdapter), new redirect_filter_1.RedirectFilter())
        .useGlobalInterceptors(new routing_interceptor_1.RoutingInterceptor());
    if (process.env.NODE_ENV !== "development") {
        await app.listen(process.env.PORT || 5000);
    }
    return app;
}
exports.viteNodeApp = bootstrap();
//# sourceMappingURL=main.js.map