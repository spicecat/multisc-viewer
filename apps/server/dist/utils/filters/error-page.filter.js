"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ErrorPageFilter = exports.ErrorPage = void 0;
const common_1 = require("@nestjs/common");
const core_1 = require("@nestjs/core");
class ErrorPage {
    props;
    constructor(props) {
        this.props = props;
    }
}
exports.ErrorPage = ErrorPage;
let ErrorPageFilter = class ErrorPageFilter extends core_1.BaseExceptionFilter {
    catch(exception, host) {
        const ctx = host.switchToHttp();
        const response = ctx.getResponse();
        if (exception instanceof common_1.HttpException) {
            const err = exception.getResponse();
            if (err instanceof ErrorPage) {
                response.status(exception.getStatus()).render("error", err.props);
                return;
            }
        }
        super.catch(exception, host);
    }
};
exports.ErrorPageFilter = ErrorPageFilter;
exports.ErrorPageFilter = ErrorPageFilter = __decorate([
    (0, common_1.Catch)()
], ErrorPageFilter);
//# sourceMappingURL=error-page.filter.js.map