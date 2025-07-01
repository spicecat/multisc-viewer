"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.RedirectFilter = exports.Redirect = void 0;
const common_1 = require("@nestjs/common");
class Redirect {
    location;
    status;
    constructor(location, status = common_1.HttpStatus.SEE_OTHER) {
        this.location = location;
        this.status = status;
    }
}
exports.Redirect = Redirect;
let RedirectFilter = class RedirectFilter {
    catch(exception, host) {
        const ctx = host.switchToHttp();
        const response = ctx.getResponse();
        response
            .status(exception.status)
            .setHeader("Location", exception.location)
            .end();
    }
};
exports.RedirectFilter = RedirectFilter;
exports.RedirectFilter = RedirectFilter = __decorate([
    (0, common_1.Catch)(Redirect)
], RedirectFilter);
//# sourceMappingURL=redirect.filter.js.map